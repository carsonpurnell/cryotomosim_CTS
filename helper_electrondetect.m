function [detect,rad] = helper_electrondetect(tilt,param)
%[detect,rad] = helper_electrondetect(tilt,param)
%simulates electron detection by camera of a CTF-convolved tiltseries. simulates electon transmission and
%interaction with poisson sampling, after estimating electrons lost by inelastic scattering.
%
%tilt - tiltseries to simulate detection through
%param - parameter structure generated by cts_param, see that function for arguments and help
%
%detect - output of simulated electron detection by counting mode camera
arguments
    tilt
    param struct
end
rad = 0;
if param.dose<=0, detect=tilt; return; end %if dose 0, skip detection and return perfect detection/original
tiltangs = param.tilt; %unfortunately similar name to tilt 

DQE = .84*.6; % gatan camera lists 84% maximum detection, so that'll work for now
%4 is arbitrary scalar to make contrast look 'normal' with the CTF modulation
%DQE should not be angle dependent, but maybe easier to  implement if merged with CTF?

switch param.tiltscheme %organize tilt ordering, sort according to split if not symmetric
    case 'symmetric'
        [~,ix] = sort(abs(tiltangs)); %straight sort, can't do grouped symmetric
    otherwise %find and sort the data from the split between tilt directions
        mdist = max(abs(tiltangs-(param.tiltscheme))); %max dist from start angle where phase will switch
        metric = abs(tiltangs-param.tiltscheme)+mdist.*(tiltangs<param.tiltscheme); %calculate sorting metric
        [~,ix] = sort(metric);
end
tiltangs = tiltangs(ix); tilt = tilt(:,:,ix); %sort tilt angles and tilts
ixr(ix) = 1:numel(ix); %generate reverse sorting index

%dose weighting/distribution
dose = param.dose; %.*param.pix^2; %convert dose in e/A^2 to e/pixel - currently deprecated
if numel(param.dose)==1
    dose = dose/size(tilt,3); %single dose distributed across tilts evenly
else
    dose = dose(ix); %for weighting just sort the weights to the tilts
end

%pix = 10;
box = param.size*param.pix; %calc from input data
if strcmp(param.tiltax,'X'), box=box([2,1,3]); end
%angles = -60:5:60;
%if strcmp(param.tiltax,'Y'), axis=1; else axis=2; end
axis = 1;
axspec = zeros(1,3); axspec(1+rem(axis,2)) = 1;
%axis = 1;
sc = [2.0,1.2]; % [perlin magnitude, surface centering
%rng(7)
surfaces = helper_surf(box,param.pix,tiltangs,axis,sc);
%histogram(surfaces{1}(:,3)); hold on
%histogram(surfaces{2}(:,3));
[newpath,gridt] = helper_thickfromsurf(surfaces,box,param.pix,tiltangs,axspec); 
%[a,b] = bounds(newpath,'all')
%histogram(newpath)

thick = param.size(3)*param.pix; %compute thickness from depth of original model
IMFP = 3100; %inelastic mean free path, average distance before inelastic electron scatter (for water/ice)
%IMFP estimated to be 350nm for water ice, is probabaly somewhat different for vitreous (higher)
%3100A (at 300keV) from https://www.sciencedirect.com/science/article/pii/S0304399122000456#b17 
%electronpath = thick*(1+abs(tand(tiltangs))); %compute the path length of electrons through ice
electronpath = thick./cosd(tiltangs); %corrected trig, slightly better appearance

%make thickscatter a map? random variation in ice thickness (should be highly smooth)
thickscatter = exp(-(electronpath*param.scatter)/IMFP); %compute electrons not inelastically/lossly scattered
newscatter = exp(-(newpath*param.scatter)/IMFP);
%change IMFP to instead be per pixel, so more electrons are lost at high density AND thickness?
%scattering map inside the loop, using pixel intensities to scale the IMFP/path?
%problematic: single map of thickness- kind of wrong because thickness is directional, not location
%potential: 2 mats, thickness above/below center (inc. model thickness) - measure thickness at angles
%   need some sort of interpolation to match together at arbitrary angles, direct array math will be wrong

radscale = .01*param.raddamage;%/param.pix^2; %damage scaling calculation to revert scaling by pixel size

dw = dose*DQE; %correct distributed dose based on maximum DQE and inelastic scattering loss
accum = 0; %initialize accumulated dose of irradiation to 0
detect = tilt.*0; rad = tilt*0; %pre-initialize output arrays
blurmap = imgaussfilt( max(tilt,[],'all')-tilt ); %2d blur each angle outside loop for speed
blurmean = imgaussfilt(tilt,0.5);
%separate subfunction for doing radiation for modularity

for i=1:size(tilt,3)
    
    if param.raddamage>0 %block for raadiation-induced noise and blurring
        accum = accum+dw*thickscatter(i); %add to accumulated dose delivered, including first tilt
        %this radiation count is inappropriate. needs to increase at higher tilt, and ignore DQE/etc.
        %use raw dose number and adjust for angle? or precompute rad scalars outside loop?
        
        %need to use the pre-CTF tilt for the rad map to avoid CTF impacts
        radmap = rescale(blurmap(:,:,i),0,sqrt(param.pix))*1; %increase noise at proteins - what is good scale?
        %bidirectional general noise - general SNR reduction
        radgauss = randn(size(radmap))*(1/2)*1; % 0-centered unstructured noise field
        % divide by pixel size to reduce impact at lower resolutions?
        % smoothed independent noise field on top?
        radclose = rand(size(radmap)).*(blurmean(:,:,i)-tilt(:,:,i))*(2/1)*1; % contrast-reducing noise add
        % accum still flips into positive intensity total, need to plateu at the mean
        addrad = randn(size(radmap))*accum*radscale.*(radmap+1)/10*0; %scaled gaussian 0-center noise field
        %additive noise biased to low density - reduce contrast/signal differentiation
        addrad = addrad+blurmap(:,:,i).*abs(rand(size(radmap)))*(param.pix)*accum*radscale/1e2;
        
        sigma = sqrt(radscale*(accum)*0.2)*1; %might need to scale filter size with pixel size
        %smoothing noise - reduce resolution and contrast
        proj = imgaussfilt(tilt(:,:,i),sigma);%,'FilterSize',filt);
        
        irad = proj*1+tilt(:,:,i)*0+accum*radscale*(radgauss+radclose);
        rad(:,:,i) = proj*1+tilt(:,:,i)*0+accum*radscale*(radgauss+radclose); %store radiation maps for review
    else
        irad = tilt(:,:,i);
    end
    %size(irad)
    %size(newscatter(:,:,i))
    detect(:,:,i) = poissrnd(irad.*newscatter(:,:,i)*dw,size(irad)); 
    %sample electrons from scaled poisson distribution
end
rad = rad(:,:,ixr);
detect = detect(:,:,ixr); %reverse the sort so the output tiltseries is a continuous rotation
end