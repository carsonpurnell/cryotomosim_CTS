function detect = helper_electrondetect(tilt,param)
%detect = helper_electrondetect(tilt,param)
%
%
%
%lazy temporary help
%tilt is a tilt series (ideally created by helper_ctf)
%param is a param struct created by cts_param, it's important
arguments
    tilt
    param struct
end
if param.dose<=0, detect=tilt; return; end %if dose 0, skip detection and return perfect detection/original
tiltangs = param.tilt; %unfortunately similar name to tilt 

%arb = 4;%/param.pix^2; %arbitrary scaling factor to make contrast look normal
%what are the new good values? is this scale working well?

DQE = .84*4; % gatan camera lists 84% maximum detection, so that'll work for now
%doesn't change with tilts, not sure how to implement fourier space falloff

switch param.tiltscheme %organize tilt ordering, sort according to split if not symmetric
    case 'symmetric'
        [~,ix] = sort(abs(tiltangs)); 
    otherwise %find and sort the data from the split between tilt directions
        mdist = max(abs(tiltangs-(param.tiltscheme))); %max dist from start angle where phase will switch
        metric = abs(tiltangs-param.tiltscheme)+mdist.*(tiltangs<param.tiltscheme); %calculate sorting metric
        [~,ix] = sort(metric);
end
tiltangs = tiltangs(ix); tilt = tilt(:,:,ix); %sort tilt angles and tilts
ixr(ix) = 1:numel(ix); %generate reverse sorting index

%dose weighting/distribution
dose = param.dose;%.*param.pix^2; %convert dose in e/A^2 to e/pixel
if numel(param.dose)==1
    dose = dose/size(tilt,3); %single dose distributed across tilts evenly
else
    dose = dose(ix); %for weighting just sort the weights to the tilts
end

thick = param.size(3)*param.pix; %compute thickness from depth of original model
IMFP = 3800; %inelastic mean free path, average distance before inelastic electron scatter (for water/ice)
%IMFP estimated to be 350nm for water ice, is probabaly somewhat different for vitreous (higher)
%electronpath = thick*(1+abs(tand(tiltangs))); %compute the path length of electrons through ice
electronpath = thick*cosd(tiltangs).^-1; %corrected trig, very slightly better appearance
thickscatter = exp(-(electronpath*param.scatter)/IMFP); %compute electrons not inelastically/lossly scattered
%change IMFP to instead be per pixel, so more electrons are lost at high density AND thickness?

radscale = .09*param.raddamage;%/param.pix^2; %damage scaling calculation to revert scaling by pixel size

dw = thickscatter.*dose*DQE;
accum = 0; %initialize accumulated dose of irradiation to 0
detect = tilt.*0; %pre-initialize output array for speed during the loop
for i=1:size(tilt,3)
    irad = tilt(:,:,i)+randn(size(tilt(:,:,i)))*accum*radscale; %radiation as 0-center noise
    %need to do some procedure to mask/weight the noise near density rather than globally
    
    accum = accum+dw(i); %add to accumulated dose delivered
    detect(:,:,i) = poissrnd(irad*dw(i),size(irad));
end

detect = detect(:,:,ixr); %reverse the sort so the output tiltseries is a continuous rotation

end