function [out,loss,noise] = helper_radiation(vol,pix,dose,rad,opt)

arguments
    vol
    pix
    dose
    rad
    opt.byslice = 1
end

rads = .0100*rad*dose*pix^1; %arbitrary scalar for parameter values to map correctly to map intensity

% quantification from https://journals.iucr.org/s/issues/2011/03/00/xh5022/xh5022.pdf
H = 8e2; % conversion constant, don't know true value from KeV to Gy
% signal = ideal * exp(-log(2)*dose*A/H);

% gaussian component
%noise = rads*10*randn(size(vol));%(rand(size(vol))-rand(size(vol)))*rads*10;
% smoothing component - run on vol separately? no, would generate weirdness
%f = 15/pix*rads; f = ceil(f)*2+1; %default 2*ceil(sigma*2)+1;
loss = zeros(size(vol)); noise = loss; %nweighted = loss;
if opt.byslice
    for i=1:size(vol,3)
        sc = (i+0)/size(vol,3);
        r = rads*sc; % 0 placeholder for pre-exposures?
        f = ceil(10/pix+r*1)*2+1;
        noise(:,:,i) = r*50*randn(size(vol,[1,2]))*sqrt(pix);
        decay = 1-(1-exp(-dose*sc*pix/H))/2;
        loss(:,:,i) = imgaussfilt3(vol(:,:,i),r,'FilterSize',f)*decay;
    end
else
    f = ceil(10/pix+rads*1)*2+1;
    noise = rads*50*randn(size(vol));
    loss = imgaussfilt3(vol,rads,'FilterSize',f)*exp(-dose*pix/H);
end
% loss: smoothing component of radiation damage
% noise: gaussian component

out = loss+noise;

%{
blurmap = imgaussfilt( max(tilt,[],'all')-tilt ); %2d blur each angle outside loop for speed
blurmean = imgaussfilt(tilt,0.5);
for i=1:size(vol,3)
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
end
%}

end