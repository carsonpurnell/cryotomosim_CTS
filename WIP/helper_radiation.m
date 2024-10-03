function [out,rad,rad2] = helper_radiation(vol,solv,dose,rad)

radscale = .01; %arbitrary scalar for parameter values to map correctly to map intensity

% gaussian component
radg = (rand(size(vol))-rand(size(vol)))*rad;
% smoothing component - run on vol separately? no, would generate weirdness
rads = imgaussfilt3(vol,rad);
% separate filling/warping component? or suitably part of others?
% solv needed as a separate component at all?

rad = 0; rad2=0;
out = vol+solv;

blurmap = imgaussfilt( max(tilt,[],'all')-tilt ); %2d blur each angle outside loop for speed
blurmean = imgaussfilt(tilt,0.5);
for i=1
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

end