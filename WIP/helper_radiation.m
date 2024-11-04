function [out,loss,noise] = helper_radiation(vol,pix,dose,rad,opt)
% [out,loss,noise] = helper_radiation(vol,pix,dose,rad,opt)
% simulates local radiation damage from low-dose electron beam exposure (no burning/bubbling/warping)

arguments
    vol
    pix
    dose
    rad
    opt.byslice = 1
end
if rad==0; out=vol; loss=0; noise=0; return; end
rads = .02*rad*dose*pix^0; % arbitrary scalar for parameter values to map correctly to intensity

% quantification from https://journals.iucr.org/s/issues/2011/03/00/xh5022/xh5022.pdf
H = 8e2; % conversion constant for electrons -  signalout = ideal * exp(-log(2)*dose*A/H);

loss = zeros(size(vol)); noise = loss;
if opt.byslice
    for i=1:size(vol,3)
        sc = (i+1)/size(vol,3); % 1 placeholder for pre-exposures?
        r = rads*sc; 
        f = ceil(10/pix+r/10)*2+1;
        noise(:,:,i) = r*pix*5*randn(size(vol,[1,2])); % gaussian component
        decay = 1-(1-exp(-dose*sc*pix/H))/2;
        loss(:,:,i) = imgaussfilt3(vol(:,:,i),r/10,'FilterSize',f)*decay;
    end
else
    f = ceil(10/pix+rads/10)*2+1;
    noise = rads*pix*5*randn(size(vol)); % gaussian component
    decay = 1-(1-exp(-dose*rads/10*pix/H))/2;
    loss = imgaussfilt3(vol,rads/10,'FilterSize',f)*decay;
end

out = loss+noise; % loss = smoothing component, noise = gaussian component

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