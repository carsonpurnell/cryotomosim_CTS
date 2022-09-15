function mem = helper_membranegen(ts)
vol = ts.vol;
pix = ts.pix;

%generate a membrane field
q = rand(size(vol)); %random field noise
%do something with an anisotropic strel? make vertical membrane structures?
sm = imgaussfilt3(q,(6+80/pix(1))); %smooth noise to simplify membrane
bin = imbinarize(rescale(sm/1)); %binarize noise into random curvature
%unpredictable white/black space, so behavior when subtracting grid it inconsistent
m = zeros(size(bin)); m(200:300,200:300,1:end) = 1;
bin = single(bin).*m;

%subtract existing (dilated) structures from binary
dil = imdilate(imbinarize(vol),strel('sphere',round(1+(50+pix(1))/pix(1))));
b = imbinarize(bin-dil);
b = imgaussfilt3(single(b),6);
b = imbinarize(b);
%need membranes to be smoother, maybe smooth and binarize again to curve things out and clean up

%maybe do the trim/pad thing, not sure if necessary. generates more membrane planes if nothing else

%edge detect and expand to correct membrane thickness
g = imgradient3(b); %
memmask = imdilate(g,strel('sphere',1));
variance = 200;
baseline = 1000;
mem = memmask.* (randn(size(vol))*variance+baseline);
mem = rescale(mem,0,1400);
%ts.vol = ts.vol+mem;

%write to struct
%ts.misc.membrane = mem; ts.vol = vol
end

% 
% %membrane shenanigans
% pix = [13.6 16];
% 
% vol = zeros(400,400,50);
% q = rand(size(vol));
% %q(10:40,10:40,:) = q(10:40,10:40,:)/2;
% %q(380:390,10:390,20:30) = q(380:390,10:390,20:30)/2;
% 
% %do something with an anisotropic strel? make vertical membrane structures?
% sm = imgaussfilt3(q,(6+100/pix(1))); %smooth noise to simplify membrane
% bin = imbinarize(rescale(sm/1)); %binarize noise into random curvature
% %unpredictable white/black space, so behavior when subtracting grid it inconsistent
% m = zeros(size(bin)); m(100:300,100:300,1:end) = 1;
% bin = single(bin).*m;
% 
% %subtract existing (dilated) structures from binary
% dil = imdilate(imbinarize(vol),strel('sphere',round(4+(50+pix(1))/pix(1))));
% b = imbinarize(bin-dil);
% b = imgaussfilt3(single(b),6);
% b = imbinarize(b);
% %need membranes to be smoother, maybe smooth and binarize again to curve things out and clean up
% 
% %maybe do the trim/pad thing, not sure if necessary. generates more membrane planes if nothing else
% 
% %edge detect and expand to correct membrane thickness
% %g = rescale(imgradient3(b),0,1400); %
% e = rescale(edge3(b,'approxcanny',0.6),0,1400);
% 
% g = imgradient3(b); %
% memmask = imdilate(g,strel('sphere',1));
% variance = 200;
% baseline = 1000;
% mem = memmask.* (randn(size(vol))*variance+baseline);
% mem = rescale(mem,0,1400);
% vol = vol+mem;

%mem = imdilate(e,strel('sphere',1));

% sm = imgaussfilt3(q,10);
% bin = imbinarize(sm);
% %bin(365:390,10:390,10:40) = 0;
% trim = bin(2:end-1,2:end-1,2:end-1);
% pad = padarray(trim,[1 1 1]);
% close = imclose(pad,strel('sphere',6));
% g = rescale(imgradient3(close),0,1300);
% b = edge3(close,'approxcanny',0.6);
% mem = imdilate(g,strel('sphere',1));

%dilate/close?
%generate random field by exclusion of existing stuff (grid etc)?
%generate initial, then binarize and close to make mostly border with a few vesicles/folds inside?