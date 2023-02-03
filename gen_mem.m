function [blob,skel] = gen_mem(sz,thick,pix,beta)
sz = sz*pix; ptvol = zeros(sz,sz,sz); 
if numel(thick)==1, thick(end+1)=10; end
%thickness as nx2 vector of min,max membrane thickness?
thickness = thick(1)+rand*thick(2); r = (thickness/2)/pix; %actually computing the radius, not the diam


n = round(sz+sqrt(beta)); 
%might need to change to using points and generating an interpoland for better scaling
%keeping things on a fixed grid makes sizing difficult

%beta 3 is pretty good for quite round membrane bodies, 1 uniform, lower is U-shaped, high increases peak
%can decrease for wilder shapes, increase for nearer sphere
minl = sz*0.25; maxl = sz*0.75; %min and max of points for use as blob origins
pts = minl+(maxl-minl).*betarnd(beta,beta,n,3); %generate beta dist points to allow center/siding
pts = round(pts);
for i=1:size(pts,1)
    x = pts(i,1); y = pts(i,2); z = pts(i,3);
    ptvol(x,y,z) = 1;
end
blobvol = bwdist(ptvol)<min(sz)/r/2;
for i=1:3
    smoothvol = imgaussfilt3(single(blobvol),r*2-i*1); %smooth out towards a rounder overall shape
    blobvol = smoothvol>0.2;
end
memvol = ctsutil('trim',smoothvol>0.1); %trim vol to save space
p = round(r*2);
memvol = padarray(memvol,[p p p]); %pad to prevent anything from clipping into edges
%trim and pad the vol to center and make sure nothing is touching the edge

skel = bwperim(memvol);
distmap = bwdist(skel); 
mask = distmap<r;

smmask = imgaussfilt(single(mask),2);
sm2 = smmask.*mask;

memnoise = rand(size(skel))*0.3.*sm2; %sliceViewer(memnoise+sm2);
%need to compute density more precisely (variable component?) - also need control variables
dens = 0.35*pix^3;
blob = (memnoise+sm2)*dens;
end