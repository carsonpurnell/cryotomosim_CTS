% matlab-only replacement for imod tiltproject

%% pick image
%[file, path] = uigetfile({'*.mrc'},'Select MRC');
fullpath = fullfile(path,file);
[img, head] = ReadMRC(fullpath);
inv = rescale(img*-1,00,255);

%{
%% imrotate3 - output size wierdness
theta = 15;
ax = [0,1,0];
rot = imrotate3(inv,theta,ax,'cubic','loose','FillValues',255);
proj = sum(rot,3);
%sliceViewer(rot);
%}

%% functionalized tilt projection
angles = -60:2:60;
ax = [1,0.0,0.0];
% how to generate spiral/circular processiong, or cumulative random walk near 0?

% what numerical scale for tilt images? not inverting to replicate current workflow
% current is inverted and scaled to original min/max before xyzproj, unscaling inversion
[tilts,rot,thick,vol] = tiltproj(inv,angles,ax);
% interpolation artifacts minor but visible - rotating atom set won't
%sliceViewer(rot)
sliceViewer(tilts);

%{
%% imwarp-based 3d rotation
theta = 15;
ax = [0.0,1,0.1];

%tformaff = randomAffine3d('rotation',[0 360]);
R = rotmataff(ax,deg2rad(theta));
%tform = affinetform3d(R)
tform = affine3d;
tform.T = R;

rout = affineOutputView(size(inv),tform,'BoundsStyle','followOutput');
%rout.XIntrinsicLimits = size(inv,1)+0.5;
%rout.YIntrinsicLimits = size(inc,2)+0.5;
outsz = size(inv); outsz(3)=rout.ImageSize(3);

%crop x/y back to original image space, only clip out in Z
xdiff = (rout.ImageSize(2)-size(inv,2))/2;
xworld = rout.XWorldLimits+[xdiff,-xdiff];
ydiff = (rout.ImageSize(1)-size(inv,1))/2;
yworld = rout.YWorldLimits+[ydiff,-ydiff];

imrefst = imref3d(outsz,xworld,yworld,rout.ZWorldLimits);

rot = imwarp(inv,tform,'linear','OutputView',imrefst,'FillValues',0);%mean(inv,'all'));
proj = sum(rot,3);
sliceViewer(rot);
%}

%% internal functions
function [tilts,rot,thick,vol] = tiltproj(vol,angles,ax)
tilts = zeros(size(vol,1),size(vol,2),numel(angles));
offeucentric = [0,0,20];
size(vol)
vol = padarray(vol,[0,0,100],'post');
size(vol)

for i=1:numel(angles)
    theta = deg2rad(angles(i));
    R = rotmataff(ax,theta);
    tform = affine3d; tform.T = R;
    
    rout = affineOutputView(size(vol),tform,'BoundsStyle','followOutput');
    %rout.XIntrinsicLimits = size(inv,1)+0.5;
    %rout.YIntrinsicLimits = size(inc,2)+0.5;
    outsz = size(vol); outsz(3)=rout.ImageSize(3);
    
    %crop x/y back to original image space, only clip out in Z
    xdiff = (rout.ImageSize(2)-size(vol,2))/2; xworld = rout.XWorldLimits+[xdiff,-xdiff];
    ydiff = (rout.ImageSize(1)-size(vol,1))/2; yworld = rout.YWorldLimits+[ydiff,-ydiff];
    
    imrefst = imref3d(outsz,xworld,yworld,rout.ZWorldLimits);
    
    % cubic and NN both seem dramatically slower
    rot = imwarp(vol,tform,'linear','OutputView',imrefst,'FillValues',0);%mean(vol,'all'));
    proj = sum(rot,3); %mean(proj,'all')
    thick = sum(rot>0,3); thick = max(thick,1);
    tilts(:,:,i) = proj;%./thick;
end
end

function t = rotmataff(ax,rad)
ax = ax/norm(ax);
x = ax(1); y = ax(2); z = ax(3);
c = cos(rad); s = sin(rad);

t1 = c + x^2*(1-c);
t2 = x*y*(1-c) - z*s;
t3 = x*z*(1-c) + y*s;
t4 = y*x*(1-c) + z*s;
t5 = c + y^2*(1-c);
t6 = y*z*(1-c)-x*s;
t7 = z*x*(1-c)-y*s;
t8 = z*y*(1-c)+x*s;
t9 = c+z^2*(1-c);

t = [t1 t2 t3 0
    t4 t5 t6 0
    t7 t8 t9 0
    0 0 0 1];
end