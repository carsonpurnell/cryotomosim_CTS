% matlab-only replacement for imod tiltproject

%% pick image
%[file, path] = uigetfile({'*.mrc'},'Select MRC');
fullpath = fullfile(path,file);
[img, head] = ReadMRC(fullpath);
inv = rescale(img*1,00,1);
param = param_simulate('pix',head.pixA); param.size = [head.mx,head.my,head.mz];

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
[tilts,rot,thick] = tiltproj(inv*1,angles,ax);
% interpolation artifacts minor but visible - rotating atom set won't
%sliceViewer(rot)
sliceViewer(tilts);

%%  dose/ctf
[detect,rad] = helper_electrondetect(tilts*1,param);
%sliceViewer(detect)
%%
%rs = rescale(tilts,min(img,[],'all'),max(img,[],'all'));
[convolved, ctf, param] = helper_ctf(detect,param);
sliceViewer(convolved*-1);

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
function [tilts,rot,thick] = tiltproj(vol,angles,ax)
tilts = zeros(size(vol,1),size(vol,2),numel(angles));
offeucentric = 50;%[0,0,20];
%%vol = padarray(vol,[0,0,50],'pre');

for i=1:numel(angles)
    volstep = padarray(vol,[0,0,offeucentric+randi(21)-11],'pre');
    theta = deg2rad(angles(i));
    R = rotmataff(ax,theta);
    tform = affine3d; 
    tform.T = R;
    
    rout = affineOutputView(size(volstep),tform,'BoundsStyle','followOutput');
    %rout.XIntrinsicLimits = size(inv,1)+0.5;
    %rout.YIntrinsicLimits = size(inc,2)+0.5;
    outsz = size(volstep); outsz(3)=rout.ImageSize(3);
    
    %crop x/y back to original image space, only clip out in Z
    xdiff = (rout.ImageSize(2)-size(volstep,2))/2; xworld = rout.XWorldLimits+[xdiff,-xdiff];
    ydiff = (rout.ImageSize(1)-size(volstep,1))/2; yworld = rout.YWorldLimits+[ydiff,-ydiff];
    
    imrefst = imref3d(outsz,xworld,yworld,rout.ZWorldLimits);
    
    % cubic and NN both seem dramatically slower
    rot = imwarp(volstep,tform,'linear','OutputView',imrefst,'FillValues',0);%mean(vol,'all'));
    %tr = [randi(21)-11,randi(21)-11,randi(21)-11]; %translation matrix
    %rot = imtranslate(rot,tr); % apply translation after rotation
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