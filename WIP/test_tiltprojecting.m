% matlab-only replacement for imod tiltproject

%% pick image
[file, path] = uigetfile({'*.mrc'},'Select MRC');
fullpath = fullfile(path,file);
[img, head] = ReadMRC(fullpath);
inv = rescale(img*-1,0,255);

%% imrotate3 - output size wierdness
theta = 15;
ax = [0,1,0];
rot = imrotate3(inv,theta,ax,'cubic','loose','FillValues',255);
proj = sum(rot,3);
%sliceViewer(rot);

%% imwarp
theta = 15;
ax = [0,1,0];


%tformaff = randomAffine3d('rotation',[0 360]);
R = rotmataff(ax,deg2rad(theta));
%tform = affinetform3d(R)
tform = affine3d;
tform.T = R;

rout = affineOutputView(size(inv),tform,'BoundsStyle','followOutput')
%rout.XIntrinsicLimits = size(inv,1)+0.5;
%rout.YIntrinsicLimits = size(inc,2)+0.5;
outsz = size(inv); outsz(3)=rout.ImageSize(3);
xdiff = rout.Imagesize(2)-size(inv,1)
ydiff
%outsz(3) = rout.ImageSize(3);
imrefst = imref3d(outsz,[0.5,size(inv,2)+0.5],rout.YWorldLimits,rout.ZWorldLimits)

rot = imwarp(inv,tform,'OutputView',imrefst,'FillValues',255);
%sliceViewer(rot);

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