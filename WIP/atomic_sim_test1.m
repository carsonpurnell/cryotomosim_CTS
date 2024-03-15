
% collect atomic model into single set of points
fn = fieldnames(split);
cen = boxsize/2;
atoms = zeros(0,4);
for i=1:numel(fn)
    atoms = [atoms;split.(fn{i})];
end

%%
angles = -60:10:60;
param = param_simulate('pix',pix,'tilt',angles);
[tilt,dtilt] = atomictiltproj(atoms,param,angles,boxsize,20);
sliceViewer(dtilt);

%% 
%{
% rotate the model to an angle (eucentric adjustment? rotate about 0?)
angle = 60;
ax = [1,0,0];
atoms(:,1:3) = (atoms(:,1:3)-cen)*rotmat(ax,deg2rad(angle))+cen;

% project a set of slices - higher resolution in Z? start with isotropy
pix = 8;
slabthick = 2;
atoms(:,3) = (atoms(:,3)-min(atoms(:,3),[],'all'))/slabthick;
sz = boxsize; sz(3) = max(atoms(:,3),[],'all');
vol = helper_atoms2vol(pix,atoms,sz);

%sim params
param = param_simulate('pix',pix,'tilt',zeros(1,size(vol,3)));
param.tilt = -40:numel(param.tilt)-41;

% get the transmission wave
d = param.dose*pix^2;
dvol = poissrnd(rescale(vol*-1)*d,size(vol));

% propogate transmission
mid = round(size(dvol,3)/2);
convolved = zeros(size(dvol));
for i=1:size(dvol,3)
    adj = (pix*slabthick*(i-mid))/1e4*3e1;
    param.defocus = -5+adj;
    param.tilt = 0;
    [convolved(:,:,i), ctf, param] = helper_ctf(dvol(:,:,i),param);
end
proj = rescale(sum(convolved,3));

imshow(proj);
%}
%% internal functs
function [tilt,dtilt] = atomictiltproj(atoms,param,angles,boxsize,slabthick)
ax = [1,0,0];
cen = boxsize/2;
%angles = param.tilt;
% get the transmission wave
d = param.dose/numel(param.tilt)*param.pix^2;

tilt = zeros(boxsize(1)/param.pix,boxsize(2)/param.pix,numel(param.tilt));
for t=1:numel(angles)
    angle = angles(t);
    atomtmp = atoms;
    atomtmp(:,1:3) = (atomtmp(:,1:3)-cen)*rotmat(ax,deg2rad(angle))+cen;
    
    % project a set of slices - higher resolution in Z? start with isotropy
    %pix = 8;
    %slabthick = 10;
    atomtmp(:,3) = (atomtmp(:,3)-min(atomtmp(:,3),[],'all'))/slabthick;
    sz = boxsize; sz(3) = max(atomtmp(:,3),[],'all');
    [vol,solv] = helper_atoms2vol(param.pix,atomtmp,sz);
    vol = vol+solv;
    %sim params
    %param = param_simulate('pix',param.pix,'tilt',zeros(1,size(vol,3)));
    %param.tilt = -40:numel(param.tilt)-41;
    
    % get the transmission wave
    %d = param.dose*param.pix^2;
    %dvol = poissrnd((vol*1)*d,size(vol)); %extremely slow with many sections - do at the end?
    
    % propogate transmission
    mid = round(size(vol,3)/2);
    convolved = zeros(size(vol));
    for i=1:size(vol,3)
        adj = (param.pix*slabthick*(i-mid))/1e4*3e1;
        param.defocus = -5+adj;
        param.tilt = 0;
        [convolved(:,:,i), ctf, param] = helper_ctf(vol(:,:,i),param);
    end
    
    tilt(:,:,t) = sum(convolved,3);
end
dtilt = poissrnd((d*rescale(tilt*1,0,1))*01,size(tilt));
end

function t = rotmat(ax,rad)
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

t = [t1 t2 t3
    t4 t5 t6
    t7 t8 t9];
end