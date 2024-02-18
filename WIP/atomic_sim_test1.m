
% collect atomic model into single set of points
fn = fieldnames(split);
cen = boxsize/2;
atoms = zeros(0,4);
for i=1:numel(fn)
    atoms = [atoms;split.(fn{i})];
end

% rotate the model to an angle (eucentric adjustment? rotate about 0?)
angle = 25;
ax = [1,0,0];
atoms(:,1:3) = (atoms(:,1:3)-cen)*rotmat(ax,deg2rad(angle))+cen;

% project a set of slices - higher resolution in Z? start with isotropy
pix = 8;
slabthick = 4;
atoms(:,3) = (atoms(:,3)-min(atoms(:,3),[],'all'))/slabthick;
sz = boxsize; sz(3) = max(atoms(:,3),[],'all');
vol = helper_atoms2vol(pix,atoms,sz);

% get the transmission wave

% propogate transmission
param = param_simulate('tilt',zeros(1,size(vol,3)));
[convolved, ctf, param] = helper_ctf(vol,param);
% all NaN for some reason



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