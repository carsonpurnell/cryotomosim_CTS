function surfaces = helper_surfmill(box,pix,angles,axis,sc)
axspec = 1+rem(axis,2);
res = pix/2; %pad = pix*4;
mtilt = tand(max(abs(angles)))*1.5;
box(axis) = round((box(axis)+box(3))*mtilt); %box = box+4*pix;
pad = 4*pix;

%res = pix/2; % oversample resolution to prevent missing values in meshgrids
%padmult = 2; pval = pix*padmult; %pad extent, need to increase for non-ordinal tilt axis
%padding = [pval,pval]; %basic minor padding to cover edge rounding
%padding(axis) = round((box(axspec)+box(3))*mtilt)+pix*2; % huge padding to cover for tilting
%[x,y] = meshgrid(pix/padmult-padding(1):res:box(1)+padding(1),pix/padmult-padding(2):res:box(2)+padding(2));
[x,y] = meshgrid(1:res:box(1)+pad,1:res:box(2)+pad);

theta = rand*pi; %180* random angle, covers all meaningful space
dof = 0; if theta>pi/2, theta=theta-pi/2; y=flip(y,2); dof = 1; end %flip if encountering negative vals
mesh = round(sin(theta)*x+cos(theta)*y);
[a,b] = bounds(mesh,'all');

adj = [0,0,0]; adj(axis) = box(axis)/3*1;

[wave] = int_noise([a,b],sc,[9:13]);
z = wave(mesh); if dof==1; z=flip(z,2); end
surfaces{1} = [x(:),y(:),z(:)+box(3)/2]-adj-[pad,pad,0]/2;

[wave] = int_noise([a,b],sc,[9:13]);
z = wave(mesh); if dof==1; z=flip(z,2); end
surfaces{2} = [x(:),y(:),z(:)-box(3)/2]-adj-[pad,pad,0]/2;

end

function [wave,layers] = int_noise(sz,amp,oct)
wave = zeros(1,sz(2)); layers = cell(1,numel(oct));
for i=oct(1):oct(2)
    d = randn(size(wave));
    for k=1:i
        d = interp1(d,1:0.5:round(sz(2)/2+1),'spline');
        d = d(1:sz(2));
    end
    d=d*(i*1.5+1);
    layers{i+oct(1)+1} = d;
    wave = wave+d;
end
wave = wave*amp;
end