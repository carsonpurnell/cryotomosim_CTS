function surfaces = helper_surf(box,pix,angles,axis,sc)
%pix = 10;
%box = [400,300,40]*pix;
%angles = -60:5:60;
mang = max(abs(angles));
axspec = 1+rem(axis,2);
mtilt = tand(mang)*1.1; % angle for initial grid padding to cover tilt area - don't know why 0.6 is enough
% make surfaces - displaced in Z
% trying with identical surf pairs first
res = pix/2; %oversample to prevent holes in data
%need to stretch coverage a lot, fill values are static and not NN or average of near edge
pad = [pix*2,pix*2]; %axis = 1; % huge padding only against tilt axis, otherwise small
pad(axis) = round((box(axspec)+box(3))*mtilt)+pix*2; % huge padding to cover for tilting
[x,y] = meshgrid(pix/2-pad(1):res:box(1)+pad(1),pix/2-pad(2):res:box(2)+pad(2));

%sc = [2.5,1.2];
surfaces{1} = gensurf(x,y,box,pix,sc);
surfaces{2} = gensurf(x,y,box,pix,-sc);

end

function gridsurf = gensurf(x,y,box,pix,sc)
[field,layers] = helper_perlin(x,pix,sc(1),6,10);
mv = mean(field,'all')*0.5; 
if mv<0&&sc(2)>0; mv=mv*2-0*sqrt(abs(mv)); end
if mv>0&&sc(2)<0; mv=mv*2-0*sqrt(abs(mv)); end
field = field-mv;
gridsurf = [x(:),y(:),field(:)+box(3)*sc(2)/2];
end

function [field,layers] = helper_perlin(gridxy,pix,mag,octaves,startoct)
%[field,layers] = helper_perlin(gridxy,pix,octaves,startoct)
arguments
    gridxy
    pix
    mag = 2
    octaves = 10
    startoct = 8
end
%sz = size(grid); %w = size(grid);
adj = round(log2(pix)); %adjustment to keep noise octaves on the same magnitude across diff pixel sizes
i = startoct-adj; %adjust frequency based on pixel size
e = i+octaves;
octaves = i:1:i+octaves;
prep = i-1; %faster prelim interp up to before first octave
%l = zeros(0); %j = i:e;%,e,e]
pad = 0; sz = size(gridxy)+pad;
layers = zeros(sz(1),sz(2),numel(octaves));
%m = zeros(size(grid)+pad);
for i=1:numel(octaves)
    oc = octaves(i);
    % prep interpolation, start with a small grid and expand to size for speed?
    %prep = 64;
    d = randn(round(sz/(2^prep)+4));
    for k=1:prep
        d = interp2(d, 'spline');
    end
    d = d(1:sz(1), 1:sz(2));
    for k=prep+1:oc %do iterative refinements and shrinking rather than 2^k expansion in one step
        d = interp2(d, 'spline');
        d = d(1:sz(1), 1:sz(2));
    end
    layers(:,:,i) = (1.4^oc) *mag* d(1:sz(1), 1:sz(2));
    %s = s + layers(:,:,i);
end
%s = s(pad+1:pad+size(grid,1),pad+1:pad+size(grid,2));
%l = l(:,:,2:end);
%s = (s - min(min(s(:,:)))) ./ (max(max(s(:,:))) - min(min(s(:,:))));
layers = layers(pad+1:pad+size(gridxy,1),pad+1:pad+size(gridxy,2),:);
field = sum(layers,3);
end