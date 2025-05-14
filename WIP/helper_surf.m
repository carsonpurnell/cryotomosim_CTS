function surfaces = helper_surf(box,pix,angles,axis,sc)
%mang = max(abs(angles)); 
axspec = 1+rem(axis,2);
mtilt = tand(max(abs(angles)))*1.0; % compute grid padding to cover whole area of high tilt

res = pix/2; % oversample resolution to prevent missing values in meshgrids
padmult = 2; pval = pix*padmult; %pad extent, need to increase for non-ordinal tilt axis
padding = [pval,pval]; %basic minor padding to cover edge rounding
padding(axis) = round((box(axspec)+box(3))*mtilt)+pix*2; % huge padding to cover for tilting
[x,y] = meshgrid(pix/padmult-padding(1):res:box(1)+padding(1),pix/padmult-padding(2):res:box(2)+padding(2));

surfaces{1} = gensurf(x,y,box,pix,sc); %generate each field of 3d points from xy coord grid
surfaces{2} = gensurf(x,y,box,pix,-sc);
end

function gridsurf = gensurf(x,y,box,pix,sc)
[field,layers] = helper_perlin(x,pix,abs(sc(1)),6,10);
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
if mag == 0; field = zeros(size(gridxy)); layers = field; return; end
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
        %d = interp2(randn(ceil((n-1)/(2^(i-1))+1),ceil((m-1)/(2^(i-1))+1)), i-1, 'spline');
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