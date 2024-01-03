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
adj = round(log2(pix));
i = startoct-adj; %adjust frequency based on pixel size
e = i+octaves;
octaves = [4-adj,i:1:i+octaves];
%l = zeros(0); %j = i:e;%,e,e]
pad = 0; sz = size(gridxy)+pad;
layers = zeros(sz(1),sz(2),numel(octaves));
%m = zeros(size(grid)+pad);
for i=1:numel(octaves)
    oc = octaves(i);
    % prep interpolation, start with a small grid and expand to size for speed?
    prep = 16;
    d = randn(round(sz/prep+4));
    for k=1:log2(prep)
        d = interp2(d, 1, 'spline');
    end
    for k=log2(prep)+1:oc %do iterative refinements and shrinking rather than 2^k expansion in one step
        d = interp2(d, 1, 'spline');
        d = d(1:sz(1), 1:sz(2));
    end
    sz
    size(layers(:,:,i))
    size(d)
    layers(:,:,i) = (1.4^oc) *mag* d(1:sz(1), 1:sz(2));
    %s = s + layers(:,:,i);
end
%s = s(pad+1:pad+size(grid,1),pad+1:pad+size(grid,2));
%l = l(:,:,2:end);
%s = (s - min(min(s(:,:)))) ./ (max(max(s(:,:))) - min(min(s(:,:))));
layers = layers(pad+1:pad+size(gridxy,1),pad+1:pad+size(gridxy,2),:);
field = sum(layers,3);
end