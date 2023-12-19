function [field,layers] = helper_perlin(gridxy,pix,octaves,startoct)

arguments
    gridxy
    pix
    octaves = 10
    startoct = 8
end
%sz = size(grid);
%s = 0;%zeros(size(m));     % Prepare output image (size: m x m)
%octaves = 12;
%w = size(grid);
i = startoct-round(log2(pix));
e = i+octaves;
octaves = i:e;
%l = zeros(0);
%j = i:e;%,e,e]
pad = 10; sz = size(gridxy)+pad;
layers = zeros(sz(1),sz(2),numel(octaves));
%m = zeros(size(grid)+pad);
%change loop to a resolution-based metric - wavelength or frequency?
for i=1:numel(octaves)%while 2 > 1 && i<e
    %jj = j(i);
    d = randn(sz);
    for k=1:octaves(i) %do iterative refinements and shrinking rather than 2^k expansion in one step
        d = interp2(d, 1, 'spline');
        d = d(1:sz(1), 1:sz(2));
    end
    %d = rescale(d,-1,1);
    layers(:,:,i) = (1.32^i) *2* d(1:sz(1), 1:sz(2));
    %s = s + layers(:,:,i);
end
%s = s(pad+1:pad+size(grid,1),pad+1:pad+size(grid,2));
%l = l(:,:,2:end);
%s = (s - min(min(s(:,:)))) ./ (max(max(s(:,:))) - min(min(s(:,:))));
layers = layers(pad+1:pad+size(gridxy,1),pad+1:pad+size(gridxy,2),:);
field = sum(layers,3);
end