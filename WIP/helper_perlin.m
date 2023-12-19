function [field,layers] = helper_perlin(grid,pix)
%sz = size(grid);
s = 0;%zeros(size(m));     % Prepare output image (size: m x m)
w = size(grid);
i = 8-round(log2(pix));
e = i+12;
l = zeros(0);
j = i:e;%,e,e]
layers = zeros(size(grid,1),size(grid,2),numel(j));
pad = 10;
m = zeros(size(m)+pad);
%i = 4;
%change loop to a resolution-based metric - wavelength or frequency?
for i=1:numel(j)%while 2 > 1 && i<e
    %i = i + 2;%1;
    jj = j(i);
    %i=j(jj);
    d = randn(size(m));
    for k=1:jj %do iterative refinements and shrinking rather than 2^k expansion in one step
        d = interp2(d, 1, 'spline');
        d = d(1:size(m,1), 1:size(m,2));
    end
    %d = rescale(d,-1,1);
    layers(:,:,i) = (1.35^i) *2* d(1:size(m,1), 1:size(m,2));
    s = s + layers(:,:,i);
    %w = w - ceil(w/2 - 1);
end
s = s(pad+1:pad+w(1),pad+1:pad+w(2));
l = l(:,:,2:end);
%s = (s - min(min(s(:,:)))) ./ (max(max(s(:,:))) - min(min(s(:,:))));
field = sum(layers,3);
end