function [out,rad,rad2] = helper_radiation(vol,solv,rad)

% gaussian component
radg = (rand(size(vol))-rand(size(vol)))*rad;
% smoothing component - run on vol separately? no, would generate weirdness
rads = imgaussfilt3(vol,rad);
% separate filling/warping component? or suitably part of others?
% solv needed as a separate component at all?

rad = 0; rad2=0;
out = vol+solv;
end