function [labelmask,fields] = helper_watershed(vol,mask)

bin = imbinarize(rescale(vol));
cl = imclose(bin,strel('sphere',1)); %close small holes to avoid chopping most particles up
%cl = bin;
%cl = abs((bwdist(cl)>1)-1);
%cl = imopen(cl,strel('sphere',2));
%cl = imclose(cl,strel('sphere',2));
%is there a better way to safely dilate, by smaller than whole pixels?

%can i use an additional input as the mask/fieldshed?
%watershedding the inverse conv did a pretty good job, imhmax helped merge some choppy fields

d = -bwdist(~cl); %calculate distances for watershed
if nargin<2
mask = imextendedmin(d,2); %generate extended minima mask from distance map (not sure why 2)
end
%mask = bwareaopen(imerode(cl,strel('sphere',2)),4); %might work better?
d2 = imimposemin(d,mask); %merge local minima to make segmentations less chopped
w = watershed(d2,26); fields = single(w);
w(~bin) = 0; %watershed the local minima and mask out non-particles

labelmask = single(w);
end