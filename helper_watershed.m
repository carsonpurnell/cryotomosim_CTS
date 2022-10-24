function labelmask = helper_watershed(vol)
%
%
%vol = riboseg;
cl = imclose(vol,strel('sphere',2)); %close small holes to avoid chopping most particles up
cl = imerode(cl,strel('sphere',1));
%cl = imclose(cl,strel('sphere',2));
%is there a better way to safely dilate, by smaller than whole pixels?

d = -bwdist(~cl); %calculate distances for watershed
mask = imextendedmin(d,1); %generate extended minima mask from distance map (not sure why 2)
d2 = imimposemin(d,mask); %merge local minima to make segmentations less chopped
w = watershed(d2); %w2 = w;
w(~cl) = 0; %watershed the local minima and mask out non-particles

labelmask = w;

end