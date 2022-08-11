function [vec,cen,points] = shapeaxis(shape,threshold)
if nargin<2, threshold = 0; end

%strip out empty planes to speed things up
shape = shape(:,any(shape ~= 0,[1 3]),:); 
shape = shape(any(shape ~= 0,[2 3]),:,:); 
shape = shape(:,:,any(shape ~= 0,[1 2]));

%generate points from the volume
[X, Y, Z] = meshgrid(1:size(shape,2), 1:size(shape,1), 1:size(shape,3));
points = [X(:) Y(:) Z(:) shape(:)];
nonzero = points(:,4)>threshold;
points = points(nonzero,1:3);

%generate  vector from points
cen = mean(points,1);
[~,~,V] = svd((points-cen),'econ'); %singular value decomposition, econ for speed
vec = V(:,1)'; %the slope vector, ready for input to imrotate3

end