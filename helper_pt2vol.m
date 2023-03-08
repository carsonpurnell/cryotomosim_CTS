function vol = helper_pt2vol(pix,pts,sz,offset)
%vol = helper_pt2vol(pix,pts,sz,offset)
%generates a volume from a given set of points depending on their density
%nx3 array uses magnitude==1 for each point, nx4 array uses the 4th column to set the magnitude of each point
%sz are the dimensions of the box to constrain the volume to. defaults to the one-way maximum from pts
%offset default to 0,0,0, and is the initial point from which the box originates. negative data points will be
%excluded (such as points centered around 0,0,0) unless offset is made sufficiently negative.
if nargin<4, offset=[0,0,0]; end
if nargin<3, sz = max(pts,[],1)+pix; end
if size(pts,2)<4, pts(:,end+1)=1; end %intensity==1 if not given by 4th column

pts(:,1:3) = round((pts(:,1:3)-offset)/pix+0.5);
emsz = floor(sz/pix); vol = zeros(emsz);
for i=1:3 %prune points that would be outside the box - looping over each dimension
    ix = pts(:,i) < emsz(i) & pts(:,i) > 1; %get points inside the box
    pts = pts(ix,:); %drop points outside the box
end
for i=1:size(pts,1) %accumulate the specified densitiy within each voxel
    x=pts(i,1); y=pts(i,2); z=pts(i,3); mag = pts(i,4); %fetch data per atom
    vol(x,y,z) = vol(x,y,z)+mag;
end

end