function norm4d = helper_volsurfnorm(skel,n)
%norm4d = helper_volsurfnorm(skel)
%creates a 4d array that contains the normal vectors across the 4th dimension for each point in skel, a 3d vol
%skel should be a 3d 'skeleton' - single-pixel thickness for best results. use bwperim to get a shape surface
%n is the number of points above/below the skeleton to use for estimating the vector. default 9 each
if nargin<2, n=9; end

perim = bwperim(bwdist(skel)<4); %dilate the skeleton
CC = bwconncomp(perim); %get the pixel arrays for each of the borders
numpixels = cellfun(@numel,CC.PixelIdxList); %count pixels in each component
[~,idx] = max(numpixels); %get the largest volume component from image
outer = perim*0;
outer(CC.PixelIdxList{idx}) = 1; %extract outer boundary
inner = perim-outer; %get inner boundary

% convert each volume to an array of points
[x,y,z] = ind2sub(size(skel),find(skel==1)); skelpts = [x,y,z];
[x,y,z] = ind2sub(size(skel),find(inner==1)); ptsin = [x,y,z];
[x,y,z] = ind2sub(size(skel),find(outer==1)); ptsout = [x,y,z];

%n = 9; %nearest n voxels on inner and outer surfaces to calculate vectors
%mskel = KDTreeSearcher(skelpts); %[ixself] = knnsearch(skelpts,skelpts,'K',5);
mdin = KDTreeSearcher(ptsin); [ixin] = knnsearch(mdin,skelpts,'K',n); 
mdout = KDTreeSearcher(ptsout); [ixout] = knnsearch(mdout,skelpts,'K',n);

norm4d = zeros(size(skel,1),size(skel,2),size(skel,3),3);
%normmat = zeros(size(idx,1),3); %normmat2=normmat1;

%vectorization to calculate vector targets for each point of skel
q = ptsin(ixin,:); q2 = reshape(q,[],3,n); win = sum(q2,3)/n;
q = ptsout(ixout,:); q2 = reshape(q,[],3,n); wout = sum(q2,3)/n;
%might need to change how points are initially generated to make it easier to reshape to something useful
ss = ixout(1:2,:)
skelpts(1,:)
size(q)
q(1:n,:)
dd = ptsout(ss,:)
%zz = reshape(dd,[],n,3)
%ff = permute(dd,[3,2,1])
%q2(1,:,:) %oriented wrong, sum across dims going between coordinate
wout(1,:)
%{
size(ptsout)
size(ixout)
size(q)
size(wout)
%}
%vectorize means by summing along a different dimension?
for i=1:size(ixin,1)
    skelcen = skelpts(i,:); incen = win(i,:); outcen = wout(i,:); 
    %outcen = mean(ptsout(ixout(i,:),:),1); %average of nearby outside points
    %incen = mean(ptsin(ixin(i,:),:),1); %average of nearby interior points - slow in loop
    long = outcen-incen; long = long/norm(long);
    under = skelcen-incen; under = under/norm(under);
    over = outcen-skelcen; over = over/norm(over); %average over and under, then average with long to refine?
    %refined = ((over+under)/2+long)/2; refined = refined/norm(refined);
    %refined = (over+under)/2; refined = refined/norm(refined);
    %refined = (refined+long)/2; refined = refined/norm(refined);
    %{
    if isnan(long)
        disp(long)
        disp(skelcen)
        disp(incen)
        disp(outcen)
    end
    %}
    
    vx = skelpts(i,1); vy = skelpts(i,2); vz = skelpts(i,3); %recover subscript data for the current point
    norm4d(vx,vy,vz,[1,2,3]) = long; %write normals to 4d storage array
end
nanchk = isnan(norm4d(:,:,:,1)); %find the NaNs
stack = nanchk+skel+outer+inner;
%sliceViewer(stack);


end