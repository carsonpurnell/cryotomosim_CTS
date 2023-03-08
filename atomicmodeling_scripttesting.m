% script to test if atomistic model-building is viable
pix = 10; clear particles;
tic
%particles(1) = helper_pdb2dat('Canhydrase_4xix_dimer.cif',pix,2,0,0);
particles(1) = helper_pdb2dat('tric__tric__6nra-open_7lum-closed.group.pdb',pix,2,0,0);
particles(2) = helper_pdb2dat('ribo__ribo__4ug0_4v6x.group.pdb',pix,2,0,0);
particles(3) = helper_pdb2dat('actin__6t1y_13x2.pdb',pix,2,0,0); %duplicate points warning
for i=1:numel(particles)
    for j=1:numel(particles(i).atomid)
        %particles(i).atomint{j} = atomdict(particles(i).atomid{j},'sc');
        com = mean(particles(i).atomcoords{j},1); %need radius from geometric, not mass center
        particles(i).atomcoords{j} = particles(i).atomcoords{j}-com;
        %size(particles(i).atomcoords{j})
        %size([0,0,0])
        %particles(i).radius{j} = max(pdist2(particles(i).atomcoords{j},single([0,0,0])));
        alphat = alphaShape(double(particles(i).atomcoords{j}),12); %surprisingly slow
        [~,p] = boundaryFacets(alphat);
        n = size(particles(i).atomcoords{j},1);
        ix = randperm(n); ix = ix(1:round(n/400));
        pi = particles(i).atomcoords{j}(ix,:);
        p = single([p;pi]); %need to add back 1-3% or so of points to prevent inside placements
        particles(i).perim{j} = unique(p,'rows');
        %clear com alphat %p pi
    end
end
toc

%{
%% load some input data - ribos
%{
[vol,sumvol,names,data] = helper_pdb2vol('ribo__ribo__4ug0_4v6x.group.mat',pix,2,1,0);
%might need to make atom id vector a row vector for memory write reasons
pts = single(data{1,2})'; atomid = (data{1,1});
pts = pts-mean(pts,1); %get and center points
%need atomdict function that accepts vector of atomic symbols and returns vector of Z/e- values
r = max(pdist2(pts,single([0,0,0]))); %find longest distance from the centroid among all points
atomint = atomdict(atomid);
alphat = alphaShape(double(pts),pix*1.2); [bf,p] = boundaryFacets(alphat);
particles(1).name = 'ribo'; particles(1).rad = r;
particles(1).pts = pts; particles(1).id = atomint; particles(1).perim = p;
%}

%% load some input data - tric
%{
[~,~,~,data] = helper_pdb2vol('tric__tric__6nra-open_7lum-closed.group.mat',pix,2,1,0);
%might need to make atom id vector a row vector for memory write reasons
pts = single(data{1,2})'; pts = pts-mean(pts,1); %get and center points
%atomid = (data{1,1}); %atomid = atomdict(atomid);
%need atomdict function that accepts vector of atomic symbols and returns vector of Z/e- values
%r = max(pdist2(pts,single([0,0,0]))); %find longest distance from the centroid among all points
alphat = alphaShape(double(pts),pix*1.2); [~,p] = boundaryFacets(alphat);
particles(2).name = 'tric'; 
particles(2).rad = max(pdist2(pts,single([0,0,0])));
particles(2).pts = pts; 
particles(2).id = atomdict(data{1,1}); 
particles(2).perim = p;
%}
%{
%% get boundary points to use for searching instead of the whole set
%so far does offer massive speed increase for the test particle - what about the rest of the model?
alphat = alphaShape(double(pts'),pix*1.2); %shape requires double for some reason
[bf,p] = boundaryFacets(alphat);
%plot(alphat)
%plot3(p(:,1),p(:,2),p(:,3),'.'); axis equal
%}
%}
%% constraint border atoms on planes
%first try just flat top/bottom z planes
%might try wavy version, definitely figure out x/y implementation as well
%obviates need for starting points in the model too
boxsize = pix*[200,300,50];
%edgedims = 3;

%% randomly add to the points and concatenate them into a list
boxsize = pix*[200,300,50];
%modelpoints = pts+boxsize/2; modelid = atomid;
%modelpoints =  modelid = 0; modelid2=single(modelid); 
dynpts = single([-100 -100 -100]); %dynid = single(0);
%modeltree = KDTreeSearcher(modelpoints');
rng(5);
tol = 2; %tolerance for overlap testing
count.s = 0; count.f = 0;
n = 1000; %ixcat = 2; %iscat 1 to erase the initial point
ixincat = 1; %index 1 to overwrite the initial preallocation point, 2 preserves it
split = cell(numel(particles)+1,2); split{1,1} = single([0,0,0]); split{1,2} = 0;
%zz = 0;
tl = tic; %ac = [];
for i=1:n
    if rem(i,n/20)==0; fprintf('%i,',i); end
    
    which=randi(numel(particles));
    sel = particles(which);
    sub = randi(numel(sel.atomcoords));
    %{
    pts = particles(which).pts;
    p = particles(which).perim;
    r = particles(which).rad;
    atomid = particles(which).id;
    name = particles(which).name;
    %}
    
    loc = rand(1,3).*boxsize;
    tform = randomAffine3d('rotation',[0 360]); 
    %tpts = tformfwd(tform,pts);
    %tpts = transformPointsForward(tform,pts')'+loc; %might need lower-level for speed, but seems fast
    %dummy = alphat; dummy.Points = ovcheck; %alpha shape modification is super slow
    %check if the models overlap, tolerance of 2A
    %precalc the maximum radius of the particle and do the model search only on points inside the radius?
    
    %modeltree = KDTreeSearcher(modelpoints'); 
    %make the model outside the loop at first, and only update after adding new points?
    %modeltree = ExhaustiveSearcher(modelpoints'); %exhaustive appears DRAMATICALLY slower
    
    %prune the large model points to only include within N radius (mode 2/sphere is ~4x faster)
    %[ix] = rangesearchpts(loc,r,modelpoints'); %preprune model points by radius to speed proximity testing
    %still a bottleneck - is send the points the slow part?
    %nested function probably speeds this up when modelpoints becomes large
    %4v7 when both used, still 4 when only using funct itself
    %1K iters: 127 v 199, funct still faster! (when not using nested nonpassing)
    % rangesearch inline implementation supposed to be faster, but is slower!
    %{
    s=0;
    for dd=1:numel(loc)
        s=s+(dynpts(:,dd)-loc(dd)).^2; %slower than rangesearch subfunct
    end
    fidx=s<sel.rad^2;
    ix5=find(fidx);
    %}
    %potential better vector solution - need to skip tons of 0,0,0 pts too from dyncat operation
    %spts = dynpts(sum(dynpts,2)~=0,:); %drop padded 0,0,0 points?
    %s = sum((dynpts(1:ixincat,:)-loc).^2,2);
    %ix5=find(s<(sel.radius{sub}+2)^2);% & s~=sum(loc.^2));
    %ixlin = s<sel.radius{sub}^2; %linear index version, much faster compute (way larger array though)
    %ix = rangesearchnest(loc,sel.radius{sub},dynpts(1:ixincat,:),tol); %how is it faster! how!!!
    
    ovcheck = transformPointsForward(tform,sel.perim{sub})+loc; %transform test points
    err = proxtest(dynpts(1:ixincat-1,:),ovcheck,tol); %prune and test atom collision
    %need to replace either with mutable quadtree or short-circuit kdtree
    %inlined: 33.5
    %{
    %idxnew = proxfilt(dynpts(1:ixincat,:),ovcheck,tol); %unfortunately slower than the radial search, but only 3x
    ix = proxfilt2(dynpts(1:ixincat-1,:),ovcheck,tol); %slower than radial, but closer search so overall faster
    %radial time: 119s 206s
    %box filt time: 99s 149s
    %vector time: 186s
    %{
    %ix = rangesearchvect(loc,sel.radius{sub},dynpts(1:ixincat,:));
    %qq = dynpts(ixlin); %linear prune speed test
    %qw = dynpts(ix5); %linear index is faster than logical overall
    %unique(round()) looks like it takes longer per run than radial search+KDT search
    
    %randomize the point order to try an exhaustive short-circuit search
    %shufpts = dynpts(randperm(size(dynpts,1)),:); %slower than searching, need to use random loop order
    
    %is there a better search version that just catches points inside of a convex hull?
    %can inshape accept a large array of points? might be problematic, points might not be INSIDE
    %easy way of making alphashape hulls just a tiny bit bigger? 
    %premake larger alphashape shell by multiplying perimeter points by ~1.1ish?
    %can alphashapes be transformed the same as points? yes, but super slow to read/write the memory
    %}
    err=0; %with n=100 exhaustive is only slightly slower than kdtree search, but progressive slowdown
    if ~isempty(ix) %this thing is taking SO VERY LONG, need more pre-optimization
        %{
        for j=randperm(size(dynpts(ix,:),1)) %shotgun rangesearch - slow
            scs = rangesearchnest(dynpts(j,:),2,ovcheck);
            if ~isempty(scs), break; end
        end
        for j=randperm(size(ovcheck,1)) %still crazy slow
            [scs,scsd] =  rangesearch(modeltree,ovcheck(j,:),5,'SortIndices',0);
            if any(scsd{:}<2), break; end
        end
        for j=randperm(size(dynpts(ix,:),1)) %shotgun pdist attempt - very slow
            scs = pdist2(ovcheck,dynpts(j,:),'euclidean','Smallest',1);
            if scs<2, break; end
        end
        %if any(kdist<2), pderr=1; else pderr=0; end %no errors but still super slow
        %}
        buck = round( size(dynpts,1)/450 );
        modeltree = KDTreeSearcher(dynpts(ix,:),'Bucketsize',buck); %67 with 1K %32 with 10K, 18 100K
        [~,d] = rangesearch(modeltree,ovcheck,tol,'SortIndices',0); %?? 1K,11.4 10K, 85 100K
        %the range search is now the slow part, ~40% of runetime for 2K iters
        d = [d{:}]; if any(d<tol), err=1; end %test if any points closer than 2A
        %if err~=pderr, zz=zz+1; end
    end
    %}
    
    %{
    %if numel(ix)>5, err=1; end
    %[ix,d] = knnsearch(modeltree,ovcheck,'K',1,'SortIndices',0); %does the sort make a difference?
    %d = [d{:}]; 
    %if any(d<2), err=1; end
    %if d(1)<10, err=1; end
    %[k,d] = dsearchn(modelpoints',tpts'); %INSANELY slow for some reason
    %[d,iix] = pdist2(modelpoints',tpts','euclidean','Smallest',5); %very very slow
    %if d(1)<10, err=1; end
    %}
    %{
    for j=1:numel(tpts) %also appears very slow, the loop is not at all optimal
        d = sum(abs(modelpoints-tpts(:,1)),1);
        if any(d<10), err=1; break; end
    end
    %}
    
    if err==0
        tpts = transformPointsForward(tform,sel.atomcoords{sub})+loc;
        %modelid = vertcat(modelid,atomid); %68.6s, slow overhead vertcat appears slower than horzcat
        %[modelid2,ixcat] = dyncathorz(ixcat,modelid2,sel.id); %slightly slower than hard cat in normal case
        %test out inline version to check if non-pass version is faster equiv to nested function
        
        %inlined dyncat - ~10x faster than cat, 15x faster than dyncat call (w 1000 iters)
        %make a bit more flexible by doubling the array size each hit? fewer expand steps, less overhead?
        %l = size(sel.atomint{sub},2); 
        
        % % inlined dyncat code % %
        l = size(ovcheck,1); e = ixincat+l-1;
        if e>size(dynpts,1)
            %dynid(:,ixincat:ixincat-1+l*10) = 0; %not used, dyn is temporary only
            %dynpts(ixincat:ixincat-1+l*10,:) = 0; %47,159
            %dynpts(ixincat:size(dynpts,1)*2,:) = 0; %44,176
            dynpts(ixincat:(size(dynpts,1)+l)*3,:) = 0; %43,153
        end
        %dynid(:,ixincat:e) = sel.atomint{sub};
        dynpts(ixincat:e,:) = ovcheck; %MUCH faster placing only perimeter points - need to prevent holing
        ixincat = ixincat+l;
        % % inlined dyncat code % %
        
        %janky offset stuff, make 'background' a particle class? or prepend the 4d vol with a zero vol?
        %ice will be the index 1 class!
        split{which+1,1} = [split{which+1,1};tpts]; %splitvol add
        split{which+1,2} = [split{which+1,2},sel.atomint{sub}];
        
        count.s=count.s+1;
    else
        count.f=count.f+1;
    end
end
toc(tl)
%dynid(ixincat:end) = []; %clear unused space from dynamic alloc vector
%dynpts(ixincat:end,:) = [];
disp(count)

%{
%% test make alphashape and boundary from the final model
%alp = alphaShape(modelpoints',pix*2);
%crazy, ludicrously slow. completely unusable with so many points and iterations.
%}
%{
%%  round pts to vol indices and add density to vol
volpts = round(modelpoints/pix+0.5); volid = modelid; %shift so things round well
emsz = floor(boxsize/pix); em = zeros(emsz);
for i=1:3
    ix = find(volpts(:,i) > emsz(i) | volpts(:,i) < 1);
    volpts(ix,:) = [];
    volid(:,ix) = [];
end
%atomint = atomdict(volid);
atomint = volid;
for i=1:numel(volid)
    %x=volpts(1,i); y=volpts(2,i); z=volpts(3,i); %fetch individual coordinates
    x=volpts(i,1); y=volpts(i,2); z=volpts(i,3);
    %c = volpts(:,i); x=c(1); y=c(2); z=c(3); %simultaneous pull is slower for some reason
    try
    em(x,y,z) = em(x,y,z)+atomint(i);
    catch
        [x,y,z], atomint(i) %#ok<NOPTS>
    end
end
%sliceViewer(em);
%}

%% functionalized volume projection
tic
%em = fnpt2vol(12,dynpts,dynid,boxsize);
%need another option for changing the location of the box from starting at 0,0,0
clear splitvol
for i=1:size(split,1)
    splitvol(:,:,:,i) = fnpt2vol(pix,split{i,1},split{i,2},boxsize,[0,0,0]);
end
[~,atlas] = max(splitvol,[],4); 
em = sum(splitvol,4); %small difference with em in a few scattered points
%atlas2 = min(atlas,proj); %zero out empty voxels to separate 0 from the first split
toc
sliceViewer(em);


%% solvation testing
rng(7)
tic
allatoms = vertcat(split{2:end,1});
distfrac = 0.5;
[solv] = gen_solvate(allatoms,boxsize,distfrac,tol);
toc
%sliceViewer(em+solv);
%{
%~33A^3 per H20 particle, rounding to 35
tic
distfrac = 1; %would need to scale by pixel size, .5 appears fine for 10a
waters = prod(boxsize/2)/35*distfrac; %number of predicted waters to fill the box
%make a very slightly larger box for solvation placement to avoid edge effects - 5-10A
solv = rand(round(waters),3,'single').*boxsize/2; %~1.2-1.4G per 100 million coords
%something to corral the ice between two slightly bending functions for an uneven surface?
allatoms = vertcat(split{2:end,1}); %concatenate all atoms into the same list
%break solv into multiple bins to search individually - generate several layers of coords/sort over loop?
%what's the optimal bucket size?
buck = round( size(allatoms,1)/450 ); %very rough early optimized value
modeltree = KDTreeSearcher(allatoms,'Bucketsize',buck); %make model for all the existing points
[idx,d] = rangesearch(modeltree,solv,tol/2,'SortIndices',0); %find waters near extant atoms - super slow
%iteratively remove waters one split at a time?

%idx is also in terms of model. need to build the model from waters, and remove nonempty IDX
%very very slow. increase bucket size to speed up? or sort and split the W into 2^n bins?
%use ~cellfun(@isempty,idx) to find proximity waters when waters are query points - but is that slow?
%almost certainly slower, but could run in serial or parallel with W bins against the same model
%alternatively vectorize all the cells from water query atoms, run unique, and remove those waters?
ol = ~cellfun(@isempty,idx); %waters that had a proximity fail
solv = solv(~ol,:); %remaining waters that are not overlapping particles
toc
%}

%% 
watervol = fnpt2vol(pix,solv,ones(1,size(solv,1))*2/distfrac,boxsize);
%sliceViewer(watervol)
sliceViewer(em+watervol);
%}

%% internal functions
function [tmp] = gen_solvate(modpts,sz,distfrac,tol)
h20vol = 35;
atomfrac = 1;
loops = 1; %appears to make no difference, may at full solvation size
buck = round( size(modpts,1)/450 ); %very rough early optimized value
modeltree = KDTreeSearcher(modpts,'Bucketsize',buck); %make model for all the existing points
%distfrac = 0.5; %would need to scale by pixel size, .5 appears fine for 10a
%so very very slow. gets too rough with denser pseudoatoms though
%shift density towards flat value, mask via atomic solv volgen?
waters = prod(sz/4)/h20vol*distfrac*atomfrac/loops; %number of predicted waters to fill the box
%pix = 10; flatd = (pix^3)/h20vol;
%make a very slightly larger box for solvation placement to avoid edge effects - 5-10A
solv = zeros(0,3);
for i=1:loops
tmp = rand(round(waters),3,'single').*sz/4; %~1.2-1.4G per 100 million coords
%pad = 20; spacing = 3/distfrac;
%[x,y,z] = meshgrid(1-pad:spacing:sz(1)/4+pad,1-pad:spacing:sz(2)/4+pad,1-pad:spacing:sz(3)/4+pad);
%solv = [x(:),y(:),z(:)]; %weird grids due to interpolation/binning shenanigans
[idx] = rangesearch(modeltree,tmp,tol/2,'SortIndices',0); %find waters near extant atoms - super slow
ol = ~cellfun(@isempty,idx); %waters that had a proximity fail
tmp = tmp(~ol,:); %keep only nontouching waters
solv = [solv;tmp];
end

end

function [v,ix] = dyncat(ix,v,b)
    l = size(b,1); e = ix+l-1;
    if e>size(v,1)
        v(ix:ix-1+l*10,:) = 0; %
    end
    v(ix:e,:) = b;
    ix = ix+l;
end

function [v,ix] = dyncathorz(ix,v,b)
    l = size(b,2); e = ix+l-1;
    if e>size(v,2)
        %fprintf('dynpad')
        v(:,ix:ix-1+l*100) = (0); %0 is by default double for some reason
    end
    v(:,ix:e) = b;
    ix = ix+l;
end

function vol = fnpt2vol(pix,pts,atomint,sz,offset)
if nargin<5, offset=[0,0,0]; end
if nargin<4, sz = max(pts,[],1)+pix; end
volpts = round((pts-offset)/pix+0.5); %shift so things round well
emsz = floor(sz/pix); vol = zeros(emsz);
for i=1:3
    ix = find(volpts(:,i) > emsz(i) | volpts(:,i) < 1);
    volpts(ix,:) = [];
    atomint(:,ix) = [];
end
for i=1:numel(atomint)
    %x=volpts(1,i); y=volpts(2,i); z=volpts(3,i); %fetch individual coordinates
    x=volpts(i,1); y=volpts(i,2); z=volpts(i,3);
    %d = volpts(i,:); %also slower
    %c = volpts(:,i); x=c(1); y=c(2); z=c(3); %simultaneous pull is slower for some reason
    vol(x,y,z) = vol(x,y,z)+atomint(i);
end
end

function idx = rangesearchvect(c,r,pts) %slower, likely memory alloc problems vectorized
s = sum((pts-c).^2,2);
idx=find(s<(r+2)^2);
end
function idx = rangesearchnest(c,r,pts,tol)
s=0;
for d=1:numel(c)
    s=s+(pts(:,d)-c(d)).^2;
end
%fidx=s<r^2; %7.63 split into two lines
%idx=find(fidx); %faster than low-level linear index
idx=find(s<(r+tol)^2); %7.723
%idx = 1:numel(fidx); %idx = idx(fidx);
%test if it is net faster to do indexing of modelpoints rather than creating idx - just return valid
%modelpoints?
end

%slower, also seems to not work properly. bad indices.
function ix = proxfilt(c,pts,tol)
%ix = ones(size(c,1),1);
idx = 1:size(c,1); %the vector to filter at the end and through loops
ix = idx; %the initial search vector
l = min(pts,[],1)-tol; %mins in each dim
h = max(pts,[],1)+tol; %maxes in each dim
for d=1:numel(h)
    ix = idx( (c(ix,d)>l(d)) & (c(ix,d)<h(d)) );
    ix = idx(ix);
end
%numel(ix)
%numel(ix>0)
end
function ix = proxfilt2(c,pts,tol)
l = min(pts,[],1)-tol; %mins in each dim
h = max(pts,[],1)+tol; %maxes in each dim
ix = c>l & c<h; %a = prod(a,2);
ix = find(sum(ix,2)>2);
%a = sum(a,2); ix = find(a);
%ix2 = find(a>0);
%{
ix = 1:size(c,1);
ix = ix.*a';
ix = ix(ix>0);
%}
end
function err = proxtest(c,pts,tol)
l = min(pts,[],1)-tol; %mins in each dim
h = max(pts,[],1)+tol; %maxes in each dim
ix = c>l & c<h; %a = prod(a,2);
ix = find(sum(ix,2)>2);
err=0; %with n=100 exhaustive is only slightly slower than kdtree search, but progressive slowdown
if ~isempty(ix) %this thing is taking SO VERY LONG, need more pre-optimization
    buck = round( size(c,1)/450 );
    modeltree = KDTreeSearcher(c(ix,:),'Bucketsize',buck); %67 with 1K %32 with 10K, 18 100K
    [~,d] = rangesearch(modeltree,pts,tol,'SortIndices',0); %?? 1K,11.4 10K, 85 100K
    %the range search is now the slow part, ~40% of runetime for 2K iters
    d = [d{:}]; if any(d<tol), err=1; end %test if any points closer than 2A
end
end

function [idx,dist]=rangesearchpts(c,r,X,mode)
% RANGESEARCH Range Search to find all points within a range.
% [index,distance]=rangesearch(c,r,X,mode) returns all points,x of X which
% are in the range of ||x-c||<r.
% Inputs:
%    c: 1 x d the querry vector
%    r: scalar defines the range
%    X: n x d array of all search points
% mode: range mode, 1 - box range, 2 (default) - radial range 
% Outputs:
%    index: indices of points in the range.
% distance: distances between the reference point and points in the range.
%
% See Also: pdist, kdtree, knnsearch
% Version 2.0 by Yi Cao at Cranfield University on 6th April 2008
%
%Examples
%Example 1: Radial range 
%{
X=rand(1000,2);
c=[0.5 0.5];
r=0.2;
idx=rangesearch(c,r,X);
t=0:0.02:2*pi;
x=c(1)+r*cos(t);
y=c(2)+r*sin(t);
subplot(211)
plot(X(:,1),X(:,2),'b.',c(1),c(2),'c+',X(idx,1),X(idx,2),'g*',x,y,'r-','linewidth',2)
title('Radial range search')
%}
%Example 2: Box range 
%{
X=rand(1000,2);
c=[0.5 0.5];
r=0.2;
idx=rangesearch(c,r,X,1);
t=0:0.02:2*pi;
x=c(1)+[-r r r -r -r];
y=c(2)+[r r -r -r r];
subplot(212)
plot(X(:,1),X(:,2),'b.',c(1),c(2),'c+',X(idx,1),X(idx,2),'g*',x,y,'r-','linewidth',2)
title('Box range search')
%}
% Example 3: Large data set search
%{
N=250000;
d=10;
X=randn(N,d);
c=randn(1,d);
r=1;
tic
idx=rangesearch(c,r,X);
toc
%}
% The time is similar to using kdtree (FEX ID 4586).
if nargin<4; mode=2; end
if mode==2
    [idx,dist]=range2(c,r,X);
else
    [idx,dist]=range1(c,r,X);
end
function [idx,dist]=range1(c,r,X)
[nPoints,nVariables]=size(X);
s=zeros(nPoints,1);
for d=1:nVariables
    x=abs(X(:,d)-c(d));
    id=x>s;
    s(id)=x(id);
end
fidx=s<r;
idx=find(fidx);
dist=s(fidx);
end
function [idx,dist]=range2(c,r,X)
nVariables=numel(c);
r2=r*r;
s=0;
for d=1:nVariables
    s=s+(X(:,d)-c(d)).^2;
end
fidx=s<r2;
idx=find(fidx);
dist=s(fidx);
end
end

function atomint = atomdict(atomid,mode)
if nargin<2, mode='z'; end
el = {'H','C','N','O','P','S','F','Na','MG','Cl','K','Ca','Mn','Fe'}; %element symbols to use for lookup
sc = [0.5288,2.5088,2.2135,1.9834,5.4876,5.1604,1.8012,4.7758,5.2078,4.8577,8.9834,9.9131,7.5062,7.1637,0];
z = [1,6,7,8,15,16,9,11,12,17,19,20,25,26,0]; %15==0 for bad entries
hp = [0,1.3,1.1,0.2,0,0.6,0,0,0,0,0,0,0,0,0]; %average hydrogens per atom
switch mode
    case 'z'
        atomint = z+z(1)*hp;
    case 'sc'
        atomint = sc+sc(1)*hp;
end
%scattering potentials computed as sum of first 5 parameters of atom form factor, holding s=0

[~,ix] = ismember(atomid,el);
%errs = find(ix<1); ix(errs) = 15;
ix(ix<1 | ix>15) = 15;
atomint = single(atomint(ix));
end