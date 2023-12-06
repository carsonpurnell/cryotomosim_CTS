function mu = mu_build(points,init,mu,opt)
% mutable and isotopic balanced octree

arguments
    points
    init = []; % init bounds or ix vector?
    mu = [] %prior tree for expansion
    opt.leafmax = 1e3
    opt.leaflen = 25
    opt.loose = 0
    %opt.bounds = []
    opt.maxdepth = 10
    opt.pad = 0
end
points = single(points); init = single(init); %avoid pdist2 warnings and compress data slightly
%bm = zeros(8,3); for i=1:8; bm(i,:) = bitget(i,1:3); end %old bitmask, slower due to eventual shape needed
bm = zeros(1,3,8); for i=1:8; bm(1,:,i) = bitget(i,1:3); end %generate bitmask
if iscell(mu) %&& isfield(mu,'opt')
    %opt = mu.opt; %currently unused backstop to import options if existing tree is used
    rootcoords = vertcat(mu{1,1}{2,:});
    %cubel = mu{1,2};
    %mutate=1; 
    if ~isempty(points) && ~isempty(init)
    init=init';
    ix = unique(init,'rows');
    for i=1:size(ix,1)
        d=ix(i,1); bl=ix(i,2); %depth and layer of the relevant leaf
        ix2 = all(init==ix(i,:),2); %logical index for which pts to place in the leaf
        tmp = points(ix2,:); % collect points for the leaf
        mu{d,1}{1,bl} = [mu{d,1}{1,bl};tmp]; % add new pts to leaf list
        [l,h] = bounds(mu{d,1}{1,bl},1);
        mu{d,1}{3,bl} = [l;h]; %update leaf bounds
        mu = leafcheck(mu,opt,bm,d,bl);
    end
    end
else
    [rootcoords,cubel] = cubify(points,init);
    %mutate=0;
    [rootsort] = rootsplit(rootcoords,cubel,points,opt.pad);
    mu{1,2} = cubel;
    mu{1,1}{3,size(rootcoords,1)} = []; %initialize root level cells
    for i=1:size(rootcoords,1)
        %tps = points(rootix==i,:); %tps = rootsort{i}; %indexing is the slow part
        mu{1,1}{1,i} = rootsort{i}; %points
        mu{1,1}{2,i} = rootcoords(i,:); %center
        %mutree{1,1}{3,i} = []; %branches - now reusing corners for that
        [l,h] = bounds(rootsort{i},1); %breaks if only single point in the branch?
        mu{1,1}{3,i} = [l;h]; %corner coverage
        %leafcheck(1,i)
        mu = leafcheck(mu,opt,bm,1,i); %pass marginally faster so far
    end
end

% todo
% change top-level organization to struct. easier handling of opts between functions
% cube length only needs stored at root level, after that easily computed from depth
% mu.tree can be a nx1 cell array, shorter code to index into it
% can points cell be reused into branch list after splitting?

%{
%find bounds from points, and bounds if given
[l,h] = bounds([points;init]); l=round(l); h=round(h);
%cen = round((h+l)/2,-1); %centroid of volume, maybe useless right now
len = round(h-l,-1); %length of each axis

%cubifying does seem slightly faster overall, not sure how much or what context is most important
%need some metric to force breaking lowest axis length into smaller increments to ensure more roots
cubel = min(len); ch = cubel/2; %cube full and half length for mesh compaction
len = ceil(len/cubel)*cubel; %length of axes in cube segmentations

%grid not quite centered, starts at min and rounds past max values
xi = l(1)+ch:cubel:len(1)-ch; yi = l(2)+ch:cubel:len(2)-ch; zi = l(3)+ch:cubel:len(3)-ch;
[x,y,z] = meshgrid(xi,yi,zi); %grid of centers for root cubes

rootcoords = single(zeros(numel(x),3)); %reshape might be more elegant but hard to figure
for i=1:numel(x)
    rootcoords(i,:) = [x(i),y(i),z(i)]; %list of coordinate centers for root cubes
end
%}

%{
%sort points into closest root - no check for padding or outside range or multi inclusion yet
%[~,rootix] = pdist2(rootcoords,points,'squaredeuclidean','Smallest',1); %squared negligibly faster?
[d,rootix] = pdist2(rootcoords,points,'euclidean','Smallest',8); %get distances to root centers
c=d<=cubel/2+opt.pad; c(1,:)=1; %index distances in range, force nearest to always positive to ensure inclusion
rootix(~c)=0; %filter to only indices of valid roots
c = any(rootix,2); %vector for throwing out most of rootix which is unused (1, sometimes 2 rows)
rootix = rootix(c,:); %fastest reassign and trim?
%{
%b = rootix.*c;
%b = rootix; b(~c,:)=[];
%rootix(~c,:)=[];
%rootix=a;
%all(a==b,'all')
%rootix
%any(rootix==1,1)
%any(rootix==2,1)
%rootix([c==1])'
%rootix(c)
%rootix = rootix'; %make column for slightly faster indexing
%}

rootsort = cell(1,size(rootcoords,1));
for i=1:size(rootcoords,1)
    rootsort{i} = points(any(rootix==i,1),:); % indexing all the pts pretty slow
end
% % lump in with cubify code? or make separate handler for search to also use?
%}

%{
%nest-share, was inconsistently faster initially
    function leafcheck(level,layer)
        %level = depth , layer = which branch in the level
        splitcheck = size(mu{level,1}{1,layer},1)>opt.leafmax &&...
            mu{level,2}>(opt.leaflen*2) && level<opt.maxdepth;
        if splitcheck==1
            branch(level,layer)
        end
    end

    function branch(level,layer)
        if size(mu,1)<level+1 %initialize next layer if needed
            mu{level+1,2} = mu{level,2}/2;
        end
        tp = mu{level,1}{1,layer};%.points;
        tc = mu{level,1}{2,layer};%.center;
        tmplen = mu{level+1,2}/4;
        %cmp = (tp>tc);
        %[~,b] = ismember(cmp,bm,'rows'); %ismember ~3x slower than in-loop array math
        %[~,f,~] = intersect(cmp,bm,'stable','rows'); %intersect even slower than ismember
        [~,bin] = max(all(bsxfun(@eq,tp>tc,bm),2),[],3); % slightly faster than loop, still sub bottleneck
        bin=bin'; %row vector is very slightly faster
        prior = size(mu{level+1,1},2); %should be 0 for new layer
        bi = prior+1:prior+8; %branches in new level
        for i2=1:8
            il = prior+i2; %adjust layer index past existing entries
            bit = bitget(i2,1:3); %get bitmask comparator to filter where points are to a child branch
            %c2=~any(bsxfun(@minus,cmp,bit),2); %match is 000, convert is >10x faster than ismember
            %c=~any(cmp-bit,2); %implicit expansion faster than bsxfun for this
            %if ~all(c2==c), 'fail', end
            
            %rr = tp(bin==i2,:); %
            mu{level+1,1}{1,il} = tp(bin==i2,:); %points   %sub bottleneck, indexing slow
            %gg = (bit*2-1)*tmplen+tc;
            mu{level+1,1}{2,il} = (bit*2-1)*tmplen+tc; %center
            %mutree{level+1,1}{3,il} = []; %branches, now replaces corners after branch
            [l,h] = bounds(mu{level+1,1}{1,il});
            mu{level+1,1}{3,il} = [l;h]; %occupancy corners
            leafcheck(level+1,il)
        end
        mu{level,1}{1,layer} = []; %size(mutree(level).layer(layer).points,1);
        mu{level,1}{3,layer} = bi; %branch listing
        %mutree{level,1}{4,layer} = []; %clear parent corners - now using 3rd cell and replaced by branches
    end
%}
%{
    function populate(level,layer) %prefill new layers with initial values
        mutree{level,1}{1,layer}.points = [];
        mutree{level,1}{1,layer}.center = [];
        mutree{level,1}{1,layer}.branches = [];
    end
%}

end

%
%sub-pass, faster overall
function mu = leafcheck(mu,opt,bm,level,layer)
%level = depth , layer = which branch in the level
splitcheck = size(mu{level,1}{1,layer},1)>opt.leafmax &&...
    mu{level,2}>(opt.leaflen*2) && level<opt.maxdepth;
if splitcheck==1
    mu = branch(mu,opt,bm,level,layer);
end
end

function mu = branch(mu,opt,bm,level,layer)
if size(mu,1)<level+1 %initialize next layer if needed
    mu{level+1,2} = mu{level,2}/2;
end
tp = mu{level,1}{1,layer};%.points;
tc = mu{level,1}{2,layer};%.center;
tmplen = mu{level+1,2}/4;
%cmp = (tp>tc);
%[~,b] = ismember(cmp,bm,'rows'); %ismember ~3x slower than in-loop array math
%[~,f,~] = intersect(cmp,bm,'stable','rows'); %intersect even slower than ismember
[~,bin] = max(all(bsxfun(@eq,tp>tc,bm),2),[],3); % slightly faster than loop, still bottleneck
bin=bin'; %row vector is very slightly faster
prior = size(mu{level+1,1},2); %should be 0 for new layer
bi = single(prior+1:prior+8); %branches in new level
for i2=1:8
    il = prior+i2; %adjust layer index past existing entries
    bit = bitget(i2,1:3); %get bitmask comparator to filter where points are to a child branch
    %c2=~any(bsxfun(@minus,cmp,bit),2); %match is 000, convert is >10x faster than ismember
    %c=~any(cmp-bit,2); %implicit expansion faster than bsxfun for this
    %if ~all(c2==c), 'fail', end
    
    %rr = tp(bin==i2,:); %
    mu{level+1,1}{1,il} = tp(bin==i2,:); %points   %sub bottleneck, indexing slow
    %gg = (bit*2-1)*tmplen+tc;
    mu{level+1,1}{2,il} = (bit*2-1)*tmplen+tc; %center
    %mutree{level+1,1}{3,il} = []; %branches, now replaces corners after branch
    [l,h] = bounds(mu{level+1,1}{1,il});
    mu{level+1,1}{3,il} = [l;h]; %occupancy corners
    mu = leafcheck(mu,opt,bm,level+1,il);
end
mu{level,1}{1,layer} = []; %size(mutree(level).layer(layer).points,1);
mu{level,1}{3,layer} = bi; %branch listing
%mutree{level,1}{4,layer} = []; %clear parent corners - now using 3rd cell and replaced by branches
end
%}

function [rootsort] = rootsplit(rootcoords,cubel,points,pad)
[d,rootix] = pdist2(rootcoords,points,'euclidean','Smallest',8); %get distances to root centers
c=d<=cubel/2+pad; c(1,:)=1; %index distances in range, force nearest to always positive to ensure inclusion
rootix(~c)=0; %filter to only indices of valid roots
c = any(rootix,2); %vector for throwing out most of rootix which is unused (1, sometimes 2 rows)
rootix = rootix(c,:); %fastest reassign and trim?

rootsort = cell(1,size(rootcoords,1));
for i=1:size(rootcoords,1)
    rootsort{i} = points(any(rootix==i,1),:); % indexing all the pts pretty slow
end
end

function [roots,cubel] = cubify(points,init)
if isempty(init)
    [l,h] = bounds(points);
else
    [l,h] = bounds(init);
end
l=round(l-1,-1); h=round(h+1,-1); %bounds from points/initial bounds
%[l;h]
len = round(h-l,-1); %length of each axis
cen = (h+l)/2; %center of space
%cen*2

split = 2;
%cubifying does seem slightly faster overall, not sure how much or what context is most important
%need some metric to force breaking lowest axis length into smaller increments to ensure more roots
cubel = min(len)/split; ch = cubel/split; %cube full and half length for mesh compaction
len = ceil(len/cubel)*cubel; %length of axes in cube segmentations
cn = len/cubel; %number of cubes per axis
%cubel = len/cn

%grid not quite centered, starts at min and rounds past max values
%x4 = l(1)+ch:cubel:len(1)-ch;
%x2 = cen(1)-cubel*(cn(1)/2-0.5):cubel:cen(1)+cubel*(cn(1)/2-0.5);
%x3 = linspace(cen(1)-cubel*(cn(1)/2-0.5),cen(1)+cubel*(cn(1)/2-0.5),cn(1));
xi = (-cn(1)/2+0.5:cn(1)/2)*cubel+cen(1);
%y4 = l(2)+ch:cubel:len(2)-ch;
%y2 = cen(2)-cubel*(cn(2)/2-0.5):cubel:cen(2)+cubel*(cn(2)/2-0.5);
%y3 = linspace(cen(2)-cubel*(cn(2)/2-0.5),cen(2)+cubel*(cn(2)/2-0.5),cn(2));
yi = (-cn(2)/2+0.5:cn(2)/2)*cubel+cen(2);
%z4 = l(3)+ch:cubel:len(3)-ch; %sometimes goes zero length
%z2 = cen(3)-cubel*(cn(3)/2-0.5):cubel:cen(3)+cubel*(cn(3)/2-0.5);
%z3 = linspace(cen(3)-cubel*(cn(3)/2-0.5),cen(3)+cubel*(cn(3)/2-0.5),cn(3));
zi = (-cn(3)/2+0.5:cn(3)/2)*cubel+cen(3); %works reliably for n==1

[x,y,z] = meshgrid(xi,yi,zi); %grid of centers for root cubes

roots = [x(:),y(:),z(:)];
%{
roots = single(zeros(numel(x),3)); %reshape might be more elegant but hard to figure
for i=1:numel(x)
    roots(i,:) = [x(i),y(i),z(i)]; %list of coordinate centers for root cubes
end
%}

end

%{
function root = rsplit(root,opt)
%root
%opt = root.properties

cmp = root.points>root.center;
for i=1:8
    mask = bitget(i,1:3); %get bitmask comparator to filter where points are to a child branch
    %c = ismember(cmp,mask,'rows'); %pretty slow
    c=~any(bsxfun(@minus,cmp,mask),2); %find 000 matches, >10x faster
    
    root.branches(i).type = 'leaf';
    root.branches(i).length = root.length/2;
    root.branches(i).center = (mask*2-1)*root.length/2+root.center;
    root.branches(i).points = root.points(c,:);
    root.branches(i).branches = [];
    
    if size(root.branches(i).points,1)>opt.leafmax && root.branches(i).length>(opt.leaflen*2)
        root.branches(i) = rsplit(root.branches(i),opt);
    end
end

root.points = []; root.type = 'branch';
end
%}