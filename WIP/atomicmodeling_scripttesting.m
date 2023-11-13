%% load input structures as atomic data
%rng(2);
pix = 8; clear particles;
input = {'actin__6t1y_13x2.pdb',...
    'tric__tric__6nra-open_7lum-closed.group.pdb',...
    'ribo__ribo__4ug0_4v6x.group.pdb','MT__6o2tx3.pdb','MT__6o2tx2.pdb'};%,...
    %'ATPS.membrane.complex.cif'};%,'a5fil.cif','a7tjz.cif'};
    %input = {'CaMK2a__5u6y.pdb','CaMK2a__5u6y.pdb','CaMK2a__5u6y.pdb','CaMK2a_3soa.pdb'};
    input = {'actin__6t1y_13x2.pdb','MT__6o2tx3.pdb','ribosome__4ug0.pdb'};
tic
%need to streamline atomic symbol to Z converter, and link into a Z to scatterval dictionary.
%extend it to work for pdb2vol as well, and the other older cts_model components.
for i=numel(input):-1:1 %backwards loop for very slightly better performance
    particles(i) = helper_pdb2dat(input{i},pix,2,1,0);
    particles(i).perim = {vertcat(particles(i).perim{:})};
    particles(i).adat = {vertcat(particles(i).adat{:})};
end
layers{1} = particles; fprintf('loaded %i structure files  ',numel(input));
toc
%{
for i=1:numel(particles)
    for j=1:numel(particles(i).atomid)
        %particles(i).atomint{j} = atomdict(particles(i).atomid{j},'sc');
        %com = mean(particles(i).atomcoords{j},1); %need radius from geometric, not mass center
        %particles(i).atomcoords{j} = particles(i).atomcoords{j}-com;
        %size(particles(i).atomcoords{j})
        %size([0,0,0])
        %particles(i).radius{j} = max(pdist2(particles(i).atomcoords{j},single([0,0,0])));
        alphat = alphaShape(double(particles(i).adat{j}(:,1:3)),12); %surprisingly slow
        [~,p] = boundaryFacets(alphat);
        n = size(particles(i).atomcoords{j},1);
        ix = randperm(n); ix = ix(1:round(n/400));
        pi = particles(i).atomcoords{j}(ix,1:3);
        p = single([p;pi]); %need to add back 1-3% or so of points to prevent inside placements
        particles(i).perim{j} = unique(p,'rows');
        %clear com alphat %p pi
    end
end

%might need to make atom id vector a row vector for memory write reasons
%need atomdict function that accepts vector of atomic symbols and returns vector of Z/e- values
%}

%% functionalized model gen part
%rng(1);
boxsize = pix*[500,400,50]*1;
[splitin.carbon,dyn] = gen_carbon(boxsize); % atomic carbon grid generator
memnum = 0;
tic; [splitin.lipid,kdcell,shapecell,dx.lipid,dyn] = modelmem(memnum,dyn,boxsize); toc;

con = helper_atomcon(boxsize,pix); % pseudonatural ice border (wavy flat, no curvature)
dyn{1}(dyn{2}:dyn{2}+size(con,1)-1,:) = con; dyn{2}=dyn{2}+size(con,1)-1;

n = 500;
%profile on
%tic; [split,dyn,mu] = fn_modelgen(layers,boxsize,n,splitin,dx,dyn); toc
tic; [split,dyn,mu] = helper_randfill_atom(layers,boxsize,n,splitin,dx,dyn); toc
%profile viewer

%% function for vol, atlas, and split generation + water solvation
outpix = pix;
[vol,solv,atlas,splitvol] = helper_atoms2vol(outpix,split,boxsize);
sliceViewer(vol*1+solv);
%{
WriteMRC(vol+solv,outpix,'solv_n1.mrc');
%WriteMRC(atlas,outpix,'radtestingmix_8a_atlas.mrc');
%}

%{
%{
%% atomic vesicle gen
%currently just a hamfisted first-pass in the modelgen. separate implementation needed? need better outputs
%mem proteins with their own modelgen? either generated with the vesicle initially, or separate pregenerator
%placed vesicle map filler: premade KDT for overlap testing, track and proxfilt only the proteins
%prefiller: run after each vesicle generation to fill, and make sure to retry placements a lot. easier
%per-type membranes (thickness exclusions or whatever) for dissimilar vesicle/membrane types
ves = 0; memhull = 0;
vesarg = {250+randi(200),[],rand*0.2+0.8, 24+randi(8)};
if ves>0
    %[splitin,memhull,dyn] = fn_modgenmembrane(ves,layers);
    lipid(1).name = 'lipid'; lipid(1).flags = 'ves';
    tic
    fprintf('generating membranes  ')
    for i=1:ves
        %[pts,perim] = vesgen_sphere(200+randi(300),18+randi(5)); %old deprec spherical version
        [pts,perim] = gen_mem(250+randi(200),[],rand*0.2+0.8, 24+randi(8)); %need fewer more intense points
        %need core shape for anchoring membrane proteins as well
        %also need hull of all of them for inside/outside checks
        %pts(:,4) = pts(:,4)/4; %288 init, 170/195 at 1/2 pts, 125 at 1/4
        lipid(1).perim{1,i} = perim;
        lipid(1).adat{1,i} = pts;
        lipid(1).modelname{i} = append('vesicle');%,string(i));
    end
    layers{2} = layers{1};
    layers{1} = lipid;
    toc
end

% need mem and memprot model generator. only way to properly have layers without superbloat.
% fill out sh level a bit more for tighter alphashape that won't blob into others as easily - also more smooth

% 470 at /5 atoms perimeter
% 102 at /50 atoms perimeter
% huge difference, need to carry forward the dynpts
%}
%{
%% functionalized model gen part
rng(1);
boxsize = pix*[400,300,50];
[splitin.carbon,dyn] = gen_carbon(boxsize); % atomic carbon grid generator
memnum = 8;
tic; [splitin.lipid,kdcell,shapecell,dx.lipid,dyn] = modelmem(memnum,dyn,boxsize); toc;
%splitin3.carbon = splitin.carbon; 
%splitin.lipid = splitin2.lipid; %super dumb temporary hackjob

n = 500;
%n = [50,3000];
%splitin.border = borderpts;
tic; [split] = fn_modelgen(layers,boxsize,n,splitin,dx,dyn); toc

%% function for vol, atlas, and split generation + water solvation
[vol,solv,atlas,splitvol] = helper_atoms2vol(pix,split,boxsize);
sliceViewer(vol+solv);
%WriteMRC(vol+solv,pix,'atomictest_fastgen1.mrc');
%}
%{
%% randomly add to the points and concatenate them into a list
boxsize = pix*[200,300,50];
%modelpoints = pts+boxsize/2; modelid = atomid;
%modelpoints =  modelid = 0; modelid2=single(modelid); 
dynpts = single([-100 -100 -100]); %dynid = single(0);
dynpts = single(zeros(0,3));
%modeltree = KDTreeSearcher(modelpoints');
rng(5);
tol = 2; %tolerance for overlap testing
count.s = 0; count.f = 0;
n = 1000; %ixcat = 2; %iscat 1 to erase the initial point
ixincat = 1; %index 1 to overwrite the initial preallocation point, 2 preserves it
split = cell(1,numel(particles)+1); split{1} = single([0,0,0,0]); %split{2} = 0;
%zz = 0; %ac = [];
tl = tic; 
for i=1:n
    if rem(i,n/20)==0; fprintf('%i,',i); end
    
    which=randi(numel(particles));
    sel = particles(which);
    sub = randi(numel(sel.adat));
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
        tpts = transformPointsForward(tform,sel.adat{sub}(:,1:3))+loc;
        tpts = [tpts,sel.adat{sub}(:,4)]; %#ok<AGROW>
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
        split{which+1} = [split{which+1};tpts]; %splitvol add
        %split{which+1,2} = [split{which+1,2},sel.atomint{sub}];
        
        count.s=count.s+1;
    else
        count.f=count.f+1;
    end
end
toc(tl); clear tl
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
    %splitvol(:,:,:,i) = fnpt2vol(pix,split{i,1}(:,1:3),split{i,1}(:,4),boxsize,[0,0,0]);
    splitvol(:,:,:,i) = helper_pt2vol(pix,split{i,1},boxsize,[0,0,0]);
end
[~,atlas] = max(splitvol,[],4); 
em = sum(splitvol,4); %small difference with em in a few scattered points
%atlas2 = min(atlas,proj); %zero out empty voxels to separate 0 from the first split
toc
sliceViewer(em);
%}
%% solvation op2: flat density+variance window
%subtract volume from each voxel, estimate waters per voxel from remaining volume with random portion
%works well, need better values and improved randomization. should probably be built into helper_pt2vol
%now is built in. kind of.

%% solvation testing with explicit water particles - way too slow for even small models
%just too slow to prune millions of points.
%{
rng(7)
tic
allatoms = vertcat(split{2:end,1});
distfrac = 0.5;
[solv] = gen_solvate(allatoms(:,1:3),round(boxsize/1),distfrac,tol);
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

%}
%}

%% internal functions
% need a single wrapper that runs carbon, membrane, memprots, allprots sequentially and carries dyn
% nesting is heavy jank and makes subcalls confusing, can't reuse code as much

function con = internal_atomcon(box,pix,n,sc)
if nargin<3
    n = 4+pix^1.5;
    sc = 2400;
end
sz = [max(box),max(box)]; 
dl = 10;
w = 1;
ptsb = internal_gen_atomborder(sz,n/2,sc*1,4)-[0,0,dl*randi(6)*w];
ptst = internal_gen_atomborder(sz,n/2,sc*1,4)+[0,0,dl*randi(6)*w+box(3)];
pts = zeros(0,3);
for i=1:4
    pts = [pts;ptsb-[0,0,dl*i]]; pts = [pts;ptst+[0,0,dl*i]];
end
con = pts;
end
function pts = internal_gen_atomborder(sz,n,sc,sep)
%n = 2.5; % noise magnitude
%sc = 500; % scale of Z noise
%sep = 3;
pad = 20; %padding - scale by input size maybe? prune afterward?
sd = max(sz)+pad*2;
[X,Y] = ndgrid(1:sep:sd,1:sep:sd);
i = min(X-1,sd-X+1); j = min(Y-1,sd-Y+1);
H = exp(-.5*(i.^2+j.^2)/n^2);
Z = real(ifft2(H.*fft2(randn(size(X))))); % 0-centered, approximately normal

pts = [X(:),Y(:),Z(:)*sc];
n = size(pts,1); ix = randperm(n); ix = ix(1:round(n/10));
pts = pts(ix,:);
pts(:,1:2) = pts(:,1:2)-pad;
%size(pts)
end

function [splitin,memhull,dyn] = fn_modgenmembrane(memnum,vesarg,layers)

%ves = number of vesicles as input
[kdcell,shapecell] = modelmem(memnum,vesarg);
%layers = the input particle layers to add to generated membranes

end

function [err,loc,tform,ovcheck,ix] = anyloc(boxsize,tperim,dyn,retry,tol,mu)
for r=1:retry
    loc = rand(1,3).*boxsize; tform = randomAffine3d('rotation',[0 360]); %random placement
    ovcheck = transformPointsForward(tform,tperim)+loc; %transform test points
    %err2 = proxtest(dyn{1}(1:dyn{2}-1,:),ovcheck,tol); %prune and test atom collision
    [err,ix] = mu_search(mu,ovcheck,tol,'short',0); %slightly faster!
    err = any(err>0);
    %if err2~=err, fprintf('%i,%i,\n',err,err2); end
    if err==0, break; end
end
end

function [pts,kdcell,shapecell,dx,dyn] = modelmem(memnum,dyn,boxsize)
dyn = {dyn,size(dyn,1)}; %convert to dyncell
kdcell = []; shapecell = [];

tol = 2; %tolerance for overlap testing
retry = 5; %retry attempts per iteration
count.s = 0; count.f = 0;
lipid{1} = zeros(0,4); lipid{2} = 1;
for i=1:memnum % simplified loop to add vesicles
    [tpts,tperim] = gen_mem(250+randi(200),[],rand*0.2+0.8, 24+randi(8));
    
    [err,loc,tform,ovcheck] = anyloc(boxsize,tperim,dyn,retry,tol);
    %{
    for r=1:retry    
        loc = rand(1,3).*boxsize; tform = randomAffine3d('rotation',[0 360]); %random placement
        ovcheck = transformPointsForward(tform,tperim)+loc; %transform test points
        err = proxtest(dyn{1}(1:dyn{2}-1,:),ovcheck,tol); %prune and test atom collision
        if err==0, break; end
    end
    %}
    
    if err==0
        tpts(:,1:3) = transformPointsForward(tform,tpts(:,1:3))+loc;
        
        [dyn] = dyncell(ovcheck,dyn);
        [lipid] = dyncell(tpts,lipid);
        
        %[dyn{1},dyn{2}] = dyncat(dyn{1},dyn{2},ovcheck);
        
        %[pts,dx] = dynsplit(tpts,pts,dx,splitname);
        %[split.(splitname),dx.(splitname)] = dyncat(split.(splitname),dx.(splitname),tpts);
        %{
        % % inlined dyncat code, split assignments % %
        tdx = dx.(splitname); %MUCH faster than hard cat, ~7x.
        l = size(tpts,1); e = tdx+l-1;
        if e>size(split.(splitname),1)
            split.(splitname)(tdx:(tdx+l)*4,:) = 0;
        end
        split.(splitname)(tdx:e,:) = tpts; dx.(splitname) = tdx+l;
        % % inlined dyncat code % %
        %}
        count.s=count.s+1;
    else
        count.f=count.f+1;
    end
    
end
pts = lipid{1}(1:lipid{2}-1,:); dx = lipid{2};
fprintf('vesicles: placed %i, failed %i  ',count.s,count.f)
end

function [dyn,ix] = dyncat(dyn,ix,pts)
l = size(pts,1); e = ix+l-1;
if e>size(dyn,1)
    dyn(ix:(ix+l)*2,:) = 0;
end
dyn(ix:e,:) = pts; 
ix = ix+l;
end
function [split,dx] = dynsplit(tpts,split,dx,splitname) %slower than inlined a bit
tdx = dx.(splitname);
l = size(tpts,1); e = tdx+l-1;
if e>size(split.(splitname),1)
    split.(splitname)(tdx:(tdx+l)*2,:) = 0;
end
split.(splitname)(tdx:e,:) = tpts; dx.(splitname) = tdx+l;
end


function [split,dyn,mu] = fn_modelgen(layers,boxsize,niter,split,dx,dyn)
if nargin<6
    dyn{1} = single(zeros(0,3)); dyn{2} = 0;
end
leaf = 1e3;
mu = mu_build(dyn{1},'leafmax',leaf,'maxdepth',2);
if nargin<4
    split = struct; %ixincat = 1; %dynpts = single(zeros(0,3));
else
    fn = fieldnames(split);
    for i=1:numel(fn) %add split into dynpts
        s = size(split.(fn{i})(:,1:3),1);
        %ix = randi(s,round(s/50),1); ix = unique(ix);
        %tmp = split.(fn{i})(ix,1:3);
        %l = size(tmp,1); %dynpts(end+1:end+l,:) = tmp;
        %dynpts = [dynpts;tmp];
        %dynpts = [dynpts;split.(fn{i})(:,1:3)];
    end
end
%ixincat = size(dynpts,1)+1; %where to start the indexing
%dynfn = dynpts; dynfnix = ixincat;
%dyn = {dynfn,dynfnix};

% available location mapping - probably not going to be faster due to needing to prune full list after
% each successful placement
%instead, after octree tech, do a basic check to see how many points are nearby
%need the dynamic octree!
%{
gridmaptol = 12;
n = prod(boxsize)/((gridmaptol*1)^3); %number of map points
locgrid = rand(n,3).*boxsize; %generate map points
gtree = KDTreeSearcher(dynpts);
[~,d] = rangesearch(gtree,locgrid,gridmaptol,'SortIndices',0);
p = zeros(1,numel(d));
for i=1:numel(d)
    if isempty(d{i}), p(i) = 1; end
end
locgrid = locgrid(logical(p),:);

%[lgridvol] = helper_atoms2vol(6,locgrid,boxsize);
%sliceViewer(lgridvol);
%}

%d = [d{:}]; %if any(d<gridmaptol), err=1; end
%locgrid = locgrid(d>gridmaptol,:);
%size(locgrid)
%plot3(locgrid(:,1),locgrid(:,2),locgrid(:,3),'.'); axis equal
%[lgridvol] = helper_atoms2vol(6,locgrid,boxsize);
%sliceViewer(lgridvol);

%tmp = fieldnames(split);
%dynpts = split;
tol = 2; %tolerance for overlap testing
retry = 4; %retry attempts per iteration
count.s = 0; count.f = 0;
%index 1 to overwrite the initial preallocation point, 2 preserves it
%split = cell(1,numel(particles)+0); %split{1} = zeros(0,4); %single([0,0,0,0]);
for i=1:numel(layers)
namelist = [layers{i}.modelname]; %slower than cell, but more consistent
for j=1:numel(namelist)
    if ~isfield(split,namelist{j})
        split.(namelist{j}) = zeros(0,4); %initialize split models of target ids
    end
    if ~isfield(dx,namelist{j})
        dx.(namelist{j}) = size(split.(namelist{j}),1)+1;
    end
end
end

for lc = 1:numel(layers)
    particles = layers{lc};
    n = niter(lc);
    
for i=1:n
    if rem(i,n/20)==0; fprintf('%i,',i); end
    
    which=randi(numel(particles));
    sel = particles(which); sub = randi(numel(sel.adat));
    tpts = sel.adat{sub}(:,1:3);
    [err,loc,tform,ovcheck,muix] = anyloc(boxsize,tpts,dyn,retry,tol,mu); % just as fast
    %{
    for r=1:retry    
        loc = rand(1,3).*boxsize; tform = randomAffine3d('rotation',[0 360]); %random placement
        
        ovcheck = transformPointsForward(tform,sel.perim{sub})+loc; %transform test points
        err = proxtest(dyn{1}(1:dyn{2}-1,:),ovcheck,tol); %prune and test atom collision
        %err = proxtest(dynpts(1:ixincat-1,:),ovcheck,tol); %prune and test atom collision
        %need to replace either with mutable quadtree or short-circuit kdtree, it's ~80% of runtime
        if err==0, break; end
    end
    %}
    
    if err==0
        tpts = sel.adat{sub};
        tpts(:,1:3) = transformPointsForward(tform,tpts(:,1:3))+loc;
        
        
        [dyn] = dyncell(ovcheck,dyn); %.15
        mu = mu_build(ovcheck,muix,mu,'leafmax',leaf,'maxdepth',2);
        [split,dx] = dynsplit(tpts,split,dx,sel.modelname{sub}); %3.6
        %{
        %[dynfn,dynfnix] = fcndyn(ovcheck,dynfn,dynfnix); % insignificantly slower than inlined
        %[dyn{1},dyn{2}] = dyncat(dyn{1},dyn{2},ovcheck); %6s
        %modnm = sel.modelname{sub};
        %[split.(modnm),dx.(modnm)] = dyncat(split.(modnm),dx.(modnm),tpts); %46 - crazy slow :(
        
        %{
        % % inlined dyncat code, dynpts % %
        l = size(ovcheck,1); e = ixincat+l-1;
        if e>size(dynpts,1)
            dynpts(ixincat:(size(dynpts,1)+l)*3,:) = 0;
        end
        dynpts(ixincat:e,:) = ovcheck; ixincat = ixincat+l;
        %}
        % % inlined dyncat code, split assignments % %
        %{
        tdx = dx.(sel.modelname{sub}); %MUCH faster than hard cat, ~7x.
        l = size(tpts,1); e = tdx+l-1;
        if e>size(split.(sel.modelname{sub}),1)
            split.(sel.modelname{sub})(tdx:(tdx+l)*4,:) = 0;
        end
        split.(sel.modelname{sub})(tdx:e,:) = tpts; dx.(sel.modelname{sub}) = tdx+l;
        % % inlined dyncat code % %
        %}
        
        %split{which} = [split{which};tpts]; %add to splitvol - cell version
        %split.(sel.modelname{sub}) = [split.(sel.modelname{sub});tpts]; %struct a bit slower :(
        %}
        count.s=count.s+1;
    else
        count.f=count.f+1;
    end
end
if lc==1
    %tic; sh=alphaShape(double(dynpts),12); toc; %plot(sh); drawnow;
end
fprintf('  placed %i, failed %i \n',count.s,count.f);
%all(dyn{1}==dynpts) %all(dynfn==dynpts)

end
sn = fieldnames(split); %trimming trailing zeros from split arrays to prevent atom2vol weirdness
for i=1:numel(sn)
    tdx = size(split.(sn{i}),1); %backstop for when the object was preexisting so there's no dx
    if isfield(dx,sn{i})
        tdx = dx.(sn{i});
    end
    split.(sn{i})(tdx:end,:) = [];
end

end

function err = proxtest(c,pts,tol)
l = min(pts,[],1)-tol; h = max(pts,[],1)+tol; %low and high bounds per dimension
ix = c>l & c<h; % compare all points against the prospective box
ix = all(ix,2); % filter to index of pts inside the box
if ~any(ix), ix=[]; end % check for early end if no points in the box
%ix = find(ix>0); %bottleneck - just too many points. mutable octree should be faster overall
err=0; %with n=100 exhaustive is only slightly slower than kdtree search, but progressive slowdown
if ~isempty(ix) %this thing is taking SO VERY LONG, need more pre-optimization
    buck = 100;%round( size(c,1)/7650 ); %very rough, is probably not linear scale
    % probably needs some sort of depth-based metric, not a flat one depth = log2 (n/leaf)
    %ot = OcTree(c(ix,:),'binCapacity',buck); %significantly slower than kdt build
    %leaf = 1e3;
    %mutree = octcubetree(c(ix,:),'leafmax',leaf); %slightly faster than kd building
    %err = mutreetest(mutree,pts); %WAYYY slower than knn search
    %err = any(err);
    %mutreec = octcubetree_cell(c(ix,:),'leafmax',leaf); %slightly faster than kd building
    %err = mutreetest_cell(mutreec,pts);
    
    modeltree = KDTreeSearcher(c(ix,:),'Bucketsize',buck); %67 with 1K %32 with 10K, 18 100K
    [~,d] = rangesearch(modeltree,pts,tol,'SortIndices',0); %?? 1K,11.4 10K, 85 100K
    d = [d{:}]; if any(d<tol), err=1; end %test if any points closer than tol
end
end


function [dyn] = dyncell(addpts,dyn)
    l = size(addpts,1); e = dyn{2}+l-1;
    if e>size(dyn{1},1)
        dyn{1}(dyn{2}:(size(dyn{1},1)+l)*3,:) = 0;
    end
    dyn{1}(dyn{2}:e,:) = addpts;
    dyn{2} = dyn{2}+l;
end
function [dynfn,ix] = fcndyn(ovcheck,dynfn,ix) % insignificantly slower than inline version
    l = size(ovcheck,1); e = ix+l-1;
    if e>size(dynfn,1)
        dynfn(ix:(size(dynfn,1)+l)*3,:) = 0;
    end
    dynfn(ix:e,:) = ovcheck;
    ix = ix+l;
end

function [pts,perim] = vesgen_sphere(r,thick)
radi = r;
rado=radi+thick;
%radi = (rand*700+150)/pix; %randomly generate inner radius of vesicle (need better range)
%rado = radi+(14+randi(14))/pix; %get outer radius from inner, should be constant something (7-9nm-ish?)
%reduced outer radius distance for pearson, skew makes it wider
%offset = round(rado+20); %centroid offset to prevent negative values
%still not sure how to do the radius and what the radial density curve should look like

w = (rado-radi)/1.5; %deviation of the membrane distribution
sf = [(rado^2)/(radi^2),(radi^2)/(rado^2)]/2; %factor to correct for excess inner density
%correction factor seems a bit off. little too much inner density still?

%fill space between radii with tons of points
%ptnum = round(radi*5*(pix^3)*pi^2); %need to actually calculate volume of shell
shellvol = pi*(rado^3-radi^3)*0.16; %volume of shell in pseudoatoms
%ptnum = round( 0.2*shellvol*1^3 )*2; %convert to angstroms, scale to some arbitrary working density
%frac = [ptnum,ptnum*sf(2),ptnum*sf(1)]; %get fractions of the total to distribute between inner and outer
rti = round(shellvol*sf(2)); %rto = ptnum-rti; %partition density between inner and outer radii
rto = round(shellvol*sf(1)); %slightly more points to balance out?
ptnum = rti+rto;
%ptrad = rand(ptnum,1)*(rado-radi)+radi; %uniform - flat monolayer
switch 1
    case 1 %mirrored pearson - relatively hard inner and outer edges
        ptrad = [pearsrnd(radi,w,0.7,3,rti,1);pearsrnd(rado,w,-0.7,3,rto,1)];
    case 2 %mirrored gamma - a bit narrower, more edge smoothing
        ptrad = radi+[betarnd(3.0,6,rti,1);betarnd(6,3.0,rto,1)]*(rado-radi)*3.0;
end
%pearson is very slow, calls beta to call gamma which takes most of the time
%need to reformulate the math so that density is hard-bound between ri/ro in angstroms

ptaz = rand(ptnum,1)*pi*2; %random circular azimuth angles
%ptel = rand(1,ptnum)*pi*2; %causes asymmetry, polar density accumulation
ptel = asin(2*rand(ptnum,1)-1); %random elevation angles, corrected for polar density accumulation

[x,y,z] = sph2cart(ptaz,ptel,ptrad); %convert spherical coords to cartesian coords
pts = [x,y,z];
%ves = 0;
n = size(pts,1); ix = randi(n,round(n/50),1);
%perimix = randperm(n); permix = perimix(1:round(n/400)); perim = pts(perimix,:);
perim = pts(ix,:);
in = ones(size(pts,1),1)*2.95;
pts = [pts,in];

end

%{
function [vol,solv] = helper_atoms2vol_i(pix,pts,sz,offset)
if nargin<4, offset=[0,0,0]; end
if nargin<3, sz = max(pts,[],1)+pix; end
%if size(pts,2)<4, pts(:,end+1)=1; end %intensity==1 if not given by 4th column
%need rough estimate of average volume for organic atoms
%very approximately 1.8a radii
%eventually might do individual vdw radii individually
avol = 4/3*pi*(1.8^3); %eyeballed volume of the average organic atom
h20 = 2.1; %overrounded number for water magnitude
pts(:,1:3) = round((pts(:,1:3)-offset)/pix+0.5);
emsz = floor(sz/pix); vol = zeros(emsz);
solv = (rand(emsz)-0.5)*1*pix^2+(pix^3);
for i=1:3
    ix = pts(:,i) < emsz(i) & pts(:,i) > 1; %get points inside the box
    pts = pts(ix,:); %drop points outside the box
end
for i=1:size(pts,1)
    x=pts(i,1); y=pts(i,2); z=pts(i,3); mag = pts(i,4); %fetch data per atom
    vol(x,y,z) = vol(x,y,z)+mag;
    solv(x,y,z) = solv(x,y,z)-avol;
end
solv = max(solv,0)/32*h20;
end
%}
function vol = ifcn_solv(pix,pts,sz,offset)
if nargin<4, offset=[0,0,0]; end
if nargin<3, sz = max(pts,[],1)+pix; end
%if size(pts,2)<4, pts(:,end+1)=1; end %intensity==1 if not given by 4th column
%need rough estimate of average volume for organic atoms
%very approximately 1.8a radii
%eventually might do individual vdw radii individually
avol = 4/3*pi*(1.8^3); %eyeballed volume of the average organic atom
h20 = 2.1; %overrounded number for water magnitude
pts(:,1:3) = round((pts(:,1:3)-offset)/pix+0.5);
emsz = floor(sz/pix); vol = (rand(emsz)-0.5)*0*pix^2+(pix^3);
for i=1:3
    ix = pts(:,i) < emsz(i) & pts(:,i) > 1; %get points inside the box
    pts = pts(ix,:); %drop points outside the box
end
for i=1:size(pts,1)
    x=pts(i,1); y=pts(i,2); z=pts(i,3); %mag = pts(i,4); %fetch data per atom
    vol(x,y,z) = vol(x,y,z)-avol*(rand*.2+0.9);
end
vol = max(vol,0)/32*h20;
end
function [tmp] = gen_solvate(modpts,sz,distfrac,tol)
h20vol = 35;
atomfrac = 1;
loops = 1; %appears to make no difference, may at full solvation size
buck = round( size(modpts,1)/450 ); %very rough early optimized value
modeltree = KDTreeSearcher(modpts,'Bucketsize',buck); %make model for all the existing points
%distfrac = 0.5; %would need to scale by pixel size, .5 appears fine for 10a
%so very very slow. gets too rough with denser pseudoatoms though
%shift density towards flat value, mask via atomic solv volgen?
waters = prod(sz)/h20vol*distfrac*atomfrac/loops; %number of predicted waters to fill the box
%pix = 10; flatd = (pix^3)/h20vol;
%make a very slightly larger box for solvation placement to avoid edge effects - 5-10A
solv = zeros(0,3);
for i=1:loops
tmp = rand(round(waters),3,'single').*sz; %~1.2-1.4G per 100 million coords
%pad = 20; spacing = 3/distfrac;
%[x,y,z] = meshgrid(1-pad:spacing:sz(1)/4+pad,1-pad:spacing:sz(2)/4+pad,1-pad:spacing:sz(3)/4+pad);
%solv = [x(:),y(:),z(:)]; %weird grids due to interpolation/binning shenanigans
[idx] = rangesearch(modeltree,tmp,tol/2,'SortIndices',0); %find waters near extant atoms - super slow
ol = ~cellfun(@isempty,idx); %waters that had a proximity fail
tmp = tmp(~ol,:); %keep only nontouching waters
solv = [solv;tmp];
end

end
function [v,ix] = dyncatold(ix,v,b)
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

%{
function vol = helper_pt2vol(pix,pts,sz,offset)
if nargin<4, offset=[0,0,0]; end
if nargin<3, sz = max(pts,[],1)+pix; end
if size(pts,2)<4, pts(:,end+1)=1; end %intensity==1 if not given by 4th column
pts(:,1:3) = round((pts(:,1:3)-offset)/pix+0.5);
emsz = floor(sz/pix); vol = zeros(emsz);
for i=1:3
    ix = pts(:,i) < emsz(i) & pts(:,i) > 1; %get points inside the box
    pts = pts(ix,:); %drop points outside the box
end
for i=1:size(pts,1)
    x=pts(i,1); y=pts(i,2); z=pts(i,3); mag = pts(i,4); %fetch data per atom
    vol(x,y,z) = vol(x,y,z)+mag;
end
end
%}
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

%{
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

function err = sskdtrange(kdt,pts,tol)
err = 0;
for i=1:size(pts,1)
    [~,d] = rangesearch(kdt,pts(i,:),tol,'SortIndices',0); %just too slow with validations
    if any([d{:}]<tol); err = 1; break; end
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
%}

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