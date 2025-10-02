% volume parameters
pix = 12;
sz = [400,300,80];
box = sz*pix;

% global membrane parameters
frac = -1;%.3; % fallback ~.1 per membrane
num = 1:12; % default 4-ish?    % change randomizer to be in the blob retension step instead of seed count?
memsz = 1;

%sphmult = 0.9; % make per-mem, and not used by voronoi for iters?
%szvar = 1; % multiplier to membrane size variation, should be membrane data, need to remove from voronoi step
thresh = 0.1; % multiplier to volume threshold (frac of mean) for inclusion
%thickvar = 3; % absolute variation in angstroms % should be part of the membrane data, not global

% individual mem params
mdict(1) = struct('class','vesicle','thick',28,'thickvar',6,'size',1,'sphericity',0.9);
%mdict(2) = struct('class','er','thick',14,'thickvar',4,'size',0.6,'sphericity',0.2);
%mdict(3) = struct('class','membrane','thick',32,'thickvar',4,'size',2,'sphericity',0.1);
%mdict(4) = struct('class','mito','thick',35,'thickvar',3,'size',3,'sphericity',0.8);
%cdict = {mdict(:).class}; % lookup from class name to index
% make a fixed dict and use a selector function to grab the target ones?
%mdat = {28,'vesicle',4;14,'er',3}; % thickness/label/thickvar in angstroms

% derived vars
% also have an auto calc fallback for not using frac? 10% per vesicle?
if numel(num)>1, num=randi(num(end)-num(1)+1)+num(1)-1; end % target range calc
% make an easier vector expansion randi selector? x:y randi(numel(num)) sort of deal?
if frac<0, frac = min(sqrt(num/10)/2,1); end % fallback computed fraction of vol
seeds = round(num/frac)+0; % number of seeds needed for given membrane number and coverage ratio
%szmult = (frac+10*(frac/num)-1/num)*memsz; % probably also needs to be a per-membrane instead of global var

%% notes blarg
% need to refactor again into segments/subfunctions: pts/cells, blob growth, pruning, final hull, vectors
% final check to nudge top/bottom points inside box? or let it stay more random & weird?


% alternate method: some sort of force field collision generator to obviate normal collision&voronoi
% generate points in region around ecnter, apply force to expand like balloon (ideally non-spherical)
% eventual size of blob determined by starting number of points - need adhesive force too for containing
% makes it easier to have separated blobs, rather than current full occupancy with missing blobs
% how to efficiently do forces? iterative pdist2, move individual poitns or entire blobs?
% have each blob record longest radial distance, pdist2 between blobs to precheck possible clashes first?
% then iterate between each potential pair to move individual points?

% maybe simpler: instead of trying to have an empty bubble w/ surface tension, fill space instead?
% repulsion should automatically fill up the area without needing anchor points
% normal-looking shapes might emerge from limited movement (time/resistance) from nonuniform initial placement
% not clear how plasma mem sheets might be generated from a volume effect, ER may also be difficult

% initial point density needs to be asymmetric to create organic shapes during force field expansion
% how to expand reliably? apply force from blob centroid to shell, then adherent force?
% adherent force might be an easy vector towards the average of N nearest points, reliable smooth surface

% ER tubules: strong attractive force to some sort of thin skeleton framework? how to generate skeleton?
% create some noise field and derive skeleton from it? smooth & dilate to create bias for skeletonization?
% mitos: nest membrane, then shrink outer layer slowly, contract inner instead of inflating?

% current cell partitioning also not easy to create nested blobs, not inherent in approach.
% force field version could do nested a bit easier, just spawn more starting seeds after inflation


% need to make initial centroids more uniformly distributed at higher sphericity to avoid squish wedges
% assign memtypes based on measured characteristics instead?
% create seeds and examine initial bins for shape, then generate with them as needed? can't prefilter though

% another option: rely on voronoi relaxation. sph value is ratio of spheres, fraction that relaxes
% other pts do not relax by iterative centroid movement in partition, so more likely to stay squish
% may get too spherical, need few iters(and movement cap) or all will be centered too well

% yet another: compute dist from edge of bin (as bwdist) rather than from centroid - more irregular
% pdist surface points (unfortunately probably very slow) to find distance to closest surface

% instead of iter growth scalar, use fractional size of the cell to limit continued expansion - might avoid
% overlaps

% how to make surfaces irregular again? flat squish works, but everything is very smooth
% also no ER-like small tubules or large planes yet - st the start to impose shape or just before finish?
% how best to create mem planes? centroid off the box, collect points within narrow range as radius & smooth?

%% set up initial volume and points
pad = max(box)/20*0+50;
isosz = ones(1,3)*max(box)+pad*2;
n = round(prod(isosz/100)/1); %approx 100-150 is reliable, higher starts to have voidless blobs
field = rand(n,3).*isosz-pad; 
% add a second larger box at lower density to feather out Z edges?
for i=1:3 % rejection loop to eliminate points outside the box to ensure isotropy
    r = field(:,i)>(box(i)+pad);
    field(r,:) = [];
end
% using a cube and using boxsize later stretches mems horizontally instead of vertically
% nontrivial to get isotropic distribution - rejection sampling only method? or randtess?
% might get most of the way there by making a cube & discarding above Z?

cen = rand(seeds,3).*[1,1,0.5].*box+[0,0,.25].*box; %centralize seeds a bit for less z clipping
classindex = randi(numel(mdict),seeds,1); % class as index to the mdict, not as class name
szmult = [mdict(classindex).size]'; % probably need to reintroduce some variability
szmult = szmult+szmult.*(rand(seeds,1)-rand(seeds,1))/2; % 0.5 to 1.5 spread?
class = {mdict(classindex).class};
blobtable = table(cen,classindex,class',szmult,...
    'VariableNames',{'centroid','classindex','class','szmult'});%,'vol'});

%{
seed = cell(seeds,2); %1 for size, 2 for centroid?
%basesize = 1;
for i=1:seeds
    seed{i,1} = 1+(rand-1/(1+szvar))*szvar; %increasing variation window, 0+
    %basesize+(rand-rand)*basesize*1*szvar; % problem: sometimes negative if scaled above 1
    % weight from this rather than separate val - derive from type sizing?
    seed{i,2} = cen(i,:);
end
%}
%{
wcen = 5; % central tendancy of weighting around 1 - use triangle instead?
weight = 1+(rand(size(cen,1),1)-1/2)/wcen; % weighting increases overlap rate
irreg = 0; % span width of randomization
weight = 1+rand(size(cen,1),1)*irreg;
% rationally get actual weights from lookup of lipid class? big for mito, small vesicles/ER?
%weight = randi(6,size(seedcen,1),1); % weights increase rate of initial overlapping
%weight = 1;
%}

%% prepartitioning
% need to functionalize voronoi relax
%{
[d] = pdist2(seedcen,field,'euclidean');%,'Smallest',1); 
d = d./weight;
[~,ix] = min(d);
% weighting causing overlaps in dense areas, need to relax more to avoid interleaving pts
% centroid might even be outside of small cells sometimes
pf = field; pf(:,4) = ix; pf(:,5) = 0; % 4th column initial closest seed, 5th column actual inclusion

% need at least one relaxation iter to reduce overlap rate
seedcen2 = seedcen*0;
for i=1:seeds
    seedcen2(i,:) = mean(pf(pf(:,4)==i,1:3));
end
[d] = pdist2(seedcen2,field,'euclidean');%,'Smallest',1); 
d = d./weight;
[~,ix] = min(d);
% weighting causing overlaps in dense areas, need to relax more to avoid interleaving pts
% centroid might even be outside of small cells sometimes
pf = field; pf(:,4) = ix; pf(:,5) = 0;
%}

iters = 2;%round(1*1/frac+sqrt(num)/2); % relaxation iters, probably not doing much anymore, set to 1-3?
%[cen,pf] = voronoirelax(field,cen,iters,[seed{:,1}]'); % was using weight, trying to obsolete
[blobtable.centroid,pf] = voronoirelax(field,blobtable.centroid,iters,blobtable.szmult);
pf(:,5) = 0;
%blobtable.centroid = cen;
%{
for i=1:size(seed,1)
    seed{i,2} = cen(i,:);
end
%}
%%
for i=1:seeds
    ix = pf(:,4)==i;
    plot3p(pf(ix,1:3),'.'); hold on;
end
axis equal


%% seed growing
%basedist = 35;  % current sphmult badly implemented, large size differences & relatively low impact
iters = 4;
for i=1:seeds
    for j=1:iters
        %tmpdist = basedist*seed{i,1};
        tmpdist = 35*memsz*blobtable.szmult(i); 
        msel = blobtable.classindex(i);
        % 35 is arbitrary scalar value to get 'good'-looking vesicles
        ix = pf(:,4)==i; ix = find(ix);
        subsel = pf(ix,1:3);
        
        % simplify growing procedure? use only from centroid for speed and simplicity
        % prorate based on # identified points, distance change in new centroid, centroid of new points
        % high sphericity low tolerance for centroid,, high for points?
        
        %search everything against the centroid to avoid weird interleaving?
        %[d] = pdist2(subsel,seed{i,2},'euclidean'); % centroid proximity catch
        [d] = pdist2(subsel,blobtable.centroid(i,:),'euclidean');
        if j>1,qq=mdict(msel).sphericity; else qq=1; end
        cendist = qq*tmpdist*15+j*30*mdict(msel).sphericity;
        ix2 = d<cendist; ix2 = ix(ix2);
        pf(ix2,5) = i;
        
        ix3 = pf(:,5)==i;
        
        ix = pf(:,4)==i; ix = find(ix);
        subsel = pf(ix,1:3);
        [d] = pdist2(subsel,pf(ix3,1:3),'euclidean'); % radiating proximity catch - OOM with few mems
        % penalize distance based on centroid to curb cornering?
        tmpix = min(d,[],2);
        proxdist = ((1-mdict(msel).sphericity)*tmpdist*1.1-j*1.5);
        tmpix2 = tmpix< proxdist;
        ix2 = ix(tmpix2);
        
        pf(ix2,5) = i;
        ix = pf(:,5)==i; ix = find(ix);
        
        blobtable.centroid(i,:) = mean(pf(ix,1:3),1);
        %seed{i,2} = mean(pf(ix,1:3),1);
    end
end
%%
minit = cell(1,seeds); v = zeros(1,seeds);
for i=1:seeds
    ixp = pf(:,5)==i;
    minit{i} = pf(ixp,1:3);
    plot3p(minit{i},'.'); hold on;
    
    %sh = alphaShape(minit{i}); 
    %sh.Alpha = criticalAlpha(sh,'one-region')*2;
    v(i) = volume(alphaShape(minit{i}));
end
axis equal
%alternate strat: prepartition, then iterate over each and have a moving centroid rather than static
% during prepartition, remove pts within X distance of other points? 
% no, do just before shelling based on thickness of the current membrane to take all avail space
% use some sort of imposed structure to force membranes to distort? too spherical currently
% shrinking area by removing outside certain planes could compress them + leave protein room
% raising initial seeds far past reasonable targets could also create squish zones
% earmark fake, space-wasting seeds with smaller initial size to squish into larger real ones?

% still need to merge centroids that are too close together to avoid micro domains

% probably need to keep all particles in the same array with 4th ID column to avoid overlap pruning loops

% also need to prune small-volume blobs to avoid micelles, also helps speed (merge into filler group?)
% easier to distance prune to non-overlap if there's a dedicated background split too
%% expanding seeds into vesicle blobs
%{
for i=1:iters
    for j=1:targs
        % initial hard distance collector - maybe does nothing useful?
        % maybe do a pre-partition of all points to the closest seed to avoid overlaps? allows per-seed size
        [d] = pdist2(rest,seed{j}(1,:),'euclidean');%,'Smallest',1);
        distparam = 60+j*4;
        %hits = unique(ix(d<distparam));
        hits = find(d<distparam);
        seed{j} = [seed{j};rest(hits,:)];
        rest(hits,:) = [];
        %secondary NN collector
        distparam = 50;%-sqrt(size(seed{j},1));
        [d,ix] = pdist2(rest,seed{j},'euclidean','Smallest',4); % pretty slow - pre-filter pts might help?
        % need some sort of limit parameter - distance from centroid/seed or sphericity or mem size
        hits = unique(ix(d<distparam));
        seed{j} = [seed{j};rest(hits,:)];
        rest(hits,:) = [];
        
        sh = alphaShape(seed{i}); sh.Alpha = criticalAlpha(sh,'one-region')*2;
        v(j) = volume(sh);
    end
end
%}
% can leave unclustered points inside of an existing cluster at lower density
%t = vertcat(seed{:});

%% plot field
%plot3p(rest.*1,'.'); axis equal; hold on; 
%plot3p(t,'.');
%{
for i=1:targs
    plot3p(seed{i}.*1,'.'); hold on;
end
axis equal
%}

%% mesh out membranes
%{
mtypes = [14,27]; % thicknesses for ER,vesicle/PM    what is mito thickness?
atoms = struct('vesicle',[],'er',[]);
mdat = {28,'vesicle';14,'er'};
thickvar = 3; % should be part of the membrane data, not global
%}
%n = targs;%round(numel(minit)/1);
runs = 1:numel(minit);
runs = runs(v>mean(v*thresh));
runs = runs(1:min(num,numel(runs)));

atoms = struct('vesicle',[],'er',[]);
normcell = cell(1,numel(minit));
for i=runs
    %msel = randi(2);
    msel = blobtable.classindex(i);
    %mempick = mdict(msel);
    id = mdict(msel).class;
    thick = mdict(msel).thick+(rand-rand)*mdict(msel).thickvar;
    %thick = mdat{msel,1}+(rand-rand)*mdat{msel,3};
    %id = mdat{msel,2};
    
    qq = vertcat(minit{[1:i-1,1+i:end]}); % scrape all other points (faster if minit itself pruned first)
    [d] = pdist2(qq,minit{i},'euclidean','smallest',1); % detect pts in cell close to pts of other cells
    %[d2] = pdist2(qq,minit{i}); %slowest
    %[d,ix] = pdist2(minit{i},qq,'euclidean','smallest',1);
    %[d2] = pdist2(minit{i},qq);
    cellpts = (d>(thick*1.6+45)); %remove pts too close to other cells - not great, common edge clipping
    % use a larger distance for membranes marked for proteins?
    % need larger retreat with weighted cells since things are more squished
    
    % instead alpha the whole cell and remove pts within distance of dense mesh? coverage is better
    cellpts = minit{i}(cellpts,:); %sometimes empty and fails alphashape
    % need to instead alphashape the whole cell, get a dense surface mesh, remove all pts within thick*~1.1
    % points are not uniform so cells might already be separated by several A, can leave more leeway
    if size(cellpts,1)<4, continue; end
    
    % apparently need to rework shape2shell/shell2dens, alphashape not behaving well and retaining surface
    % simultaneous creation of both, so regeneration the alphashape doesn't break surface detail between?
    
    % too many iterations here, need simplification into only necessary components
    %sh = alphaShape(cellpts,1000); %1000 doesn't break at lower
    %[~,tmp1] = boundaryFacets(sh);
    %sh.Alpha = criticalAlpha(sh,'one-region')*10 % 500-1500ish?
    
    sh1 = alphaShape(cellpts,250);
    tmp2 = randtess(0.5,sh1,'s'); % sometimes has wacky infills from incomplete internal tesselation
    
    %{
    %sh2 = alphaShape(cellpts,40);  %400
    %sh2.Alpha = criticalAlpha(sh,'one-region')*5 % approx 465?
    %[~,tmp] = boundaryFacets(sh);
    
    %tmp = randtess(.5,sh,'s');
    sh2 = alphaShape(tmp2,600); %600
    % interblobs are introduced here, despite changing alpha & no issues in tesselation
    %sh2.Alpha = criticalAlpha(sh2,'one-region')*1.5; % 300-800
    
    tmp3 = randtess(1.0,sh2,'s'); % raising improves collision, but slower normal vectors
    sh3 = alphaShape(tmp3); % minor bottleneck
    sh3.Alpha = criticalAlpha(sh3,'one-region')*1.2; % still required to avoid interblobs
    
    initshape{i} = sh3; % not enough memory to store both
    %shape{i} = randtess(4,sh,'s');
    shell = shape2shell(sh3,thick*.9); % new bottleneck
    [tmp,head,tail] = shell2pts(shell,pix/4,thick*.2);
    %}
    
    initshape{i} = sh1;
    [tmp,head,tail,shell,mesh] = shape2mem(sh1,thick,pix/4);
    % currently very wiggly, quite possibly too wiggly
    % denser mesh to reduce the wiggle? or smiter iters in already very round mems?
    
    atoms.(id) = [atoms.(id);tmp];
    % second round of trimming edges from actual membrane shape (or only round, drop initial prune?)
    
    normcell{i} = shapenorm(mesh,initshape{i}); % calculate normal vectors for mesh points
    % bottleneck, need to reduce point counts
end

[vol,~,~,splitvol] = helper_atoms2vol(pix,atoms,box);
sliceViewer(vol);

%plot3p(mesh,'.'); hold on; plot3p(mesh+normcell{end}*10,'.');

%% normal vector calcs
% needs functionalizing
%normcell = cell(1,numel(initshape));
%n=4;
%{
for i=runs
    %{
    %skelpts = [skelpts;shell{i}.Points];
    tmpskel = initshape{i}.Points;
    % this pdist is now the bottleneck
    [~,ix] = pdist2(tmpskel,tmpskel,'squaredeuclidean','Smallest',n); % nearest n pts for speed in fewer iters
    %  se 71.3
    % fse needs version 2023
    
    p1 = tmpskel(ix(2,:),:); % hardpoints for plane/normal vector
    p2 = tmpskel(ix(3,:),:);
    p3 = tmpskel(ix(4,:),:);
    ntmp = cross(p1-p3,p2-p3);
    normt = ntmp./vecnorm(ntmp,2,2);
    qw = ntmp./sqrt(sum(ntmp.^2,2)); % identical within precision
    qwe = qw-normt;
    % vectors still randomly oriented, need to create inverted copy and inShape test them to get only out vecs
    in = inShape(initshape{i},initshape{i}.Points+normt*5); %test if pts inside shape
    % why is nothing in shape?
    invecs = normt(in,:);
    normt(in,:) = normt(in,:)*-1; % invert all inward-facing normal vectors
    %}
    normt = shapenorm(initshape{i});
    normcell{i} = normt;
end
%}
%plot3p(initshape{1}.Points,'.'); hold on; plot3p(initshape{1}.Points+normcell{1}*5,'.');

%% normal from hardpoint plane (pts 2-4)
% cross function should be vectorized, should be able to compute all normals at once after prepping nx3 arrays
%{
pts = tmp(2:4,:);
dif1 = pts(1,:)-pts(3,:);
dif2 = pts(2,:)-pts(3,:);
nvec = cross(dif1,dif2);
nvec = nvec/norm(nvec); % appears to be accurate normal vector, need to check with double norm
planevec = cross(dif1,nvec);
planevec = planevec/norm(planevec); % single planar vector, visibly perpendicular
plot3(tmp(:,1),tmp(:,2),tmp(:,3),'.'); hold on;
plot3(pts(:,1),pts(:,2),pts(:,3),'.'); hold on;
qq = [tmp(1,:);tmp(1,:)+1*nvec;tmp(1,:)-1*nvec];
plot3(qq(:,1),qq(:,2),qq(:,3),'o'); axis equal
ww = [tmp(1,:);tmp(1,:)+5*planevec;tmp(1,:)-5*planevec];
plot3(ww(:,1),ww(:,2),ww(:,3),'-o'); axis equal
%}

%% get plane and normal vectors - super annoying for fitted normal plane (no hardpoints)
%{
x = tmp(:,1); y = tmp(:,2); zz = tmp(:,3);
xx = min(x):max(x); yy = min(y):max(y);
[xx,yy] = meshgrid(xx,yy);
%[xx yy] = meshgrid(x,y);
%zz = C(1)*xx+C(2)*yy + C(3) + 2*randn(size(xx));
plot3(x(:),y(:),zz(:),'.')

XC = mean(tmp,1);
Y=tmp-XC;
[~,~,V2]=svd(Y,0);
C = V2(:,end); % same as other V

P=[mean(x),mean(y),mean(zz)];
[U,S,V]=svd([x-P(1),y-P(2),zz-P(3)],0);
N=-1/V(end,end)*V(:,end);
A=N(1); B=N(2); C=-P*N;
C = [A,B,C]; % plane coefficients

O = ones(length(x),1);
C2 = [x y O]\zz; C2=C2'; % also plane coeffifients, nonidentical but close

zzft = C(1)*xx+C(2)*yy + C(3);
hold on;
nv = V(:,3)'; %nv = nv([2,1,3]); % not xy intervetd
qq = [tmp(1,:);tmp(1,:)+2*nv;tmp(1,:)-2*nv];
plot3(qq(:,1),qq(:,2),qq(:,3),'o');
surf(xx,yy,zzft,'edgecolor','none')
hold on;
zz2 = C2(1)*xx+C2(2)*yy + C2(3);
%surf(xx,yy,zz2,'edgecolor','none')
grid on
%}

%% garbled tangent calculations
%{
n=9; % closest pts
mdin = KDTreeSearcher(ptsin); [ixin] = knnsearch(mdin,skelpts,'K',n)'; 
mdout = KDTreeSearcher(ptsout); [ixout] = knnsearch(mdout,skelpts,'K',n)';
normvec = zeros(size(skelpts)); % might need to loop over individual membranes - especially if combo memprots

for i=1:size(skelpts,1)
    skelcen = skelpts(i,:);
    long = outcen-incen; long = long/norm(long);
    under = skelcen-incen; under = under/norm(under);
    over = outcen-skelcen; over = over/norm(over); %average over and under, then average with long to refine?
    refined = ((over+under)/2+long)/2; refined = refined/norm(refined);
    %refined = (over+under)/2; refined = refined/norm(refined);
    %refined = (refined+long)/2; refined = refined/norm(refined);
    
    vx = skelpts(i,1); vy = skelpts(i,2); vz = skelpts(i,3); %recover subscript data for the current point
    norm4d(vx,vy,vz,[1,2,3]) = refined; %write normals to 4d storage array
end

%% 
%[x,y,z] = ind2sub(size(skel),find(skel==1)); skelpts = [x,y,z];

%n = 9; %nearest n voxels on inner and outer surfaces to calculate vectors
%mskel = KDTreeSearcher(skelpts); %[ixself] = knnsearch(skelpts,skelpts,'K',5);
mdin = KDTreeSearcher(ptsin); [ixin] = knnsearch(mdin,skelpts,'K',n)'; 
mdout = KDTreeSearcher(ptsout); [ixout] = knnsearch(mdout,skelpts,'K',n)';

norm4d = zeros(size(skel,1),size(skel,2),size(skel,3),3);
%normmat = zeros(size(idx,1),3); %normmat2=normmat1;

%vectorization to calculate vector targets for each point of skel
%q = ptsin(ixin,:); q2 = reshape(q,[],3,n); win = sum(q2,3)/n;
%q = ptsout(ixout,:); q2 = reshape(q,[],3,n); wout = sum(q2,3)/n;
%need to refactor this vectorization to be more streamlined, probably doable in many fewer steps
q = ptsin(ixin(:),:); q2 = reshape(q',3,n,[]); win = permute(mean(q2,2),[1,3,2])';
q = ptsout(ixout(:),:); q2 = reshape(q',3,n,[]); wout = permute(mean(q2,2),[1,3,2])';

%{
%giant pile of testing jank from figuring out how to trick reshape into working
%might need to change how points are initially generated to make it easier to reshape to something useful
%instead get the xyz coords individually to make reshaping easier?
%size(ixout)
%ixout(1:n*2)
%gg = ixout(:,1:2)
%reshape(gg',[],1)
%q(1:n*2,:)
%reshape(q(1:n*2,:)',3,n,[])
%q2(1:2,:,1:4)
%q2(
%size(ixout) %numel*n, reshape to linear order to make it easy?
%ixout(1:2,:)
%reshape(ixout(1:2,:)',[],1)
%size(q), size(q2)
%ss = ixout(1:2,:)
%skelpts(1,:)
%size(q)
%q(1:n*2,:)
%reshape(q(1:n*2,:)',3,n,[])
%xx = q(1:n*2,:)' %q transposed to make order work the same as ss - breaks again if multiple rows
%aa = reshape(xx,3,n,[])
%cc = permute(mean(aa,2),[1,3,2])
%dd = ptsout(ss,:)' %ss already works as column 
%dd is running through indices in column order/linearly, q runs through the indices in row order
%q2(:,:,1)
%ptsout(ixout(1,:),:)'
%ptsout(ss,:)'
%zz = reshape(dd,[],n,3)
%ff = permute(dd,[3,2,1])
%q2(1,:,:) %oriented wrong, sum across dims going between coordinate
%size(wout), size(ixin), size(skelpts)
%wout(1:10,:)
%win(1:10,:)
%}
for i=1:size(skelpts,1)
    skelcen = skelpts(i,:);
    incen = win(i,:);
    %{
    try
        incen = win(i,:); %index out of bounds error line
    catch
        %size(q), size(q2), size(win)
    end
    %}
    outcen = wout(i,:);
    long = outcen-incen; long = long/norm(long);
    under = skelcen-incen; under = under/norm(under);
    over = outcen-skelcen; over = over/norm(over); %average over and under, then average with long to refine?
    refined = ((over+under)/2+long)/2; refined = refined/norm(refined);
    %refined = (over+under)/2; refined = refined/norm(refined);
    %refined = (refined+long)/2; refined = refined/norm(refined);
    
    vx = skelpts(i,1); vy = skelpts(i,2); vz = skelpts(i,3); %recover subscript data for the current point
    norm4d(vx,vy,vz,[1,2,3]) = refined; %write normals to 4d storage array
end
%}


%% internal functs
function smpts = smiter(pts,iters,nearn)
smpts = zeros(size(pts));
for j=1:iters
    [~,ix] = pdist2(pts,pts,'euclidean','Smallest',nearn);
    ix = ix(2:end,:);
    for i=1:size(ix,2)
        mm = mean(pts(ix(:,i),:));
        smpts(i,:) = mm;
    end
    pts = unique(smpts,'rows');
end
end

function [shell,meshpts] = shape2shell(shape,thick)
meshpts = randtess(1,shape,'s'); % might need raised (higher resolution?) if holes prevent memvec computing
vec = randn(size(meshpts)); vec = thick*vec./vecnorm(vec,2,2);
shell = alphaShape(meshpts+vec,30+thick); % hopefully works across pixel/membrane sizes
%shell.Alpha = criticalAlpha(shell,'one-region')*1.2; % slow, so using hard value instead
end

function [pts,head,tail] = shell2pts(shell,atomfrac,surfvar)
tail = randtess(0.028/atomfrac,shell,'v'); % need better reference ratios
head = randtess(24.0/atomfrac,shell,'s'); % head domain layers across shell

vec = randn(size(head)); % random displacement directions
spd = (rand(size(vec,1),1)-rand(size(vec,1),1))*surfvar; % triangular random displacement distances
vec = vec./vecnorm(vec,2,2).*spd; % displacement vectors
head=head+vec; pts = [head;tail];
pts(:,4) = 6.0/2 *atomfrac; % magnitude of pseudoatoms
end

function [atoms,head,tail,shell,mesh] = shape2mem(shape,thick,atomfrac)
mesh = randtess(0.2,shape,'s'); % might need raised (higher resolution?) if holes prevent memvec computing
vec = randn(size(mesh)); vec = 0.9*thick*vec./vecnorm(vec,2,2);
shell = alphaShape(mesh+vec,30+thick*2); % hopefully works across pixel/membrane sizes

tail = randtess(0.013/atomfrac,shell,'v'); % need better reference ratios
head = randtess(12.5/atomfrac,shell,'s'); % head domain layers across shell

vec = randn(size(head)); % random displacement directions for head density
spd = (rand(size(vec,1),1)-rand(size(vec,1),1))*(thick*0.2); % triangular random displacement distances
vec = vec./vecnorm(vec,2,2).*spd; % displacement vectors
head=head+vec; atoms = [head;tail];
atoms(:,4) = 6.0/1 *atomfrac; % magnitude of pseudoatoms
end

function [cen,pts] = voronoirelax(pts,cen,iters,weight)
if nargin<4, weight = ones(1,size(cen,1)); end
if iters==0; iters=1; skip=1; else skip=0; end
for j=1:iters
    [d] = pdist2(cen,pts(:,1:3),'euclidean'); % distances for each point to each centroid
    d = d./weight; % weight distances by centroid
    [~,ix] = min(d); % identify cell membership by lowest weighted distance
    pts(:,4) = ix;
    if skip<1
        for i=1:size(cen,1)
            cen(i,:) = mean(pts(pts(:,4)==i,1:3));
        end
    end
end
end

function nvecs = shapenorm(pts,sh) % compute normal vectors from shape core mesh
[~,ix] = pdist2(pts,pts,'squaredeuclidean','Smallest',4); % nearest n pts for speed in fewer iters

p1 = pts(ix(2,:),:); % hardpoints for plane/normal vector
p2 = pts(ix(3,:),:);
p3 = pts(ix(4,:),:);
ntmp = cross(p1-p3,p2-p3); % vectors normal to local facet triangle
nvecs = ntmp./vecnorm(ntmp,2,2); % normalize vectors to unit length
in = inShape(sh,pts+nvecs*5); %test if pts inside shape
nvecs(in,:) = nvecs(in,:)*-1; % invert pts inside the shape to make all point outward from surface
end