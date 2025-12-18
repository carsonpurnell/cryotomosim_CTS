function [memdat] = gen_mem_atom(sz,pix,param)
% [memdat] = gen_mem_atom(sz,pix,param)
%
%
%
% [vol,~,~,splitvol] = helper_atoms2vol(pix,memdat.atoms,sz*pix); % review output example
arguments
    sz
    pix
    param.num = 2:6 % 1 still empty, alphashape failing
    param.frac = -1
    param.memsz = 1
    param.prior = []
end
box = sz*pix;

if numel(param.num)>2
    param.num = param.num(randi(numel(param.num))); % random member of set
elseif numel(param.num)==2
    param.num = sort(param.num); % sort to make sure no negatives
    param.num=randi(param.num(end)-param.num(1)+1)+param.num(1)-1; % random number between two integers
end

% individual mem definitions (mostly garbage still)
mdict(1) = struct('class','vesicle','thick',27,'thickvar',6,'size',1,'sphericity',0.9,'bare',0.2,'protfrac',1);
%mdict(2) = struct('class','er','thick',13,'thickvar',4,'size',0.6,'sphericity',0.2,'memprot',0.2);
%mdict(3) = struct('class','membrane','thick',32,'thickvar',4,'size',2,'sphericity',0.1,'memprot',0.5);
%mdict(4) = struct('class','mito','thick',35,'thickvar',3,'size',3,'sphericity',0.8,'memprot',0.0);
% make a fixed dict and use a selector function to grab the target ones?

% derived vars
% also have an auto calc fallback for not using frac? 10% per vesicle?
%if numel(num)>1, num=randi(num(end)-num(1)+1)+num(1)-1; end % target range calc
% make an easier vector expansion randi selector? x:y randi(numel(num)) sort of deal?
if param.frac<0, param.frac = min(sqrt(param.num/10)/2,1); end % fallback computed fraction of vol
param.seeds = round(param.num/param.frac)+0; 
% number of seeds needed for given membrane number and coverage ratio

%param = struct('num',num,'frac',frac,'memsz',memsz,'seeds',seeds);

[minit,blobtable] = voronoiblobcells(box,param,mdict,param.prior);
% problem: single membrane always in the middle, needs to be random

[atoms,memcell,normcell] = blob2mem(minit,blobtable,mdict,pix,param.prior);

%memdat = struct('atoms',atoms,'memcell',memcell,'normcell',normcell); %get split annoyingly
memdat.atoms = atoms;
memdat.memcell = memcell; % need to store alphashapes for in/out of vesicle placements
memdat.normcell = normcell;
memdat.table = blobtable;
end


%% internal functions
function [minit,blobtable] = voronoiblobcells(box,memparam,mdict,prior)

%memparam = struct('num',num,'frac',frac,'memsz',memsz,'seeds',seeds);
num = memparam.num;
frac = memparam.frac;
memsz = memparam.memsz;
seeds = memparam.seeds;

% set up initial volume and points
pad = max(box)/20*0+50;
isosz = ones(1,3)*max(box)+pad*2;
n = round(prod(isosz/110)/1); %approx 100-150 is reliable, higher starts to have voidless blobs
field = rand(n,3).*isosz-pad; 
% add a second larger box at lower density to feather out Z edges?
for i=1:3 % rejection loop to eliminate points outside the box to ensure isotropy
    r = field(:,i)>(box(i)+pad);
    field(r,:) = [];
end
kdt = KDTreeSearcher(prior);
[~,d] = knnsearch(kdt,field,'K',1,'SortIndices',0);
prox = d>150;
field = field(prox,:);
%plot3p(prior,'.'); hold on; plot3p(field,'.'); axis equal

% using a cube and using boxsize later stretches mems horizontally instead of vertically
% nontrivial to get isotropic distribution - rejection sampling only method? or randtess?
% might get most of the way there by making a cube & discarding above Z?

cen = rand(seeds,3).*[1,1,0.5].*box+[0,0,.25].*box; %centralize seeds a bit for less z clipping
classindex = randi(numel(mdict),seeds,1); % class as index to the mdict, not as class name
szmult = [mdict(classindex).size]'; % probably need to reintroduce some variability
szmult = szmult+szmult.*(rand(seeds,1)-rand(seeds,1))/2; % 0.5 to 1.5 spread?
class = {mdict(classindex).class};
bare = [mdict(classindex).bare]';
protfrac = [mdict(classindex).protfrac]';
blobtable = table(cen,classindex,class',szmult,bare,protfrac,...
    'VariableNames',{'centroid','classindex','class','szmult','bare','protfrac'});%,'vol'});

% can't subindex into table to change it, multiple flags per me difficult. go to struct?

% prepartitioning cells
iters = 2;%round(1*1/frac+sqrt(num)/2); % relaxation iters, probably not doing much anymore, set to 1-3?
%[cen,pf] = voronoirelax(field,cen,iters,[seed{:,1}]'); % was using weight, trying to obsolete
[blobtable.centroid,pf] = voronoirelax(field,blobtable.centroid,iters,blobtable.szmult);
pf(:,5) = 0;

% seed growing
% needs functionalized, allow for alternate methods that all feed into the same atomic meshing system
iters = 4;
minit = cell(1,seeds); v = zeros(1,seeds);

for i=1:seeds
    % randomly make 50% bare for test
    %if rand<.5, blobtable{i,3} = {'bare'}; end
    %protfrac = 
    for j=1:iters
        tmpdist = 35*memsz*blobtable.szmult(i); % 35 arbitrary, might need to change
        msel = blobtable.classindex(i);
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
    end
    ixp = pf(:,5)==i;
    minit{i} = pf(ixp,1:3);
    tmp = alphaShape(minit{i}); % no points with single membrane?
    v(i) = volume(tmp); % measure volume of local blob
end
thresh = 0.1; % to prevent excessively tiny membranes
blobtable.vol = v';
runs = 1:numel(minit);
runs = runs(v>mean(v*thresh));
runs = runs(1:min(num,numel(runs)));
% prune down to only relevant blobs
minit = minit(runs);
blobtable = blobtable(runs,:);
% do run selection here and prune down to only those blobs, leave next section for hull/mesh only
end


function [atoms,memcell,normcell] = blob2mem(minit,blobtable,mdict,pix,prior)

atoms = struct;%('vesicle',[],'er',[]);
for i=1:numel(mdict)
    atoms.(mdict(i).class) = [];
end
normcell = cell(1,numel(minit)); memcell = normcell;
for i=1:numel(minit)
    msel = blobtable.classindex(i);
    id = mdict(msel).class;
    thick = mdict(msel).thick+(rand-rand)*mdict(msel).thickvar;
    
    if numel(minit)>1
        qq = vertcat(minit{[1:i-1,1+i:end]}); % scrape all other points (faster if minit itself pruned first)
        [d] = pdist2(qq,minit{i},'euclidean','smallest',1); % detect pts in cell close to pts of other cells
        cellpts = (d>(thick*1.9+40+numel(minit)*2)); 
        %remove pts too close to other cells - not great, common edge clipping
        % use a larger distance for membranes marked for proteins?
        % need larger retreat with weighted cells since things are more squished
    else
        cellpts = 1:size(minit{i},1);
    end
    % instead alpha the whole cell and remove pts within distance of dense mesh? coverage is better
    cellpts = minit{i}(cellpts,:); %sometimes empty and fails alphashape
    
    kdt = KDTreeSearcher(prior);
    [~,d] = knnsearch(kdt,cellpts,'K',1,'SortIndices',0);
    prox = d>50;
    cellpts = cellpts(prox,:);
    
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
    
    %initshape{i} = sh1;
    [tmp,head,tail,shell,mesh] = shape2mem(sh1,thick,pix/1.5);
    % currently very wiggly, quite possibly too wiggly
    % denser mesh to reduce the wiggle? or smiter iters in already very round mems?
    
    atoms.(id) = [atoms.(id);tmp];
    
    pmesh = mesh(1:round(end/2.5),:); % preprune the mesh to speed vecnorms & clean weird facets
    memcell{i} = pmesh; % core points mesh for placing proteins
    normcell{i} = shapenorm(pmesh,sh1); % calculate normal vectors for mesh points
    % these vectors are generally correct, but a minority are significantly off normal
    % use a slower reliable method? or do more triangles and pool them?
end
f = fieldnames(atoms);
for i=1:numel(f)
    if isempty(atoms.(f{i}))
        atoms = rmfield(atoms,f{i});
    end
end
end 

function [atoms,head,tail,shell,mesh] = shape2mem(shape,thick,atomfrac)
mesh = randtess(0.5,shape,'s'); % might need raised (higher resolution?) if holes prevent memvec computing
vec = randn(size(mesh)); vec = 0.9*thick*vec./vecnorm(vec,2,2);
shell = alphaShape(mesh+vec,30+thick*2); % hopefully works across pixel/membrane sizes

tail = randtess(0.014/atomfrac,shell,'v'); % need better reference ratios
head = randtess(15/atomfrac,shell,'s'); % head domain layers across shell

vec = randn(size(head)); % random displacement directions for head density
spd = (rand(size(vec,1),1)-rand(size(vec,1),1))*(thick*0.2); % triangular random displacement distances
vec = vec./vecnorm(vec,2,2).*spd; % displacement vectors
head=head+vec; atoms = [head;tail];
atoms(:,4) = 5.5/1 *atomfrac; % magnitude of pseudoatoms
end

function [cen,pts] = voronoirelax(pts,cen,iters,weight)
if nargin<4, weight = ones(1,size(cen,1)); end
if iters==0||size(cen,1)<2; iters=1; skip=1; else skip=0; end

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

p0 = pts(ix(1,:),:);
p1 = pts(ix(2,:),:); % hardpoints for plane/normal vector
p2 = pts(ix(3,:),:);
p3 = pts(ix(4,:),:);

ntmp = cross(p1-p0,p2-p0)+cross(p2-p0,p3-p0)+cross(p3-p0,p1-p0); % 3 facet normals
%ntmp = ntmp+cross(p4-p6,p5-p6)+cross(p1-p5,p3-p5);
nvecs = ntmp./vecnorm(ntmp,2,2); % normalize vectors to unit length
in = inShape(sh,pts+nvecs*2); %test if pts inside shape
nvecs(in,:) = nvecs(in,:)*-1; % invert pts inside the shape to make all point outward from surface
end