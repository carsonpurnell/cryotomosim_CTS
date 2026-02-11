%% placing memprots with new atomic membrane structure
% now generally functional - do still need vesicle/cytosol placement, doable in normal placer
pix = 8;
targ = {'ATPS.membrane.complex.mat'};
%targ = {'ATPS.membrane.mat'};
%targ = {'GABAar.membrane.complex.mat'};
targ = {'GABAar.membrane.complex.mat','ATPS__flip.6j5i.membrane.cif'...%}; %flip.mat not flipped?
    'COVID19_spike.membrane.complex.pdb','ETC1_huge__6zkq.membrane.cif','kchannel__1bl8.membrane.pdb'};
targ = {'ATPS__flip.6j5i.membrane.cif'};

%targ = {'ATPS__flip.6j5i.membrane.cif','tric_6nra-open_7lum-closed.group.mat','GABAar.membrane.complex.mat'};
targ = {'ATPS__flip.6j5i.membrane.mat';'GABAar.membrane.complex.mat';'tric_6nra-open_7lum-closed.group.mat'};

pmod = param_model(pix,'layers',3,'mem',5:10);
%%
sz = [600,600,100];
%
if true
    [carbon,perim] = gen_carbon(sz*pix);
else
    carbon = []; perim = [];
end
% irregular carbon overlaps from C-shape membranes closing over the carbon edge sometimes
% add some randomization in for the z-height of the carbon

%skips if mem==1?
memdat = gen_mem_atom(sz,pix,'num',pmod.mem,'frac',1,'prior',perim);%,'memsz',1,'frac',-1); % needs carbon exclusion and input
% needs a bit more work, a few vectors (probably due to corners) are not well-oriented - denser mesh?
% alternate method - dense surface mesh of expanded membrane hull, remove inner points, get nearest?
% would need to be very dense. but could average with the near-3 result to cover most cases?

%[c] = c+helper_atoms2vol(pix,carbon,sz*pix); %sliceViewer(c);

%{
for i=1:100
[carbon,perim] = gen_carbon(sz*pix);
[c] = c+helper_atoms2vol(pix,carbon,sz*pix);
end
sliceViewer(c); % cumulative carbon projection
%}

%%
[split,dx,dyn] = helper_randfill_atom_mem(memdat,pmod,perim,sz);
% need to add mem perims to dyn at the end to feed into normal placement engine

split.carbon = carbon; % move into membrane placement instead?
dx.carbon = size(carbon,1);

[vol,solv,atlas] = helper_atoms2vol(pix,split,sz*pix);
sliceViewer(vol);

%%
%{
tmp = struct2cell(split);
atoms = vertcat(tmp{:}); % rediculously slow - reverse search order?

%tic; kdt = KDTreeSearcher(memdat.atoms.vesicle(:,1:3)); toc %11 12 15
tic; kdt2 = KDTreeSearcher(atoms(:,1:3)); toc %77 68 67
%
%tic; [ix,d] = knnsearch(kdt,atoms(:,1:3),'K',1,'SortIndices',0); toc %554 119 121
tic; [~,d2] = knnsearch(kdt2,memdat.atoms.vesicle(:,1:3),'K',1,'SortIndices',0); toc %147 201 168
%
ixf = d2>4.0; % possibly a bit high
mematoms = memdat.atoms.vesicle(ixf,:);

% KDT seems to take way too long to be worth it, even the faster method is 4 mins.
% alt 1: ad-hoc in atoms2vol (or variant for handling mems) in simulator - janky and annoying
% alt 2: per-membrane KDT for lot fewer searches at a time - after placement loop?

%
split.carbon = carbon; % after membrane pruning for a bit of speed

%split.vesicle = mematoms;
[vol,solv,atlas] = helper_atoms2vol(pix,split,sz*pix);
%mvol = helper_atoms2vol(pix,mematoms,sz*pix);
%sliceViewer(max(vol,mvol));
sliceViewer(vol);

%% internal prep stuff
split = struct; dx = struct; list = struct;
layers = pmod.layers;
for i=1:numel(layers)
namelist = [layers{i}.modelname]; %slower than cell, but more consistent
for j=1:numel(namelist)
    if ~isfield(split,namelist{j})
        split.(namelist{j}) = zeros(0,4); %initialize split models of target ids
        list.(namelist{j}) = zeros(0,3);
    end
    if ~isfield(dx,namelist{j})
        dx.(namelist{j}) = size(split.(namelist{j}),1)+1;
    end
end
end
%memat = vertcat(memdat.memcell{:});
%memat(:,4) = 3;
%split.mem = memdat.atoms.vesicle;

if isempty(perim)
    dyn{1} = single(zeros(0,3)); dyn{2} = 1;
else
    dyn{1} = perim; dyn{2} = size(perim,1)+1;
end
leaf = 1e3;
mu = mu_build(dyn{1},[0,0,0;sz*pix],'leafmax',leaf,'maxdepth',2);

init = [0,0,1]; % origin vector for rotation calculations
%sel = pmod.layers{1}(1);
tol = 2;
%% placement loop
% currently each membrane is acting as a layer, new set of extra layers would be too messy
% use existing layers, randomly (or from membrane definition?) select layer to use?
retry = 5;
count.s = 0; count.f = 0;
prog = 0; progdel = ''; % initialize starting vals for progress bar
for j=1:numel(memdat.memcell)
    prog = prog + 100/numel(memdat.memcell); %progress update block
    progstr = sprintf('progress %3.0f, membrane %i of %i', prog,j,numel(memdat.memcell));
    fprintf([progdel, progstr]);
    progdel = repmat(sprintf('\b'), 1, length(progstr));
    
memsel = j;%randi(numel(memdat)); % select membrane to place on

%memflags = memdat.table.class(j); if any(matches(memflags,'bare')); continue; end
if rand<memdat.table.bare(j); continue; end % skip if membrane doesn't hit dictionary fraction

%qq = vertcat(memdat.memcell{[1:j-1,1+j:end]});
kdt = KDTreeSearcher(vertcat(memdat.memcell{[1:j-1,1+j:end]}));

% placement attempt iterations, based on meshpts available and class fraction
iters = size(memdat.memcell{memsel},1)*.02*memdat.table.protfrac(j);

for i=1:iters 
    sel = pmod.layers{1}(randi(numel(pmod.layers{1})));
    % inner loop: random axial rotation, rotation to transmembrane vector, collision test
    if any(ismember(sel.flags,'complex'))
        sel.sumperim = vertcat(sel.perim{:});
        subsel = 0;
    else
        subsel = randi(numel(sel.id));
        sel.sumperim = sel.perim{subsel};
    end
    
    
    for v=1:retry % % % start of placement test block for building a subfunct
    % need to funct out placement testing, and allow multiple attempts per iteration
    
    % coordinates are not spatially ordered, so can be selected randomly or linearly
    meshsel = randi(size(memdat.memcell{memsel},1));
    
    memloc = memdat.memcell{memsel}(meshsel,:); % selected coordinate
    surfvec = memdat.normcell{memsel}(meshsel,:); % normal vector to surface at coordinate
    
    rotax=cross(init,surfvec); %compute the normal axis from the rotation angle
    theta = -acos( dot(init,surfvec) ); % angle between ori vec and surface
    
    spinang = rand*2*pi;
    rot1 = sel.sumperim*rotmat(init,spinang); % apply random axial rotation
    rot2 = rot1*rotmat(rotax,theta)+memloc; % apply rotation to surface vector and translate
    %diagori = init*rotmat(rotax,theta);
    
    %
    [err,muix] = mu_search(mu,rot2,tol,'short',0);
    err = any(err>0);
    if err==1; continue; end
    %
    % mu first marginally faster - might be more so with more atoms (carbon, etc)
    % mu first stays faster with increasingly large iterations/accumulated atoms
    
    if err==0
        [ix,d] = knnsearch(kdt,rot2,'K',1,'SortIndices',0); % sort false might be faster
        %if any(d<15), er2=1; else er2=0; end % hard switch since no base value for er2
        if any(d<16), err=1; else err=0; end
    end
    
    if err==0; break; end
    end % % % end of block for test placement subfunct
    
    %{
    if err==0
    [err,muix] = mu_search(mu,rot2,tol,'short',0);
    err = any(err>0);
    end
    %}
    
    % if no collision, switch to place subunits as needed after replicating rotations
    if err==0
        [dyn] = dyncell(rot2,dyn); % write to partial list for fast lookups (obsolete under mu?)
        % no muix yet, not testing
        mu = mu_build(rot2,muix,mu,'leafmax',leaf,'maxdepth',2);
        
        % make write block into a generic subfunct and propogate to other atomfills?
        if subsel==0 % if complex, write each individual submodel after transforming
            tmpix = 1:numel(sel.id);
            % for assembly, instead construct tmpix with the randomized models needed
            for u=tmpix
                tmp = sel.adat{u};
                %if er2==1; tmp(:,4)=tmp(:,4)*2; end % diag collision test against membranes
                tmp(:,1:3) = tmp(:,1:3)*rotmat(init,spinang);
                tmp(:,1:3) = tmp(:,1:3)*rotmat(rotax,theta)+memloc;
                
                [split,dx] = dynsplit(tmp,split,dx,sel.id{u});
                list.(sel.id{u})(end+1,:) = memloc;
            end
        else % if not complex, use subsel index to write only that model
            tmp = sel.adat{subsel};
            %if er2==1; tmp(:,4)=tmp(:,4)*2; end % diag collision test against membranes
            tmp(:,1:3) = tmp(:,1:3)*rotmat(init,spinang);
            tmp(:,1:3) = tmp(:,1:3)*rotmat(rotax,theta)+memloc;
            
            [split,dx] = dynsplit(tmp,split,dx,sel.id{subsel});
            list.(sel.id{subsel})(end+1,:) = memloc;
        end
        
        count.s=count.s+1;
    else
        count.f=count.f+1;
    end
end

end
%% cleanup and output stuff
fprintf('   membrane embedding complete, placed %i, failed %i \n',count.s,count.f);

% remove trailing zeros from atom registry and sparse tracker
f = fieldnames(split);
for i=1:numel(f)
    split.(f{i})(dx.(f{i}):end,:) = [];
end
dyn{1}(dyn{2}:end,:) = [];


split.carbon = carbon;
[vol,solv,atlas] = helper_atoms2vol(pix,split,sz*pix);
mvol = helper_atoms2vol(pix,memdat.atoms,sz*pix);
sliceViewer(max(vol,mvol));

%diagpts = [init;rotax;surfvec];
% plot3p(list.ATPS_head,'o'); hold on; plot3p(list.ATPS_head+memdat.normcell{1}(1:50,:)*100,'.'); % diag vecs
% plot3p(dyn{1}(1:dyn{2}-1,:),'.'); hold on; plot3p(memdat.memcell{1},'.'); % diag placements
%}

%% internal functions

function [split,dx,dyn] = int_fill_mem(memdat,carbon,perim,pmod,sz)
% internal prep stuff
split = struct; dx = struct; list = struct;
pix = pmod.pix;
layers = pmod.layers;
for i=1:numel(layers)
namelist = [layers{i}.modelname]; %slower than cell, but more consistent
for j=1:numel(namelist)
    if ~isfield(split,namelist{j})
        split.(namelist{j}) = zeros(0,4); %initialize split models of target ids
        list.(namelist{j}) = zeros(0,3);
    end
    if ~isfield(dx,namelist{j})
        dx.(namelist{j}) = size(split.(namelist{j}),1)+1;
    end
end
end
%memat = vertcat(memdat.memcell{:});
%memat(:,4) = 3;
%split.mem = memdat.atoms.vesicle;

if isempty(perim)
    dyn{1} = single(zeros(0,3)); dyn{2} = 1;
else
    dyn{1} = perim; dyn{2} = size(perim,1)+1;
end
leaf = 1e3;
mu = mu_build(dyn{1},[0,0,0;sz*pix],'leafmax',leaf,'maxdepth',2);

init = [0,0,1]; % origin vector for rotation calculations
%sel = pmod.layers{1}(1);
tol = 2;

% placement loop
% currently each membrane is acting as a layer, new set of extra layers would be too messy
% use existing layers, randomly (or from membrane definition?) select layer to use?
retry = 5;
count.s = 0; count.f = 0;
prog = 0; progdel = ''; % initialize starting vals for progress bar
for j=1:numel(memdat.memcell)
    prog = prog + 100/numel(memdat.memcell); %progress update block
    progstr = sprintf('progress %3.0f, membrane %i of %i', prog,j,numel(memdat.memcell));
    fprintf([progdel, progstr]);
    progdel = repmat(sprintf('\b'), 1, length(progstr));
    
memsel = j;%randi(numel(memdat)); % select membrane to place on

%memflags = memdat.table.class(j); if any(matches(memflags,'bare')); continue; end
if rand<memdat.table.bare(j); continue; end % skip if membrane doesn't hit dictionary fraction

%qq = vertcat(memdat.memcell{[1:j-1,1+j:end]});
kdt = KDTreeSearcher(vertcat(memdat.memcell{[1:j-1,1+j:end]}));

% placement attempt iterations, based on meshpts available and class fraction
iters = size(memdat.memcell{memsel},1)*.02*memdat.table.protfrac(j);

for i=1:iters 
    sel = pmod.layers{1}(randi(numel(pmod.layers{1})));
    % inner loop: random axial rotation, rotation to transmembrane vector, collision test
    if any(ismember(sel.flags,'complex'))
        sel.sumperim = vertcat(sel.perim{:});
        subsel = 0;
    else
        subsel = randi(numel(sel.id));
        sel.sumperim = sel.perim{subsel};
    end
    
    
    for v=1:retry % % % start of placement test block for building a subfunct
    % need to funct out placement testing, and allow multiple attempts per iteration
    
    % coordinates are not spatially ordered, so can be selected randomly or linearly
    meshsel = randi(size(memdat.memcell{memsel},1));
    
    memloc = memdat.memcell{memsel}(meshsel,:); % selected coordinate
    surfvec = memdat.normcell{memsel}(meshsel,:); % normal vector to surface at coordinate
    
    rotax=cross(init,surfvec); %compute the normal axis from the rotation angle
    theta = -acos( dot(init,surfvec) ); % angle between ori vec and surface
    
    spinang = rand*2*pi;
    rot1 = sel.sumperim*rotmat(init,spinang); % apply random axial rotation
    rot2 = rot1*rotmat(rotax,theta)+memloc; % apply rotation to surface vector and translate
    %diagori = init*rotmat(rotax,theta);
    
    %
    [err,muix] = mu_search(mu,rot2,tol,'short',0);
    err = any(err>0);
    if err==1; continue; end
    %} 
    % mu first marginally faster - might be more so with more atoms (carbon, etc)
    % mu first stays faster with increasingly large iterations/accumulated atoms
    
    if err==0
        [ix,d] = knnsearch(kdt,rot2,'K',1,'SortIndices',0); % sort false might be faster
        %if any(d<15), er2=1; else er2=0; end % hard switch since no base value for er2
        if any(d<16), err=1; else err=0; end
    end
    
    if err==0; break; end
    end % % % end of block for test placement subfunct
    
    %{
    if err==0
    [err,muix] = mu_search(mu,rot2,tol,'short',0);
    err = any(err>0);
    end
    %}
    
    % if no collision, switch to place subunits as needed after replicating rotations
    if err==0
        [dyn] = dyncell(rot2,dyn); % write to partial list for fast lookups (obsolete under mu?)
        % no muix yet, not testing
        mu = mu_build(rot2,muix,mu,'leafmax',leaf,'maxdepth',2);
        
        % make write block into a generic subfunct and propogate to other atomfills?
        if subsel==0 % if complex, write each individual submodel after transforming
            tmpix = 1:numel(sel.id);
            % for assembly, instead construct tmpix with the randomized models needed
            for u=tmpix
                tmp = sel.adat{u};
                %if er2==1; tmp(:,4)=tmp(:,4)*2; end % diag collision test against membranes
                tmp(:,1:3) = tmp(:,1:3)*rotmat(init,spinang);
                tmp(:,1:3) = tmp(:,1:3)*rotmat(rotax,theta)+memloc;
                
                [split,dx] = dynsplit(tmp,split,dx,sel.id{u});
                list.(sel.id{u})(end+1,:) = memloc;
            end
        else % if not complex, use subsel index to write only that model
            tmp = sel.adat{subsel};
            %if er2==1; tmp(:,4)=tmp(:,4)*2; end % diag collision test against membranes
            tmp(:,1:3) = tmp(:,1:3)*rotmat(init,spinang);
            tmp(:,1:3) = tmp(:,1:3)*rotmat(rotax,theta)+memloc;
            
            [split,dx] = dynsplit(tmp,split,dx,sel.id{subsel});
            list.(sel.id{subsel})(end+1,:) = memloc;
        end
        
        count.s=count.s+1;
    else
        count.f=count.f+1;
    end
end

end
% cleanup and output stuff
fprintf('   membrane embedding complete, placed %i, failed %i \n',count.s,count.f);

% remove trailing zeros from atom registry and sparse tracker
f = fieldnames(split);
for i=1:numel(f)
    split.(f{i})(dx.(f{i}):end,:) = [];
end
dyn{1}(dyn{2}:end,:) = [];

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

function [dyn] = dyncell(addpts,dyn)
    l = size(addpts,1); e = dyn{2}+l-1;
    if e>size(dyn{1},1)
        dyn{1}(dyn{2}:(size(dyn{1},1)+l)*3,:) = 0;
    end
    dyn{1}(dyn{2}:e,:) = addpts;
    dyn{2} = dyn{2}+l;
end
function [split,dx] = dynsplit(tpts,split,dx,splitname) %slower than inlined a bit
tdx = dx.(splitname);
l = size(tpts,1); e = tdx+l-1;
if e>size(split.(splitname),1)
    split.(splitname)(tdx:(tdx+l)*2,:) = 0;
end
split.(splitname)(tdx:e,:) = tpts; dx.(splitname) = tdx+l;
end


function t = rotmat(ax,rad)
ax = ax/norm(ax);
x = ax(1); y = ax(2); z = ax(3);
c = cos(rad); s = sin(rad);

t1 = c + x^2*(1-c);
t2 = x*y*(1-c) - z*s;
t3 = x*z*(1-c) + y*s;
t4 = y*x*(1-c) + z*s;
t5 = c + y^2*(1-c);
t6 = y*z*(1-c)-x*s;
t7 = z*x*(1-c)-y*s;
t8 = z*y*(1-c)+x*s;
t9 = c+z^2*(1-c);

t = [t1 t2 t3
    t4 t5 t6
    t7 t8 t9];
end