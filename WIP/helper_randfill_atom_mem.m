function [split,dx,dyn] = helper_randfill_atom_mem(memdat,pmod,perim,sz)
%% prep stuff
split = struct; dx = struct; list = struct;
%pix = pmod.pix;

% combine filtering to only/no membranes into a subfunct and add to randfill
% maybe also combine with a split setup sub? or move split setup to here only?

for i=1:numel(pmod.layers) % filter to only membrane entries in each layer
    ix = zeros(1,numel(pmod.layers{i}));
    for j=1:numel(pmod.layers{i})
        ix(j) = any(contains([pmod.layers{i}(j).flags],'membrane'));
    end
    layers{i} = pmod.layers{i}(logical(ix));
end
ix = cellfun(@numel,layers);
layers = layers(logical(ix)); % prune out layers without membranes

%bail if no membranes
if numel(layers)<1, dyn = {perim,1};disp('no membrane structs to embed, skipping');return, end

% probably need cleanup dx/dyn/split stuff

% need to test multiple different layers
% each layer a different mixture - loop over each with the membrane loop

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
mu = mu_build(dyn{1},[0,0,0;sz*pmod.pix],'leafmax',leaf,'maxdepth',2);

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
    progstr = sprintf('progress %3.0f, membrane %i of %i',prog,j,numel(memdat.memcell));
    fprintf([progdel, progstr]);
    progdel = repmat(sprintf('\b'), 1, length(progstr));
    
memsel = j;%randi(numel(memdat)); % select membrane to place on

%memflags = memdat.table.class(j); if any(matches(memflags,'bare')); continue; end
if rand<memdat.table.bare(j); continue; end % skip if membrane doesn't hit dictionary fraction

%qq = vertcat(memdat.memcell{[1:j-1,1+j:end]});
kdt = KDTreeSearcher(vertcat(memdat.memcell{[1:j-1,1+j:end]}));

% placement attempt iterations, based on meshpts available and class fraction
iters = size(memdat.memcell{memsel},1)*.03*memdat.table.protfrac(j);

memcycle = randi(numel(layers));

for i=1:iters 
    sel = layers{memcycle}(randi(numel(layers{memcycle})));
    
    
    % inner loop: random axial rotation, rotation to transmembrane vector, collision test
    if any(ismember(sel.flags,'complex'))
        sel.sumperim = vertcat(sel.perim{:});
        subsel = 0;
    else
        subsel = randi(numel(sel.id));
        sel.sumperim = sel.perim{subsel};
    end
    
    % % % start of placement test block (might subfunct)
    for v=1:retry 
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
    
    [err,muix] = mu_search(mu,rot2,tol,'short',0);
    err = any(err>0);
    if err==1; continue; end
    if err==0
        [~,d] = knnsearch(kdt,rot2,'K',1,'SortIndices',0); % sort false might be faster
        if any(d<16), err=1; else err=0; end
    end
    if err==0; break; end
    end
    % % % end of block for test placement
    
    
    % if no collision, switch to place subunits as needed after replicating rotations
    if err==0
        [dyn] = dyncell(rot2,dyn); % write to partial list for fast lookups
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

% % % major bottleneck below here, dramatically slowing overall progress
f = fieldnames(split); % remove trailing zeros in atom registry
for i=1:numel(f)
    split.(f{i})(dx.(f{i}):end,:) = [];
end
[dyn] = dyncell(vertcat(memdat.memcell{:}),dyn); % add mems to dyn for collision detection
%dyn{1}(dyn{2}:end,:) = []; % remove trailing zeros in sparse tracker

%% cleanup membrane, remove lipid pseudoatoms replaced by membrane proteins
tmp = struct2cell(split); atoms = vertcat(tmp{:});

%kdt = KDTreeSearcher(memdat.atoms.vesicle(:,1:3)); %toc %11 12 15
kdt2 = KDTreeSearcher(atoms(:,1:3)); %toc %77 68 67

f = fieldnames(memdat.atoms);
for i=1:numel(f)
    %[ix,d] = knnsearch(kdt,atoms(:,1:3),'K',1,'SortIndices',0); %toc %554 119 121
    % bottleneck right here
    [~,d2] = knnsearch(kdt2,memdat.atoms.(f{i})(:,1:3),'K',1,'SortIndices',0); %toc
    
    ixf = d2>4.0; % appears to reliably prevent overlap without gaps
    mematoms = memdat.atoms.(f{i})(ixf,:);
    split.(f{i}) = mematoms;
end


end


%% internal functions
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