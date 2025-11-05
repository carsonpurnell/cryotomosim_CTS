%% placing memprots with new atomic membrane structure
pix = 12;
targ = {'ATPS.membrane.complex.mat'};
targ = {'GABAar.membrane.complex.mat'};
pmod = param_model(pix,'layers',targ);

sz = [300,300,60];
memdat = gen_mem_atom(sz,pix);
% needs a bit more work, a few vectors (probably due to corners) are not well-oriented - denser mesh?

%%
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

for i=1:numel(pmod.layers)
    for j=1:numel(pmod.layers{i}.id)
        atoms.(pmod.layers{i}.id{j}) = zeros(0,4);
    end
end
dyn{1} = single(zeros(0,3)); dyn{2} = 1;
leaf = 1e3;
mu = mu_build(dyn{1},[0,0,0;sz*pix],'leafmax',leaf,'maxdepth',2);

%pick a membrane
memsel = randi(numel(memdat));
init = [0,0,-1];

sel = pmod.layers{1}(1);
if any(ismember(sel.flags,'complex'))
    sel.sumperim = vertcat(sel.perim{:});
    complex = 1;
else
    sel.sumperim = sel.perim{1};
    complex = 0;
end
% start off preselecting coords from it or just start running through them? they are not spatially ordered
iters = 50;
for i=1:iters
    % inner loop: random axial rotation, rotation to transmembrane vector, collision test
    % mem vector and figuring stuff
    memloc = memdat.memcell{memsel}(i,:); % selected coordinate
    surfvec = memdat.normcell{memsel}(i,:); % normal vector to surface at coordinate
    
    rotax=cross(init,surfvec); %compute the normal axis from the rotation angle
    theta = acos( dot(init,surfvec) ); % angle between ori vec and surfae
    
    spinang = rand*2*pi;
    rot1 = sel.sumperim*rotmat(init,spinang); % random axial rotation
    rot2 = rot1*rotmat(rotax,theta)+memloc;
    
    
    
    err = 0; % force placement, no collision test yet
    
    % if no collision, switch to place subunits as needed after replicating rotations
    if err==0
        [dyn] = dyncell(rot2,dyn); % write to partial list for fast lookups (obsolete under mu?)
        % no muix yet, not testing
        % mu = mu_build(rot2,muix,mu,'leafmax',leaf,'maxdepth',2);
        
        if complex ==1
            for u=1:numel(sel.id)
                tmp = sel.adat{u};
                tmp(:,1:3) = tmp(:,1:3)*rotmat(init,spinang);
                tmp(:,1:3) = tmp(:,1:3)*rotmat(rotax,theta)+memloc;
                
                [split,dx] = dynsplit(tmp,split,dx,sel.id{u});
                list.(sel.id{u})(end+1,:) = memloc;
            end
        else
            
        end
        %tpts = sel.adat{sub};
        %tpts(:,1:3) = transformPointsForward(tform,tpts(:,1:3))+loc;
        
        
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
        %count.s=count.s+1;
    else
        %count.f=count.f+1;
        disp('fail, something borked')
    end
end
[vol,solv] = helper_atoms2vol(pix,split,sz*pix);
mvol = helper_atoms2vol(pix,memdat.atoms,sz*pix);
sliceViewer(max(vol,mvol));


%% internal functions


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