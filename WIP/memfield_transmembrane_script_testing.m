%% placing memprots with new atomic membrane structure
pix = 12;
targ = {'ATPS.membrane.complex.mat'};
%targ = {'GABAar.membrane.complex.mat'};
pmod = param_model(pix,'layers',targ);

sz = [200,200,100];
memdat = gen_mem_atom(sz,pix,'num',3);
% needs a bit more work, a few vectors (probably due to corners) are not well-oriented - denser mesh?
% alternate method - dense surface mesh of expanded membrane hull, remove inner points, get nearest?
% would need to be very dense. but could average with the near-3 result to cover most cases?
rng(5)

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

init = [0,0,1];
sel = pmod.layers{1}(1);

for j=1:numel(memdat.memcell)
memsel = j;%randi(numel(memdat)); % select membrane to place on

qq = vertcat(memdat.memcell{[1:j-1,1+j:end]});
kdt = KDTreeSearcher(qq);

if any(ismember(sel.flags,'complex'))
    sel.sumperim = vertcat(sel.perim{:});
    subsel = 0;
else
    sel.sumperim = sel.perim{1};
    subsel = randi(numel(sel.id));
end
% start off preselecting coords from it or just start running through them? 
% they are not spatially ordered
iters = 150;
count.s = 0; count.f = 0;
for i=1:iters
    % inner loop: random axial rotation, rotation to transmembrane vector, collision test
    % mem vector and figuring stuff
    memloc = memdat.memcell{memsel}(i,:); % selected coordinate
    surfvec = memdat.normcell{memsel}(i,:); % normal vector to surface at coordinate
    
    % need to funct out placement testing, and allow multiple attempts per iteration
    
    rotax=cross(init,surfvec); %compute the normal axis from the rotation angle
    theta = -acos( dot(init,surfvec) ); % angle between ori vec and surface
    
    spinang = rand*2*pi;
    rot1 = sel.sumperim*rotmat(init,spinang); % random axial rotation
    rot2 = rot1*rotmat(rotax,theta)+memloc;
    
    %diagori = init*rotmat(rotax,theta);
    
    %
    tol = 2;
    [err,muix] = mu_search(mu,rot2,tol,'short',0);
    err = any(err>0);
    %} 
    %mu first faster overall, mu is faster per iteration so saves time avoiding knn
    %
    
    %if err==0
        [ix,d] = knnsearch(kdt,rot2,'K',1,'SortIndices',0); % sort false might be faster
        %if any(d<15), er2=1; else er2=0; end % hard switch since no base value for er2
        if any(d<15), err=1; else err=0; end
    %end
    
    %{
    if err==0
    tol = 2;
    [err,muix] = mu_search(mu,rot2,tol,'short',0);
    err = any(err>0);
    end
    %}
    %20/6 21/6 mu first
    %11/19 11/17 10/17 mu second
    
    % if no collision, switch to place subunits as needed after replicating rotations
    if err==0
        [dyn] = dyncell(rot2,dyn); % write to partial list for fast lookups (obsolete under mu?)
        % no muix yet, not testing
        mu = mu_build(rot2,muix,mu,'leafmax',leaf,'maxdepth',2);
        
        if subsel==0
            for u=1:numel(sel.id)
                tmp = sel.adat{u};
                %if er2==1; tmp(:,4)=tmp(:,4)*2; end % diag collision test against membranes
                tmp(:,1:3) = tmp(:,1:3)*rotmat(init,spinang);
                tmp(:,1:3) = tmp(:,1:3)*rotmat(rotax,theta)+memloc;
                
                [split,dx] = dynsplit(tmp,split,dx,sel.id{u});
                list.(sel.id{u})(end+1,:) = memloc;
            end
        else
            % use subsel index to place that item alone
        end
        %tpts = sel.adat{sub};
        %tpts(:,1:3) = transformPointsForward(tform,tpts(:,1:3))+loc;
        
        count.s=count.s+1;
    else
        count.f=count.f+1;
        %disp('fail, something borked')
    end
end

end

% remove trailing zeros from atom registry and sparse tracker
f = fieldnames(split);
for i=1:numel(f)
    split.(f{i})(dx.(f{i}):end,:) = [];
end
dyn{1}(dyn{2}:end,:) = [];

disp(count)
[vol,solv] = helper_atoms2vol(pix,split,sz*pix);
mvol = helper_atoms2vol(pix,memdat.atoms,sz*pix);
sliceViewer(max(vol,mvol));

diagpts = [init;rotax;surfvec];
% plot3p(list.ATPS_head,'o'); hold on; plot3p(list.ATPS_head+memdat.normcell{1}(1:50,:)*100,'.'); % diag vecs
% plot3p(dyn{1}(1:dyn{2}-1,:),'.'); hold on; plot3p(memdat.memcell{1},'.'); % diag placements
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