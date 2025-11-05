%% placing memprots with new atomic membrane structure
pix = 10;
targ = {'ATPS.membrane.complex.mat'};
pmod = param_model(pix,'layers',targ);

sz = [200,200,50];
memdat = gen_mem_atom(sz,pix);
% needs a bit more work, a few vectors (probably due to corners) are not well-oriented

%%
for i=1:numel(pmod.layers)
    for j=1:numel(pmod.layers{i}.id)
        atoms.(pmod.layers{i}.id{j}) = zeros(0,4);
    end
end
dyn{1} = single(zeros(0,3)); dyn{2} = 0;
leaf = 1e3;
mu = mu_build(dyn{1},[0,0,0;sz*pix],'leafmax',leaf,'maxdepth',2);

%pick a membrane
memsel = randi(numel(memdat));

sel = pmod.layers{1}(1);
if any(ismember(sel.flags,'complex'))
    sel.sumperim = vertcat(sel.perim{:});
    complex = 1;
else
    complex = 0;
end
% start off preselecting coords from it or just start running through them? they are not spatially ordered
iters = 15;
for i=1:iters
    % inner loop: random axial rotation, rotation to transmembrane vector, collision test
    r1 = rand*2*pi;
    rot1 = sel.sumperim*rotmat([0,0,1],r1); % random axial rotation
    
    
    
    % if no collision, switch to place subunits as needed after replicating rotations
    if err==0
        tpts = sel.adat{sub};
        tpts(:,1:3) = transformPointsForward(tform,tpts(:,1:3))+loc;
        
        [dyn] = dyncell(ovcheck,dyn);
        mu = mu_build(ovcheck,muix,mu,'leafmax',leaf,'maxdepth',2);
        [split,dx] = dynsplit(tpts,split,dx,sel.modelname{sub});
        list.(sel.modelname{sub})(end+1,:) = loc;
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
    end
end



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