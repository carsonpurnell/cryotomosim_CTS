function [vol,solv,atlas,splitvol,acount,split] = cts_model_atomic(box,input,param)

arguments
    box
    input = 'gui'
    %pix = 10
    param = {10}
end

if strcmp(input,'gui')
    filter = '*.pdb;*.pdb1;*.mrc;*.cif;*.mmcif;*.mat';
    input = util_loadfiles(filter);
end

if iscell(param), param = param_model(param{:},'layers',input); end
pix = param.pix;
for i=1:numel(input)
    particles(i) = helper_pdb2dat(input{i},pix,2,0,0);
end


layers{1} = particles; fprintf('loaded %i structure files  ',numel(input));

%% functionalized model gen part
%rng(7); 
con = 1;
boxsize = pix*box*1;
[splitin.carbon,dyn] = gen_carbon(boxsize); % atomic carbon grid generator
memnum = 5;
memnum = param.mem;
tic; [splitin.lipid,kdcell,shapecell,dx.lipid,dyn] = modelmem(memnum,dyn,boxsize); toc;

if con==1
    con = helper_atomcon(boxsize,pix,0,0); % pseudonatural ice border (wavy flat, no curvature)
    dyn{1}(dyn{2}:dyn{2}+size(con,1)-1,:) = con; dyn{2}=dyn{2}+size(con,1)-1;
end

n = 100;
n = param.iters(1);
%profile on
%tic; [split,dyn,mu] = fn_modelgen(layers,boxsize,n,splitin,dx,dyn); toc
tic; [split,dyn,mu] = helper_randfill_atom(layers,boxsize,n,splitin,dx,dyn); toc
%profile viewer

%% function for vol, atlas, and split generation + water solvation
outpix = pix;
%profile on
[vol,solv,atlas,splitvol,acount] = helper_atoms2vol(outpix,split,boxsize);
%profile viewer
%sliceViewer(vol*1+solv);

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

function [dyn] = dyncell(addpts,dyn)
    l = size(addpts,1); e = dyn{2}+l-1;
    if e>size(dyn{1},1)
        dyn{1}(dyn{2}:(size(dyn{1},1)+l)*3,:) = 0;
    end
    dyn{1}(dyn{2}:e,:) = addpts;
    dyn{2} = dyn{2}+l;
end

function [pts,kdcell,shapecell,dx,dyn] = modelmem(memnum,dyn,boxsize)
dyn = {dyn,size(dyn,1)}; %convert to dyncell
pad = [0,0,0;boxsize];
dyn{1} = [dyn{1};pad];
kdcell = []; shapecell = [];
mu = mu_build(dyn{1},'leafmax',1e3,'maxdepth',2);

%memparam = {zeros(3,2)};
memparam = {[300,500;.4,.6;24,6],1};

tol = 2; %tolerance for overlap testing
retry = 5; %retry attempts per iteration
count.s = 0; count.f = 0;
lipid{1} = zeros(0,4); lipid{2} = 1;
for i=1:memnum % simplified loop to add vesicles
    %[tpts,tperim] = gen_mem(300+randi(500),[],rand*0.6+0.4, 24+randi(6));%.3/.7 and 24/8
    [tpts,tperim] = gen_mem(memparam{1,1});
    
    [err,loc,tform,ovcheck,muix] = anyloc(boxsize,tperim,dyn,retry,tol,mu);
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
        mu = mu_build(ovcheck,muix,mu,'leafmax',1e3,'maxdepth',2);
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