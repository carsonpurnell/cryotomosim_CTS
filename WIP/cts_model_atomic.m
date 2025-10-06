function [cts,vol,solv,atlas,splitvol,acount,split,dat,list,outfile] = cts_model_atomic(box,input,param,opt)
% [vol,solv,atlas,splitvol,acount,split,dat,list] = cts_model_atomic(box,input,param,opt)
% ex
% [vol,solv,atlas] = cts_model_atomic([300,400,50],'gui',{8,'mem',1,'iters',400});
%

arguments
    box
    input = 'gui' %DEPRECATED
    %pix = 10
    param = {10}
    opt.suffix = ''
    opt.con = 1
    opt.ermem = 0
    opt.outdir = []
end
%{
if strcmp(input,'gui')
    filter = '*.pdb;*.pdb1;*.mrc;*.cif;*.mmcif;*.mat';
    input = util_loadfiles(filter);
end
%}

if strcmp(opt.outdir,'gui') % custom output locations
    outdir = uigetdir(); cd(outdir)
elseif isempty(opt.outdir)
    outdir = fullfile(getenv('HOME'),'tomosim');
    cd(getenv('HOME'));
    if ~isfolder('tomosim'), mkdir tomosim; end, cd tomosim
else
    outdir = opt.outdir; cd(outdir)
end

if iscell(param), param = param_model(param{:}); end
pix = param.pix;
%{
%if numel(param.layers)>0
    for i=1:numel(param.layers)
        filter = '*.pdb;*.pdb1;*.mrc;*.cif;*.mmcif;*.mat';
        input = util_loadfiles(filter);
        particles(i) = helper_pdb2dat(input,pix,2,0,0);
    end
%end
layers{1} = particles; 
%}
fprintf('loaded %i layers of particles  \n',numel(param.layers));


% functionalized model gen part 
%con = 0;
boxsize = pix*box*1;
% need a working toggle to setup split/dyn without carbon or mem
% do constraint first to prep dyn? would need to rejigger modelmem placement testing
if param.grid~=0
    [splitin.carbon,dyn] = gen_carbon(boxsize); % atomic carbon grid generator
else
    dyn = zeros(1,3);
end
%memnum = 5;
memnum = param.mem;
tic; 

%[splitin,kdcell,shapecell,dx,dyn] = modelmem(memnum,dyn,boxsize,opt.ermem); 
dyn = {dyn,size(dyn,1)}; %convert to dyncell
pad = [0,0,0;boxsize];
dyn{1} = [dyn{1};pad];

% new mem stuff testing, need better parsing and passing args from model_param etc
[memdat] = gen_mem_atom(box,pix,'num',memnum);

splitin = memdat.atoms;
f = fieldnames(memdat.atoms);
for i=1:numel(f)
    tmp = memdat.atoms.(f{i}); n = size(tmp,1);
    ix = randperm(n); ix = ix(1:round(n/2)); % 1/6 of atoms for collision detection later
    dx.(f{i}) = size(tmp,1)+1;
    dyn{1} = [dyn{1};tmp(ix,1:3)]; dyn{2} = numel(ix)+1;
end

toc;
%size(dyn)
%dyn
%dx
%splitin has entries, but are all 0?!
%[mi1,ma1] = bounds(splitin.lipid,1)
%[mi2,ma2] = bounds(splitin.ERmem,1)
%return 
if opt.con==1
    con = helper_atomcon(boxsize,pix,0,0); % pseudonatural ice border (wavy flat, no curvature)
    dyn{1}(dyn{2}:dyn{2}+size(con,1)-1,:) = con; dyn{2}=dyn{2}+size(con,1)-1;
end

%n = 100;
n = param.iters;
%profile on

%tic; [split,dyn,mu] = fn_modelgen(layers,boxsize,n,splitin,dx,dyn); toc
%for i=1:numel(param.layers)
tic; [split,dyn,mu,list] = helper_randfill_atom(param.layers,boxsize,n,splitin,dx,dyn); toc
% list will be broken
%end
%profile viewer


%% function for vol, atlas, and split generation + water solvation
% preprune split to eliminate empty bins (membrane/carbon)
f = fieldnames(split);
for i=1:numel(f)
    if size(split.(f{i}),1)==0
        split = rmfield(split,f{i});
    end
    %[mi,ma] = bounds(split.(f{i}),1)
end
%split has stuff in it, despite not contributing to vol
[vol,solv,atlas,splitvol,acount] = helper_atoms2vol(pix,split,boxsize);
%sliceViewer(vol*1+solv);

%% folder and file generation stuff
time = string(datetime('now','Format','yyyy-MM-dd''t''HH.mm')); %timestamp
ident = char(strjoin(fieldnames(split),'_')); %combine target names to one string
if length(ident)>60, ident=ident(1:60); end %truncation check to prevent invalidly long filenames
if ~strncmp('_',opt.suffix,1), opt.suffix = append('_',opt.suffix); end
foldername = append('model_',time,'_',ident,'_pixelsize_',string(pix),opt.suffix); %combine info for folder name

%move to output directory in user/tomosim
%cd(getenv('HOME')); if ~isfolder('tomosim'), mkdir tomosim; end, cd tomosim
mkdir(foldername); cd(foldername);

WriteMRC(vol,pix,append(ident,opt.suffix,'.mrc'))
dat.box = boxsize;
dat.data = split;
save(append(ident,opt.suffix,'.atom.mat'),'dat','-v7.3')
cts.vol = vol+solv; cts.splitmodel = splitvol; cts.param.pix = pix;
cts.model.particles = vol; cts.model.ice = solv;
cts.list = list;
cts.param = param;

outname = append(ident,opt.suffix,'.mat');
outfile = fullfile(getenv('HOME'),'tomosim',foldername,outname);
save(outname,'cts','-v7.3')


f = fieldnames(cts.list);
for i=1:numel(f)
    coordmatrix = cts.list.(f{i});
    if ~isempty(coordmatrix)
        writematrix(coordmatrix,append('zcoords_',f{i},'.csv'));
    end
end
%output text file of input informations?
cd(userpath)

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

function [pts,kdcell,shapecell,dx,dyn] = modelmem(memnum,dyn,boxsize,ermem)
dyn = {dyn,size(dyn,1)}; %convert to dyncell
pad = [0,0,0;boxsize];
dyn{1} = [dyn{1};pad];
kdcell = []; shapecell = [];
mu = mu_build(dyn{1},'leafmax',1e3,'maxdepth',2);

%memparam = {zeros(3,2)};
memparam = {[300,500;.4,.6;28,6],'lipid'}; %size, sp, thick?)
lipid.(memparam{1,2}){1} = zeros(0,4); lipid.(memparam{1,2}){2} = 1;
if ermem==1
    memparam(2,:) = {[200,500;.1,.2;12,2],'ERmem'};
    lipid.(memparam{2,2}){1} = zeros(0,4); lipid.(memparam{2,2}){2} = 1;
end

tol = 2; %tolerance for overlap testing
retry = 5; %retry attempts per iteration
count.s = 0; count.f = 0;
%lipid{1} = zeros(0,4); lipid{2} = 1;
for i=1:memnum % simplified loop to add vesicles
    mtype = randi(size(memparam,1));
    %[tpts,tperim] = gen_mem(300+randi(500),[],rand*0.6+0.4, 24+randi(6));%.3/.7 and 24/8
    [tpts,tperim] = gen_mem(memparam{mtype,1});
    
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
        %[miQ,maQ] = bounds(tpts,1) NOT full of zeros
        
        [dyn] = dyncell(ovcheck,dyn);
        [lipid.(memparam{mtype,2})] = dyncell(tpts,lipid.(memparam{mtype,2}));
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
%lipid
%[mi1,ma1] = bounds(lipid.lipid{1},1)
%[mi2,ma2] = bounds(lipid.ERmem{1},1)
for i=1:size(memparam,1)
    tmp = lipid.(memparam{i,2}){2};
    lipid.(memparam{i,2}) = lipid.(memparam{i,2}){1}(1:tmp-1,:);
    %lipid.(memparam{i,2}){1}(tmp:end,:) = [];
    %lipid.(memparam{i,2})(2) = [];
    dx.(memparam{i,2}) = tmp;
    %lipid.(memparam{i,2}){1}(1:lipid{2}-1,:);
end
%lipid
%dx
pts = lipid; % all vals reduced to 0 at this point
%[mi3,ma3] = bounds(pts.lipid,1)
%[mi4,ma4] = bounds(pts.ERmem,1)
%pts = lipid{1}(1:lipid{2}-1,:); dx = lipid{2};
%dx.lipid = lipid.lipid{2};
fprintf('vesicles: placed %i, failed %i  ',count.s,count.f)
end