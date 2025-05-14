function [split,dyn,mu,list] = helper_randfill_atom(layers,boxsize,niter,split,dx,dyn)
if nargin<6
    dyn{1} = single(zeros(0,3)); dyn{2} = 0;
end
leaf = 1e3;
mu = mu_build(dyn{1},[0,0,0;boxsize],'leafmax',leaf,'maxdepth',2); %speed drops DRAMATICALLY at depth 3

if nargin<4
    split = struct; %ixincat = 1; %dynpts = single(zeros(0,3));
else
    fn = fieldnames(split);
    for i=1:numel(fn) %add split into dynpts
        %s = size(split.(fn{i})(:,1:3),1);
        %ix = randi(s,round(s/50),1); ix = unique(ix); %tmp = split.(fn{i})(ix,1:3);
        %l = size(tmp,1); %dynpts(end+1:end+l,:) = tmp; %dynpts = [dynpts;tmp];
    end
end
%ixincat = size(dynpts,1)+1; %deprecated counter for starting point for dyncat insertions

tol = 2; %tolerance for overlap testing
retry = 4; %retry attempts per iteration
count.s = 0; count.f = 0; %success and failure counters
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

for lc = 1:numel(layers)
    particles = layers{lc};
    n = niter(lc);
    
for i=1:n
    if rem(i,n/20)==0; fprintf('%i,',i); end
    which=randi(numel(particles));
    sel = particles(which); 
    sub = randi(numel(sel.adat));
    [err,loc,tform,ovcheck,muix] = anyloc(boxsize,sel.perim{sub},dyn,retry,tol,mu); % just as fast
    
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
        count.s=count.s+1;
    else
        count.f=count.f+1;
    end
end
%{
if lc==1
    %tic; sh=alphaShape(double(dynpts),12); toc; %plot(sh); drawnow;
end
%}
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
