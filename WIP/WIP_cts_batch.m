% testing batch simulation runs
%% parameter setups
n = 3; % number to mod-sim
%pix = 8.5; % currently fixed, should be easy to implement as variable
sz = [400,400,50]; % side lengths
% separate vector or second row to indicate variation in model size?
batchname = 'ATPs_mem';
targs = {'ATPS__flip.6j5i.membrane.cif'};%'tubulin__1tub.distract.mat','act1-A2.distract.mat'...
    %'1trv_thioredoxin.distract.pdb','ribo__ribo__4ug0_4v6x.group.mat','7b5s.distract.mat'};

%targets = x; % probably most complicated, might need its own entry function
% for now, use only one fixed layer set to avoid even more complexity
% can a whole struct be shoved in as a parameter in place of several name-value pairs?

% how to parse so many inputs? how to randomize across the parameter space?
% if 1x1 nonrandom, if 1x2, flat x:y randomizer, if 2x1 use x+y*(rng), if >2 pick from the set?
%pmod = param_model(8,'layers',targs,'mem',[3,12]);
%batchmod = batchparam('layers',targs,'mem',[3,12],)
%psim = param_simulate('pix',pix);
%%
% probably should split batch, special handle targs, maybe pixel size?
%batchmod = batchparam(n,1,'layers',{targs},'pix',[8,9],'iters',[200,1000],'mem',[2,12]);
%batchsim = batchparam(n,2,'dose',[60,150],'defocus',[-3,-5],'scatter',[0.5,1.5],'tilt',-60:3:60);

%batchmod = batchparam_mod(n,'layers',{targs},'pix',[8,9],'iters',[200,1000],'mem',[2,12]);
%batchsim = batchparam_sim(n,'dose',[60,150],'defocus',[-3,-5],'scatter',[0.5,1.5],'tilt',-60:3:60);

batchmod = batchparam(n,'layers',{targs},'pix',[8,9],'iters',[200,1000],'mem',[2,12]);
batchsim = batchparam(n,'dose',[60,150],'defocus',[-3,-5],'scatter',[0.5,1.5],'tilt',-60:3:60);

opt.ideal = param_simulate('dose',500,'defocus',-4,'raddamage',0,'scatter',0.5,'tilt',-80:2:80);
opt.ideal = 0;

%% execute runs

for i=1:n
    %pmod.pix = batchmod{i}.pix; pmod.iters(2) = batchmod{i}.iters;
    suf = append(batchname,'_',string(i));
% model

[cts,~,~,~,~,~,~,~,~,outfile] = cts_model_atomic(sz,batchmod{i},'suffix',suf,'dynamotable',1);
[path,name,ext] = fileparts(outfile);
outfile = fullfile(path,append(name,'.atom.mat')); %bake into sim function?

% simulation
batchsim{i}.pix = cts.param.pix;
%tmp = namedargs2cell(batchsim{i});
%tsim = param_simulate(tmp{:});
cts_simulate_atomic(outfile,batchsim{i},'suffix',append('sim_',string(i)));
% ideal sim run
    if isstruct(opt.ideal) % run ideal sim if argument given
        %should already be a consolidated param? or allow a cell array of values?
        isim = opt.ideal; isim.pix = cts.param.pix;
        cts_simulate_atomic(outfile,isim,'suffix',append('ideal_',string(i)));
    end
    fprintf('CTS batch: finished %i of %i runs\n',i,n)
end
fprintf('done batch of %i runs\n',n)

%% output/display?



%% internal functions
% the actual cts_batchrunner parts

function param = batchparam_mod(n,varargin)
if rem(numel(varargin),2)==1, error('CTS batch params: bad number of args'); end
param = cell(n,1);
for i=1:n
    tmp = batchrand(varargin);
    %{
    param{i} = struct(varargin{:}); % place the initial parameters
    f = fieldnames(param{i});
    for j=1:numel(f)
        if isequal(size(param{i}.(f{j})),[1,2]) && ~iscell(param{i}.(f{j}))
            % hideous randomization, subfunct for it?
            param{i}.(f{j}) = param{i}.(f{j})(1)+rand*(param{i}.(f{j})(2)-param{i}.(f{j})(1));
            param{i}.(f{j}) = round(param{i}.(f{j}),2,'significant');
        else
            %param{i}.(f{j}) = param{i}.(f{j}); % do nothing same val
        end
    end
    tmp = namedargs2cell(param{i});
    %}
    ix = find(strcmp(tmp,'pix'));
    pix = tmp{ix+1};
    param{i} = param_model(pix,tmp{:});
end
end

function param = batchparam_sim(n,varargin)
if rem(numel(varargin),2)==1, error('CTS batch params: bad number of args'); end
param = cell(n,1);
for i=1:n
    tmp = batchrand(varargin);
    %{
    param{i} = struct(varargin{:}); % place the initial parameters
    % the block below could be a generic local funct for each batcher to do randomization
    f = fieldnames(param{i});
    for j=1:numel(f)
        if isequal(size(param{i}.(f{j})),[1,2]) && ~iscell(param{i}.(f{j}))
            % hideous randomization, subfunct for it?
            param{i}.(f{j}) = param{i}.(f{j})(1)+rand*(param{i}.(f{j})(2)-param{i}.(f{j})(1));
            param{i}.(f{j}) = round(param{i}.(f{j}),2,'significant');
        else
            %param{i}.(f{j}) = param{i}.(f{j}); % do nothing same val
        end
    end
    tmp = namedargs2cell(param{i});
    %}
    param{i} = param_simulate(tmp{:});
end
end

function param = batchparam(n,varargin)
if rem(numel(varargin),2)==1, error('CTS batch params: bad number of args'); end
param = cell(n,1);
ix = find(strcmp(varargin,'pix'));
for i=1:n
    tmp = batchrand(varargin);
    if ix>0
        pix = tmp{ix+1};
        param{i} = param_model(pix,tmp{:});
    else
        param{i} = param_simulate(tmp{:});
    end
end
end

function paramcell = batchrand(vars)
paramcell = struct(vars{:}); % place the initial parameters
f = fieldnames(paramcell);
for j=1:numel(f)
    if isequal(size(paramcell.(f{j})),[1,2]) && ~iscell(paramcell.(f{j}))
        % hideous randomization, subfunct for it?
        paramcell.(f{j}) = paramcell.(f{j})(1)+rand*(paramcell.(f{j})(2)-paramcell.(f{j})(1));
        paramcell.(f{j}) = round(paramcell.(f{j}),2,'significant');
    else
        %param{i}.(f{j}) = param{i}.(f{j}); % do nothing same val
    end
end
paramcell = namedargs2cell(paramcell);
end

function param = batchparam_both(n,ptype,varargin)
if rem(numel(varargin),2)==1, error('CTS batch params: bad number of args'); end
%{
if ptype==1
    % match which cells have pix/layers string and add 1 for the value
    ix = find(strcmp(varargin,'pix'));
    pix = varargin{ix+1};
    ix = find(strcmp(varargin,'layers'));
    layers = varargin{ix+1}{:}
    modpreload = param_model(pix,'layers',layers);
    varargin{ix+1} = modpreload.layers;
end
%}

param = cell(n,1);
for i=1:n
    param{i} = struct(varargin{:}); % place the initial parameters
    f = fieldnames(param{i});
    for j=1:numel(f)
        if isequal(size(param{i}.(f{j})),[1,2]) && ~iscell(param{i}.(f{j}))
            %disp('switch det')
            % hideous randomization, subfunct for it?
            param{i}.(f{j}) = param{i}.(f{j})(1)+rand*(param{i}.(f{j})(2)-param{i}.(f{j})(1));
            param{i}.(f{j}) = round(param{i}.(f{j}),2,'significant');
        else
            %param{i}.(f{j}) = param{i}.(f{j}); % do nothing same val
        end
    end
    % run param functions here to generate full params beforehand? need flag or separate mod/sim functs
    tmp = namedargs2cell(param{i});
    
    if ptype==1
        param{i} = param_model(param{i}.pix,tmp{:});
    else
        param{i} = param_simulate(tmp{:});
    end
end

end
% test for structs, it does work even with missing arguments if converted to cell first
%{
clear test
test.a = 5; test.b = 'q';
test = namedargs2cell(test);
tfunct(test{:})
function tfunct(opt)
arguments
    opt.a = 0
    opt.b = 0
    opt.c = 0
end
opt
end
%}
