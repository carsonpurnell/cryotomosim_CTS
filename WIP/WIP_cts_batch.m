% testing batch simulation runs
%% parameter setups
n = 3; % number to mod-sim
pix = 8.5; % currently fixed, should be easy to implement as variable
sz = [400,400,60]; % side lengths
% separate vector or second row to indicate variation in model size?
batchname = 'ATPs_mem';
targs = {'ATPS__flip.6j5i.membrane.cif'};

%targets = x; % probably most complicated, might need its own entry function
% for now, use only one fixed layer set to avoid even more complexity
% can a whole struct be shoved in as a parameter in place of several name-value pairs?

% how to parse so many inputs? how to randomize across the parameter space?
% if 1x1 nonrandom, if 1x2, flat x:y randomizer, if 2x1 use x+y*(rng), if >2 pick from the set?
pmod = param_model(pix,'layers',targs,'mem',[3,12]);
%batchmod = batchparam('layers',targs,'mem',[3,12],)
psim = param_simulate('pix',pix);
%%
batchsim = batchparam(n,'pix',pix,'dose',[60,150],'defocus',[-3,-5],'tilt',-60:3:60);

%% execute runs

for i=1:n
% model
suf = append(batchname,'_',string(i));
[cts,~,~,~,~,~,~,~,~,outfile] = cts_model_atomic(sz,pmod,'suffix',suf,'dynamotable',1);
[path,name,ext] = fileparts(outfile);
outfile = fullfile(path,append(name,'.atom.mat')); %bake into sim function?

% simulation
tmp = namedargs2cell(batchsim{i});
tsim = param_simulate(tmp{:});
cts_simulate_atomic(outfile,tsim,'suffix',append('sim_',string(i)));
end
% how to specify iterative simulations of the same model? cell array of input parameters?
% more likely extra option to input another set of simulation parameters


%% output/display?

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

%% internal functions
% the actual cts_batchrunner parts

function param = batchparam(n,varargin)
if rem(numel(varargin),2)==1, error('CTS batch params: bad number of args'); end

param = cell(n,1);
for i=1:n
    param{i} = struct(varargin{:}); % place the initial parameters
    f = fieldnames(param{i});
    for j=1:numel(f)
        if isequal(size(param{i}.(f{j})),[1,2])
            %disp('switch det')
            % hideous randomization, subfunct for it?
            param{i}.(f{j}) = param{i}.(f{j})(1)+rand*(param{i}.(f{j})(2)-param{i}.(f{j})(1));
            param{i}.(f{j}) = round(param{i}.(f{j}),2,'significant');
        else
            %param{i}.(f{j}) = param{i}.(f{j}); % do nothing same val
        end
    end
    % run param functions here to generate full params beforehand? need flag or separate mod/sim functs
end

param;
end

