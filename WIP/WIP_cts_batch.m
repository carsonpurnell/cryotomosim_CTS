% testing batch simulation runs
%% parameter setups
n = 3; % number to mod-sim
pix = 8.5; % currently fixed, should be easy to implement as variable
sz = [500,500,70]; % side lengths
% separate vector or second row to indicate variation in model size?
batchname = 'ATPs_mem';
targs = {'ATPS__flip.6j5i.membrane.cif'};

%targets = x; % probably most complicated, might need its own entry function
% for now, use only one fixed layer set to avoid even more complexity
% can a whole struct be shoved in as a parameter in place of several name-value pairs?

% how to parse so many inputs? how to randomize across the parameter space?
% if 1x1 nonrandom, if 1x2, flat x:y randomizer, if 2x1 use x+y*(rng), if >2 pick from the set?
pmod = param_model(pix,'layers',targs,'mem',[3,12]);
psim = param_simulate('pix',pix);

%% execute runs

for i=1:n
% model
suf = append(batchname,'_',string(i));
[cts,~,~,~,~,~,~,~,~,outfile] = cts_model_atomic(sz,pmod,'suffix',suf,'dynamotable',1);
[path,name,ext] = fileparts(outfile);
outfile = fullfile(path,append(name,'.atom.mat')); %bake into sim function?

% simulation
cts_simulate_atomic(outfile,psim,'suffix',append('sim_',string(i)));
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



