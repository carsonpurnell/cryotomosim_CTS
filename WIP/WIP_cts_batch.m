% testing batch simulation runs
%% parameter setups
n = 3; % number to mod-sim
pix = 10; % currently fixed, should be easy to implement as variable
sz = [400,400,60]; % side lengths
% separate vector or second row to indicate variation in model size?

%targets = x; % probably most complicated, might need its own entry function
% for now, use only one fixed layer set to avoid even more complexity
% can a whole struct be shoved in as a parameter in place of several name-value pairs?

% how to parse so many inputs? how to randomize across the parameter space?
% if 1x1 nonrandom, if 1x2, flat x:y randomizer, if 2x1 use x+y*(rng), if >2 pick from the set?
pmod = param_model(8);
psim = param_simulate;

%% execute runs

% model
% simulation
% how to specify iterative simulations of the same model? cell array of input parameters?
% more likely extra option to input another set of simulation parameters


%% output/display?

% test for structs, it does work even with missing arguments if converted to cell first
clear test
test.a = 5; test.b = 'q';
test = namedargs2cell(test);
tfunct(test{:})

%% internal functions
% the actual cts_batchrunner parts



function tfunct(opt)
arguments
    opt.a = 0
    opt.b = 0
    opt.c = 0
end
opt
end