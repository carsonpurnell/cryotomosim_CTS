% script for generating some of the figures in the CTS segmentation paper

%% model parameters
sz = [400,400,50]; % model size
pix = 10; % pixel size
%{
targets = {'ribo__ribo__4ug0_4v6x.group.mat',... % two ribo variations
    'cofilactin__cofilactin__actin__actin__x3-x4_long.bundle.pdb'... % two lengths each
    'MT__6o2tx3.mat',... % long MT stretch
    'tric__tric__6nra-open_7lum-closed.group.mat'};
%}
targets = {'ribo__ribo__4ug0_4v6x.group.mat'}; % currently only using ribos for clean&fast testing
%targets  1; % for using GUI to load a single layer instead

%% SNR/dose
n = 2; % number of training samples, one extra will be created as independent test data
varname = 'dose'; % name of the simulation variable to modify
variters = [10,25,75,100,200]; % values of the variable to iterate over for each sample

pmod = param_model(pix,'layers',targets,'iters',500,'mem',0,'grid',0);
psim = param_simulate(); % currently all default values

digits = numel(num2str(n)); fspec = append('%0',num2str(digits),'i'); % naming stuff
for i=1:n+1
    suf = append('var_',varname,'_',num2str(i,fspec));
    [cts,outfile] = cts_model(zeros(sz),pmod,'suffix',suf);
    
    for j=1:numel(variters)
        psim.(varname) = variters(j); % get iteration value for variable
        sufsim = append('var_',varname,'_',string(variters(j)));
        [~,~,~,atlas] = cts_simulate(outfile,psim,'suffix',sufsim);
    end
end


%% defocus
n = 2; % number of training samples, one extra will be created as independent test data
varname = 'defocus'; % name of the simulation variable to modify
variters = [-1,-3,-5,-7,-9]; % values of the variable to iterate over for each sample

pmod = param_model(pix,'layers',targets,'iters',500,'mem',0,'grid',0);
psim = param_simulate(); % currently all default values

digits = numel(num2str(n)); fspec = append('%0',num2str(digits),'i'); % naming stuff
for i=1:n+1
    suf = append('var_',varname,'_',num2str(i,fspec));
    [cts,outfile] = cts_model(zeros(sz),pmod,'suffix',suf);
    
    for j=1:numel(variters)
        psim.(varname) = variters(j); % get iteration value for variable
        sufsim = append('var_',varname,'_',string(variters(j)));
        [~,~,~,atlas] = cts_simulate(outfile,psim,'suffix',sufsim);
    end
end

%% pixel size
n = 2; % number of training samples, one extra will be created as independent test data
varname = 'pix'; % name of the simulation variable to modify
variters = [8,10,12,14,16]; % values of the variable to iterate over for each sample

pmod = param_model(pix,'layers',targets,'iters',500,'mem',0,'grid',0);
psim = param_simulate(); % currently all default values

digits = numel(num2str(n)); fspec = append('%0',num2str(digits),'i'); % naming stuff
for i=1:n+1
    suf = append('var_',varname,'_',num2str(i,fspec));
    [vol,~,~,~,~,~,~,~,outfile] = cts_model_atomic(sz,targets,pmod,'suffix',suf); % need reframed params
    % outfile needs to be .atom rather than .mat
    [path,file,ext] = fileparts(outfile); outfile = fullfile(path,append(file,'.atom.mat')); % fix for atomic file
    for j=1:numel(variters)
        psim.(varname) = variters(j); % get iteration value for variable
        sufsim = append('var_',varname,'_',string(variters(j)));
        [~,~,atlas] = cts_simulate_atomic(outfile,psim,'suffix',sufsim);
    end
end


%% density

% more complicated, needs lots of models rather than simulation replicates
% no separate test dataset, just N replicates of each combination of layers?
% need at least 1 separate fully cellular versions as a potential test against all mixtures
targets = {'ribo__ribo__4ug0_4v6x.group.mat'}; % actual target class separate from others
% or MTs and also use real dataset? not poorly aligned one though
glob = 1;
fil = 1;
distract = 1;