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

%% SNR
n = 3; % number of training samples, one extra will be created as independent test data
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


%% pixel size



%% density