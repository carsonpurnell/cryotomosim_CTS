%% model params
pix = 10;
pmod = param_model(pix,'iters',1000,'density',0.8,'layers',4,'mem',20);
%% run model
[cts] = cts_model(zeros(500,400,60),pmod,'suffix','cytoall_distract');

%% get the model file - manual!
modfile = util_loadfiles;
%modfile = 'C:\cygwin\home\exels\tomosim\model_2023-10-30t10.01_MT_mip_cofilactin_actin_tric_ribo_pixelsize_10\MT_mip_cofilactin_actin_tric_riboBASENAME.mat';
%modfile = 'C:\cygwin\home\exels\tomosim\atomicmods_manual\radtestingmix_8a.mrc';%camk2_radtest.mrc';
%modfile = 'C:\cygwin\home\exels\tomosim\model_2023-10-05t11.04\ribo_cofilactin_actinLTEST.mat';

%% sim params
psim = param_simulate('defocus',-5,'raddamage',0); %initial values for simulation - currently default
lvar = 'dose'; % name of variable to vary (named in the psim struct)
entries = [10,50,100,500,1000,5000]; % vector of values to simulate at
%% loop modify sim params and run the sim
for i=1:numel(entries)
    suf = append('tiltregressvol','_',num2str(lvar),'_',num2str(entries(i)));
    psim.(lvar) = entries(i); %replace variable of interest with target value
    [detected,conv] = cts_simulate(modfile,psim,'suffix',suf); % run the simulation instance
end