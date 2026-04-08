% % example script for using CTS batch mode to generate numerous samples

%% set general batch inputs
n = 5; % number of total runs in the batch=
sz = [400,400,50]; % size of the samples, in pixels
batchname = 'ATPs_mem'; % suffix appended to model folder names
% list of structure files for modeling in all runs: 
targs = {'ATPS__flip.6j5i.membrane.cif','tubulin__1tub.distract.mat','act1-A2.distract.mat'...
    '1trv_thioredoxin.distract.pdb','ribo__ribo__4ug0_4v6x.group.mat','7b5s.distract.mat'};


%% generate the parameters for the batch runs
batchmod = param_batch(n,'pix',[8,9],'layers',{targs},'iters',[200,1000],'mem',[2,12]);
batchsim = param_batch(n,'dose',[60,150],'defocus',[-3,-5],'scatter',[0.5,1.5],'tilt',-60:3:60);

ideal = param_simulate('dose',500,'ice',0,'defocus',-3,'raddamage',0,'scatter',0.5,'tilt',-80:2:80);
%ideal = 0;

%% run the batch process
cts_batch(sz,batchmod,batchsim,'method','atom','batchname',batchname,'ideal',ideal);
