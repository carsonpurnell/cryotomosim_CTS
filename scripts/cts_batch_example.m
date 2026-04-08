% % example script for using CTS batch mode to generate numerous samples

n = 5; % number to mod-sim
%pix = 8.5; % currently fixed, should be easy to implement as variable
sz = [400,400,50]; % side lengths
% separate vector or second row to indicate variation in model size?
batchname = 'ATPs_mem';
targs = {'ATPS__flip.6j5i.membrane.cif','tubulin__1tub.distract.mat','act1-A2.distract.mat'...
    '1trv_thioredoxin.distract.pdb','ribo__ribo__4ug0_4v6x.group.mat','7b5s.distract.mat'};

%targs = {'memtest__8h9u.membrane.cif','ATPS__flip.6j5i.membrane.cif'}; batchname = 'ATPs2';
% calmodulin superres attempt?
%targs = {'tri__9K8D.cif'}; batchname = 'tri_test';

batchmod = param_batch(n,'pix',[8,9],'layers',{targs},'iters',[200,1000],'mem',[2,12]);
batchsim = param_batch(n,'dose',[60,150],'defocus',[-3,-5],'scatter',[0.5,1.5],'tilt',-60:3:60);

ideal = param_simulate('dose',500,'ice',0,'defocus',-3,'raddamage',0,'scatter',0.5,'tilt',-80:2:80);
%ideal = 0;

cts_batch(sz,batchmod,batchsim,'method','atom','batchname',batchname,'ideal',ideal);
