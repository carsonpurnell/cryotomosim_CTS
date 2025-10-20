
fname = 'ribo_4a';
pmod = param_model(4,'iters',200,'mem',0,'grid',0);
[cts,~,~,~,~,~,~,~,~,outfile] = cts_model_atomic([800,800,120],1,pmod,'suffix',append(fname,'_9'),'dynamotable',1);
[path,name,ext] = fileparts(outfile);
outfile = fullfile(path,append(name,'.atom.mat')); %bake into sim function?

%%
% 80/2 is best, 80/1 has OOM on laptop
psim = param_simulate('pix',4,'defocus',-2,'dose',800,'raddamage',0,'tilt',[-80,2,80]);
cts_simulate_atomic(outfile,psim,'suffix','test_ideal_80-2_def-2');