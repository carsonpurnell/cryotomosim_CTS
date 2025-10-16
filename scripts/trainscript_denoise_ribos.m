
fname = 'ribo_4a';
pmod = param_model(4,'iters',200,'mem',0,'grid',0);
[cts,~,~,~,~,~,~,~,~,outfile] = cts_model_atomic([800,800,120],1,pmod,'suffix',append(fname,'_9'),'dynamotable',1);
[path,name,ext] = fileparts(outfile);
outfile = fullfile(path,append(name,'.atom.mat')); %bake into sim function?

%%
a = 80;
psim = param_simulate('pix',4,'defocus',-3,'dose',a*5,'raddamage',0,'tilt',[-a,1,a]);
cts_simulate_atomic(outfile,psim,'suffix','test_ideal_80-1');