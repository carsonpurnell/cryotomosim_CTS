
fname = 'ribo_4a';
pmod = param_model(4,'iters',200,'mem',0,'grid',0);
[cts,~,~,~,~,~,~,~,~,outfile] = cts_model_atomic([800,800,120],1,pmod,'suffix',append(fname,'_8'),'dynamotable',1);
[path,name,ext] = fileparts(outfile);
outfile = fullfile(path,append(name,'.atom.mat')); %bake into sim function?

%%
psim = param_simulate('pix',4,'defocus',-3,'dose',300,'raddamage',0.01,'tilt',[-40,1,40]);
cts_simulate_atomic(outfile,psim,'suffix','test_ideal_40-1');