
fname = 'ribo_4a';
pmod = param_model(4,'iters',200,'mem',0,'grid',0);
[cts,~,~,~,~,~,~,~,~,outfile] = cts_model_atomic([800,800,120],1,pmod,'suffix',append(fname,'_8'),'dynamotable',1);

%%
psim = param_simulate('pix',4,'defocus',-3,'dose',600,'raddamage',0.002,'tilt',[-80,2,80]);
cts_simulate_atomic(oufile,psim,'suffix','test_ideal_80-2');