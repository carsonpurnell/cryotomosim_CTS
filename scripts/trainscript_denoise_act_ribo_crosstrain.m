fname = 'cross_5a';
ribo = {'ribosome__4ug0.mat'};
actin  = {'actin__6t1y_13x1.mat','actin__6t1y_13x2.mat'};

pmod_r = param_model(5,'layers',ribo,'iters',600,'mem',0,'grid',0);
pmod_a = param_model(5,'layers',actin,'iters',600,'mem',0,'grid',0);
sz = [500,500,100];
psim = param_simulate('pix',5,'defocus',-4,'dose',100,'raddamage',1,'tilt',[-60,2,60]);
psim_i = param_simulate('pix',5,'defocus',-3,'dose',600,'raddamage',0,'tilt',[-84,2,84]);

%% ribo
[cts,~,~,~,~,~,~,~,~,outfile] = cts_model_atomic(sz,1,pmod_r,'suffix',append(fname,'_ribo3'));
[path,name,ext] = fileparts(outfile);
outfile = fullfile(path,append(name,'.atom.mat')); %bake into sim function?

cts_simulate_atomic(outfile,psim,'suffix','ribo_sim3');
cts_simulate_atomic(outfile,psim_i,'suffix','ribo_ideal3');

%% actin
[cts,~,~,~,~,~,~,~,~,outfile] = cts_model_atomic(sz,1,pmod_a,'suffix',append(fname,'_actin3'));
[path,name,ext] = fileparts(outfile);
outfile = fullfile(path,append(name,'.atom.mat')); %bake into sim function?

cts_simulate_atomic(outfile,psim,'suffix','actin_sim3');
cts_simulate_atomic(outfile,psim_i,'suffix','actin_ideal3');