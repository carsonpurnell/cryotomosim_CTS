
fname = 'ribo_4a';
targs = {'ribosome__4ug0.mat'};
pix = 4;
pmod = param_model(pix,'layers',targs,'iters',200,'mem',0,'grid',0);
psim1 = param_simulate('pix',pix,'defocus',-4,'dose',60,'raddamage',2,'tilterr',1,'tilt',[-60,2,60]);
psim2 = param_simulate('pix',pix,'defocus',-4,'dose',120,'raddamage',0.5,'tilterr',0.5,'tilt',[-60,2,60]);
psim3 = param_simulate('pix',pix,'defocus',-3,'dose',300,'raddamage',0.1,'scatter',0,'tilterr',0,'tilt',[-60,2,60]);
pideal = param_simulate('pix',pix,'defocus',-2,'dose',500,'raddamage',0,'scatter',0,'tilt',[-80,2,80]);
% not enough difference - bad ice implementation during tilts eroding signal?

for i=1:4
[cts,~,~,~,~,~,~,~,~,outfile] = cts_model_atomic([800,800,120],1,pmod,'suffix',append(fname,'_',string(i)),'dynamotable',1);
[path,name,ext] = fileparts(outfile);
outfile = fullfile(path,append(name,'.atom.mat')); %bake into sim function?

cts_simulate_atomic(outfile,psim1,'suffix',append('s1_',string(i)));
cts_simulate_atomic(outfile,psim2,'suffix',append('s2_',string(i)));
cts_simulate_atomic(outfile,psim3,'suffix',append('s3_',string(i)));
cts_simulate_atomic(outfile,pideal,'suffix',append('id_',string(i)));
end

% how many runs, and what levels of denoising? simulate at poor, mid, high quality and non-extreme ideal?
% might show extreme background flattening happening only at certain relative noise of original