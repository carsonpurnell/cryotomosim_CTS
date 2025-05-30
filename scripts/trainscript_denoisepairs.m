randset = rand*1e5; % starting seed for random number generation
rng(randset) % static random number generator for replication

n = 10; % number of different simulations
vs = 400; % side length of model volume

% hard set target list for first layer
targets = {'ribo__ribo__4ug0_4v6x.group.mat',... % two ribo variations
    'cofilactin__cofilactin__actin__actin__x3-x4_long.bundle.pdb'... % two lengths each
    'MT__6o2tx3.mat',... % long MT stretch
    'tric__tric__6nra-open_7lum-closed.group.mat'}; % two tric variations

% distractor pool
distractors = {'actin_monomer_2q0u.distract.mat',... % 1
    '1tub_tubulin_dimer.distract.mat',...
    '2cg9_HSP90.distract.mat',...
    '2vz6_CaMK2A.distract.pdb',...
    '7sgm_clathrin.distract.mat',...                 % 5
    'actinin_1sjj.distract.mat',...
    '7b5s_E3ligase.distract.mat',...
    '6lfm_gprotein.distract.membrane.mat',...
    'GABAar.distract.membrane.mat'};


ptable = table; % initialize table of parameters, one row per run

% modeling params
ptable.pix(1:n) = 5+randi(9,n,1); %6-14 angstrom size/scale variation
ptable.thick(1:n) = 40+round(40*rand(n,1)); % 40-80 pixel thickness for SNR/orientation variation
ptable.iters(1:n) = 100+10*round(70*rand(n,1)); % 100-800 (increment 10) iterations
ptable.mem(1:n) = (randi(4,n,1)-1)*4; %0-12 membranes
ptable.beads(1:n) = (randi(4,n,1)-1)*4; %0-12 beads

% simulation params
ptable.dose(1:n) = 80+1*round(80*rand(n,1)); %80-160 dose, uniform distribution
ptable.defocus(1:n) = -2-round(10*rand(n,1))/5-ptable.pix/5; % -2 to -4, minus pix/5
% radiation, tilting?
ptable.mill(1:n) = (randi(2,n,1)-1).*(rand(n,1)*2/5+0.8); % 0.8 to 1.2 mill, 50% of the time
% randomize use of milling surface artifacts

% z thickness from thin plane to thick low SNR (don't have variable layer thickness yet though)
% same density cap for all layers for simplicity

digits = numel(num2str(n));
fspec = append('%0',num2str(digits),'i'); %formatspec for suffixes
for i=1:n
    rng(randset+i);
    vol = zeros([vs,vs,ptable.thick(i)]);
    
    %distix = logical(randi(2,1,numel(distractors))-1);
    %if all(distix==false), distix(randi(numel(distix)))=1; end %ensure at least one distractor
    %ptable.distractors(i,1:numel(distix)) = single(distix);
    %linput = {targets,distractors(distix)};
    %tl = helper_input(linput,ptable.pix(i));
    
    dx = unique(randi(numel(distractors),1,2+randi(5)));
    ptable.distractors(i) = join(string(dx));
    dis = distractors(dx);
    dparam = param_model(ptable.pix(i),'layers',dis);
    modelparam = param_model(ptable.pix(i),'layers',targets,'iters',ptable.iters(i),...
        'mem',ptable.mem(i),'beads',ptable.beads(i));
    modelparam.layers{2} = dparam.layers{1}; %add distractors into layer
    modelparam.iters(2) = modelparam.iters(1)*2; modelparam.density(2) = modelparam.density(1);
    % distractor iters too?
    
    suf = append('denoise_batch_',num2str(i,fspec));
    rng(randset+i);
    [cts,outfile] = cts_model(vol,modelparam,'suffix',suf);
    
    %atomic or volumetric model? atomic membrane proteins still not ready
    simparam_noisy = param_simulate('dose',ptable.dose(i),'defocus',ptable.defocus(i),...
        'tilterr',1,'raddamage',1,'scatter',1,'mill',ptable.mill(i));
    rng(randset+i);
    sufsim = append(string(i),'_in');
    [~,~,~,atlas] = cts_simulate(outfile,simparam_noisy,'suffix',sufsim);
    
    simparam_quality = param_simulate('dose',ptable.dose(i)*0,'defocus',ptable.defocus(i)/2,...
        'raddamage',0,'tilt',-85:1:85,'ctfoverlap',2,'scatter',0,'mill',0);
    rng(randset+i);
    sufsim = append(string(i),'_out');
    [det,~,~,~] = cts_simulate(outfile,simparam_quality,'suffix',sufsim);
    % re-load recon for analysis, compute metrics like SSI between image pair?
    % also reorganize outputs better - into a struct of some sort? or save the mrc and dump?
end
disp(ptable) % display table of parameters used
