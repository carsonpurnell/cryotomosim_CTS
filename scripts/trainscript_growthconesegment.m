%rng(0) % static random number generator for replication

% hardcode input files
% preproc into layers, rearrange as needed rather than reload
targets = {'ribo__ribo__4ug0_4v6x.group.mat',... % two ribo variations
    'cofilactin__cofilactin__actin__actin__x3-x4_long.bundle.pdb'... % two lengths each
    'MT__6o2tx3.mat',... % long MT stretch
    'tric__tric__6nra-open_7lum-closed.group.mat'}; % two tric variations
%anything else?

% need all members, randomly add 1-2 replicates to shift occupancies?

distractors = {'act1-A2.distract.mat','1tub_Distractor.mat'};
% pool of distractors to pull random sample from to use as second/third? layer in modeling
%act/tub monomers, tiny mem decorations, 

n = 3; % number of different simulations
ptable = table; % initialize table of parameters, one row per run

% modeling params
ptable.pix(1:n) = 12+round(rand(n,1)*20)/10; %12-14 angstrom size/scale variation
ptable.thick(1:n) = 40+round(20*rand(n,1)); % 45-55 pixel thickness for SNR/orientation variation
ptable.iters(1:n) = 100+10*round(50*rand(n,1)); % 100-600 (increment 10) iterations
ptable.mem(1:n) = randi(12,n,1); %1-12 membranes
%randomized beads?

% simulation params
ptable.dose(1:n) = 70+round(80*rand(n,1)); %70-150 dose, uniform distribution
ptable.defocus(1:n) = -4-round(10*rand(n,1))/5; % -4 to -6 defocus
% radiation, tilting?

% z thickness from thin plane to thick low SNR (don't have variable layer thickness yet though)
% small range of xy sizes, mostly to have differently oriented rectangles? or extraneous? complex reporting
% same density cap for all layers for simplicity
% 2 target layers of fixed components + distractor layer of randomized grit

digits = numel(num2str(n));
fspec = append('%0',num2str(digits),'i'); %formatspec for suffixes
for i=1:n
    vol = zeros([400,400,ptable.thick(i)]);
    distix = logical(randi(2,1,numel(distractors))-1);
    if all(distix==false), distix(randi(numel(distix)))=1; end %ensure at least one distractor
    ptable.distractors(i,1:numel(distix)) = single(distix);
    linput = {targets,distractors(distix)};
    %tl = helper_input(linput,ptable.pix(i));
    dparam = param_model(ptable.pix(i),'layers',distractors(distix));
    modelparam = param_model(ptable.pix(i),'layers',targets,'iters',ptable.iters(i),...
        'mem',ptable.mem(i),'beads',10);
    modelparam.layers{2} = dparam.layers{1}; %add distractors into layer
    modelparam.iters(2) = modelparam.iters(1); modelparam.density(2) = modelparam.density(1);
    % distractor iters too?
    
    suf = append('batch_',num2str(i,fspec));
    [cts,outfile] = cts_model(vol,modelparam,'suffix',suf);
    
    %atomic or volumetric model? atomic membrane proteins still not ready
    simparam = param_simulate('dose',ptable.dose(i),'defocus',ptable.defocus(i));
    [~,~,~,atlas] = cts_simulate(outfile,simparam,'suffix','bulksim');
    ptable.classes(i) = numel(unique(atlas)); % class count to ensure full coverage
end
disp(ptable) % display table of parameters used

