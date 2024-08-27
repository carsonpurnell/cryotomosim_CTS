%rng(0) % static random number generator for replication

% hardcode input files
% preproc into layers, rearrange as needed rather than reload
targets = {'ribo__ribo__4ug0_4v6x.group.mat',... % two ribo variations
    'cofilactin__cofilactin__actin__actin__x3-x4_long.bundle.pdb'... % two lengths each
    'MT__6o2tx3.mat',... % long MT stretch
    'tric__tric__6nra-open_7lum-closed.group.pdb'}; % two tric variations
%anything else?

%distractors = {}
% pool of distractors to pull random sample from to use as second/third? layer in modeling
%act/tub monomers, tiny mem decorations, 

n = 3; % number of different simulations
ptable = table; % initialize table of parameters, one row per run

% modeling params
ptable.pix(1:n) = 12+round(rand(n,1)*20)/10; %12-14 angstrom size/scale variation
ptable.thick(1:n) = 45+round(10*rand(n,1)); % 45-55 pixel thickness for SNR/orientation variation
ptable.iters(1:n) = 200+10*round(40*rand(n,1)); % 200-600 (increment 10) iterations
%500+5*(round(50*(rand(n,1)-rand(n,1)))); %500+-250 (50 increments) iterations
ptable.mem(1:n) = randi(10,n,1); %1-6 membranes

% simulation params
ptable.dose(1:n) = 80+round(80*rand(n,1)); %80-120 dose, uniform distribution
ptable.defocus(1:n) = -4-round(10*rand(n,1))/5; % -4 to -6 defocus
% radiation, tilting?

% z thickness from thin plane to thick low SNR (don't have variable layer thickness yet though)
% small range of xy sizes, mostly to have differently oriented rectangles? or extraneous? complex reporting
% same density cap for all layers for simplicity
% 2 target layers of fixed components + distractor layer of randomized grit

digits = numel(num2str(n));
fspec = append('%0',num2str(digits),'i');
for i=1:n
    vol = zeros([400,400,ptable.thick(i)]);
    linput = targets;
    %tl = helper_input(linput,ptable.pix(i));
    modelparam = param_model(ptable.pix(i),'layers',linput,'iters',ptable.iters(i),'mem',ptable.mem(i));
    
    suf = append('batch_',num2str(i,fspec));
    [cts,outfile] = cts_model(vol,modelparam,'suffix',suf);
    
    %atomic or volumetric model? atomic membrane proteins still not ready
    simparam = param_simulate('dose',ptable.dose(i),'defocus',ptable.defocus(i));
    [~,~,~,atlas] = cts_simulate(outfile,simparam,'suffix','bulksim');
    ptable.classes(i) = numel(unique(atlas)); % class count to ensure full coverage
end
disp(ptable) % display table of parameters used

