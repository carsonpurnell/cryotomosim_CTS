%rng(0) % static random number generator for replication

% hardcode input files
% preproc into layers, rearrange as needed rather than reload
targets = {'ribo__ribo__4ug0_4v6x.group.mat',... two ribo variations
    'cofilactin__cofilactin__actin__actin__x3-x4_long.bundle.pdb'...
    'MT__6o2tx3.mat',... % long MT stretch
    'tric__tric__6nra-open_7lum-closed.group.pdb'}; % two tric variations
%anything else?

%distractors = {}

% pool of distractors to pull random sample from
%act/tub monomers, tiny mem decorations, 

n = 3; % number of different simulations
ptable = table;
% randomize outside the loop, put random vals in a table for reference during generator loop
ptable.pix(1:n) = 12+round(rand(n,1)*20)/10; %12-14 angstrom size/scale variation
ptable.thick(1:n) = 45+round(10*rand(n,1)); % 45-55 pixel thickness for SNR/orientation variation
ptable.iters(1:n) = 200+10*round(50*rand(n,1)); % 200-700 (increment 10) iterations
%500+5*(round(50*(rand(n,1)-rand(n,1)))); %500+-250 (50 increments) iterations
ptable.mem(1:n) = 1+randi(6,n,1); %3-8 membranes

ptable.dose(1:n) = 80+round(80*rand(n,1)); %80-120 dose, uniform distribution
ptable.defocus(1:n) = -4-round(10*rand(n,1))/5; % -4 to -6 defocus

% randomize parameters - size, shape, densities, ordering (Mt/fils mainly), ice density
% z thickness from thin plane to thick low SNR (don't have variable layer thickness yet though)
% small range of xy sizes, mostly to have differently oriented rectangles? or extraneous? complex reporting
% same density for all layers for simplicity
% 2 target layers of fixed components + distractor layer of randomized grit
% run models

digits = numel(num2str(n));
fspec = append('%0',num2str(digits),'i');
for i=1:n
    vol = zeros([300,300,ptable.thick(i)]);
    linput = targets;
    %tl = helper_input(linput,ptable.pix(i));
    modelparam = param_model(ptable.pix(i),'layers',linput,'iters',ptable.iters(i),'mem',ptable.mem(i));
    
    suf = append('batch_',num2str(i,fspec));
    [cts,outfile] = cts_model(vol,modelparam,'suffix',suf);
    
    %atomic or volumetric model? atomic membrane proteins still not ready
    [~,~,~,atlas] = cts_simulate(outfile,{'dose',ptable.dose(i),'defocus',ptable.defocus(i)},'suffix','bulksim');
    ptable.classes(i) = numel(unique(atlas)); %count classes in atlas to make sure nothing got missed
end
disp(ptable) % display table of parameters used

%randomize params - slight pixel size, defocus, tilts? dose, radiation? tilt axis?
% run simulations

% check to make sure that the atlas/list has every member for DF DL training

% 10-20