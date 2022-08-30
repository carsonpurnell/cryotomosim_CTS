function [ts] = cts_model(particleset,vol,pix,opt)
%[ts] = tomosim_model(particleset,vol,pix,opt)
%generates model information for a single tomographic acquisition, stored in output struct ts
%works by iteratively placing input particles at random orientations in random locations without overlap
%
%Inputs
%
%particles      required, preferred input = 'gui'
%cell array of filenames to load as the target particles, .pdb or .mrc, processed by helper_input
%'gui' uses a GUI to select input files (uses uipickfiles from the FEX if possible)
%if not using 'gui', need to specify files on path or full file paths as strings in a cell array
%when loading, the class ID of an input is the leading string up to the first '__' or '.'
%different objects with the same ID will get placed in the same class of splitmodel
%
%vol            required, preferred input = zeros(x,y,z)
%3d array to fill with the model. standard is empty but accepts non-empty arrays (needs voids to work though)
%
%pix            required
%pixel size of the generated model and output .mrc, in angstroms
%
%opt - name-value pair arguments in the format ...'name',value... (or in 2021+ ...name=value...)
%density        default 0.4
%    cutoff for occupancy of model-filling steps. rarely encountered unless particles are highly compact
%constraint     default 'sides' (i.e. top and bottom of z)
%    must be 'sides','box','tube', or 'none'. determines how many sides have borders that prevent clipping
%distract       default 'none'
%    usage is identical to particles, places additional non-target objects after the initial target model
%beads          default [0 50], format [number radius1 radius2... radiusn]
%    number of beads to place, from a set of beads generated based on the input radii (default 50A)
%grid           default [0 0], author uses [15 2000]
%    [thickness radius] of carbon film and grid hole, in nm.
%ice            default 1
%     change to 0 to not generate vitreous ice in the model
%mem            default 0
%    1 to generate a super janktastic not at all realistic fascimile of a cell membrane
%graph          default 0
%    1 to have a plot continuously update with particle placement success/failure numbers
%
%Output
%a folder will be generated with a name of the input particles, with a .mrc of vol and a .mat of ts
%ts is a multilayer struct organizing several inputs, outputs, and intermediates
%vol has the final output model, containing all objects 
%pix retains the input pixel size
%model retains the separate classes of components - the isolated grid, the constraint borders, targets, or
%distractors
%targets is a struct array of each input particle grouping, including filename, id, volume, and type
%distractors is equivalent to targets, but for any input distractors
%splitmodel has fields of each target id, containing a volume of only those particles from the model

arguments
    particleset
    vol (:,:,:) double
    pix (1,1) double
    
    opt.density = 0.4
    opt.constraint string {mustBeMember(opt.constraint,{'none','box','tube','sides'})} = 'sides'
    opt.distract = 'none'
    opt.beads = 0 %new beads [number radius1 radius2 ... radiusn]
    %opt.beads = 0 %have beads use [number diamenter] and diam defaults to 10nm/100A.
    opt.grid = [15 2000] %[thick radius] both in nm [15 2000] our real grids
    opt.ice = 1 %0 to not add ice
    opt.mem = 0 %is a super jankfest that might not work right now
    opt.graph = 0
end

runtime = numel(vol)/60*1.2e-5;
fprintf('Estimated model generation time: %g minutes\n',runtime)
if runtime>30
    txt = input('Runtime is long, verify inputs. ctlr+C to end, or enter "proceed" to run anyway: ','s');
    if ~strcmp(txt,'proceed')
        ts = NaN; fprintf('Model generation declined, process aborted.\n')
        return
    end
end

%initialize the struct so the order is invariant and fill with input information
ts = struct('vol',vol,'pix',pix,'model',[],'particles',[],'splitmodel',[],'inputs',[]);
ts.inputs.pix = pix;
ts.inputs.density = opt.density; ts.inputs.constraint = opt.constraint;
ts.inputs.beads = opt.beads; ts.inputs.grid = opt.grid; ts.inputs.ice = opt.ice; 
%inputs.mem = opt.mem;

%load target particles and write to struct
if isstruct(particleset)
    ts.particles.targets = particleset;
else
    [ts.particles.targets] = helper_input(particleset,pix);
end

if opt.grid(1)~=0 % new carbon grid and hole generator
    fprintf('Generating carbon film ')
    [ts.vol] = gen_carbongrid(vol,pix,opt.grid);
    ts.model.grid = ts.vol; fprintf('\n')
end

%need a new mem gen that uses generated points between radii of two spheres
if opt.mem~=0 %membrane generation, likely broken with struct changes but old junk anyway
    mem = helper_membranegen(ts);
    ts.model.mem = mem; ts.vol = ts.vol+mem; %vol = ts.vol;
end

constraint = zeros(size(ts.vol)); %constraints are a big ugly mess right now
switch opt.constraint %write constraints to initial starting volume
    case 'none'
    case 'box' %intensity is ^2.3 to better match protein and prevent bad binarizations/overlap
        constraint(1:end,1:end,[1 end]) = pix^2.5; %constraint(1:end,1:end,end) = 1; %z end panes
        constraint(1:end,[1 end],1:end) = pix^2.5; %constraint(1:end,end,1:end) = 1; %y end panes
        constraint([1 end],1:end,1:end) = pix^2.5; %constraint(end,1:end,1:end) = 1; %x end panes
        disp('Warning: with a complete box, some particles may be impossible to place')
        ts.model.constraintbox = constraint;
    case 'tube'
        constraint(1:end,1:end,[1 end]) = pix^2.5; %constraint(1:end,1:end,end) = 1; %z end panes
        constraint(1:end,[1 end],1:end) = pix^2.5; %constraint(1:end,end,1:end) = 1; %y end panes
        ts.model.constrainttube = constraint;
    case 'sides'
        constraint(1:end,1:end,[1 end]) = pix^2.5; %constraint(1:end,1:end,end) = 1; %z end panes
        ts.model.constraintsides = constraint;
end

%generate model and add (in case input vol had stuff in it)
iters = round(ts.pix(1)*sqrt(numel(ts.vol))/30); %modeling iters, maybe simplify
[fill, ts.splitmodel] = helper_randomfill(ts.vol+constraint,ts.particles.targets,iters,...
    opt.density,'type','target','graph',opt.graph); 
ts.vol = ts.vol+fill; 
ts.model.targets = fill;
ts.model.particles = ts.vol;

if ~strcmp(opt.distract,'none') %DISTRACTORS
[ts.particles.distractors] = helper_input(opt.distract,pix,3); %load distractor particle set

%generated distraction filler iterations and add to volume to generate the sample
iters = round( iters*sqrt(numel(ts.particles.distractors(1,:))) ); %distractor iters
[fill] = helper_randomfill(ts.vol+constraint,ts.particles.distractors,iters,opt.density,'type','distractor');
ts.model.distractors = fill;
ts.vol = fill + ts.vol; 
ts.model.particles = ts.vol;
end

if opt.beads~=0 %bead generation and placement block
    beadstrc = gen_beads(pix,opt.beads(2:end)); %external generation of varied beads
    ts.particles.beads = beadstrc;
    [fill, ~] = helper_randomfill(ts.vol+constraint,beadstrc,opt.beads(1),opt.density,'type','bead');
    ts.model.beads = fill;
    ts.vol = fill + ts.vol; 
    ts.model.particles = ts.vol;
end

if ~opt.ice==0 % vitreous ice generator, randomized molecular h2o throughout the volume
    fprintf('Generating vitreous ice')
    [iced, ice] = gen_ice(ts.vol,pix);
    ts.model.ice = ice; ts.vol = iced; fprintf('\n')
end

%folder and file generation stuff
time = string(datetime('now','Format','yyyy-MM-dd''t''HH.mm')); %timestamp
ident = strjoin(fieldnames(ts.splitmodel),'_'); %combine target names to one string
foldername = append('model_',time,'_',ident,'_pixelsize_',string(pix)); %combine info for folder name
%filename needs to use inputs like size and other options after particle info

%move to output directory in user/tomosim
cd(getenv('HOME')); if ~isfolder('tomosim'), mkdir tomosim; end, cd tomosim
mkdir(foldername); cd(foldername);
WriteMRC(ts.vol,ts.pix(1),append(ident,'.mrc'))
save(append(ident,'.mat'),'ts','-v7.3')

%output text file of input informations?

cd(userpath)
end