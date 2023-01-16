function [cts] = cts_model(vol,param,opt)
%[cts] = cts_model(vol,pix,param,opt)
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
    %targets
    vol (:,:,:) double
    %pix (1,1) double
    param = param_model
    opt.suffix = ''
    
    %{
    opt.density = 0.4
    opt.constraint string {mustBeMember(opt.constraint,{'none','box','tube','sides'})} = 'sides'
    opt.distract = 'none'
    opt.beads = 0 %new beads [number radius1 radius2 ... radiusn]
    opt.grid = [15 2000] %[thick radius] both in nm [15 2000] our real grids
    opt.ice = 1 %0 to not add ice
    opt.mem = 0 %is now working, mostly! value is now the number of vesicles randomly generated
    opt.iters = 0
    %}
    %opt.layers = 0
    opt.graph = 0
    %suffix or other indicator string
    %opt to save incremental models for each layer and component of model building
    %save the splitmodels in another file to reduce bloat in cts?
end
if iscell(param), param = param_model(param{:}); end
pix = param.pix;
%{
runtime = numel(vol)/60*1.2e-5; %for my laptop, doesn't really apply to anything else
%need to compute by iterations too, vol alone not that relevant
fprintf('Estimated model generation time with hamster laptop: %g minutes\n',runtime)
if runtime>30 %if >30 mins force manual input start
    txt = input('Runtime is long, verify inputs. ctlr+C to end, or enter "proceed" to run anyway: ','s');
    if ~strcmp(txt,'proceed')
        cts = NaN; fprintf('Model generation declined, process aborted.\n')
        return
    end
end
%}

%initialize the struct so the order is invariant and fill with input information
cts = struct('vol',vol,'pix',pix,'model',[],'particles',[],'splitmodel',[]);%,'inputs',[]);
% cts.inputs.pix = pix;
% cts.inputs.density = opt.density; cts.inputs.constraint = opt.constraint;
% cts.inputs.beads = opt.beads; cts.inputs.grid = opt.grid; cts.inputs.mem = opt.mem;
% cts.inputs.ice = opt.ice; 

%{
if isempty(param.iters) || param.iters==0
    param.iters = round(cts.pix(1)*sqrt(numel(cts.vol))/30); %modeling iters, maybe simplify
else
    param.iters = param.iters;
end
%}

% block placeholder for loop input of particle set layers
%{
layers = cell(1,param.layers);
for i=1:param.layers
    fprintf('Loading layer %i structures \n',i)
    %need a check or something for when a layer is already parsed files? or just let input return by itself?
    %more likely need to be able to load a saved list of layers
    layers{i} = helper_input('gui',pix); %load layer
    %param.iters(i) = param.iters(min(i,end));
    %param.density(i) = param.density(min(i,end));
end
%}
%disp(param)
% block placeholder for loop input of particle set layers

%load input targets
%[cts.particles.targets] = helper_input(param.targets,pix); %load target particles
cts.particles.targets = param.layers;%{1};
% if isempty(param.iters) || param.iters==0
%     param.iters = round(cts.pix(1)*sqrt(numel(cts.vol))/30); %modeling iters, maybe simplify
% else
%     param.iters = param.iters;
% end


if param.grid(1)~=0 % new carbon grid and hole generator
    fprintf('Generating carbon film ')
    [cts.model.grid] = gen_carbongrid(vol,pix,param.grid);
    cts.vol = cts.model.grid+cts.vol; fprintf('   complete \n')
end

if param.mem~=0 %new membrane gen, makes spherical vesicles and places randomly
    fprintf('Generating vesicular membranes ')
    [cts.model.mem,count,~,vescen,vesvol] = gen_vesicle(cts.vol,round(param.mem),pix);
    cts.vol = cts.model.mem+cts.vol;
    fprintf('   complete,  %i placed, %i failed \n',count.s,count.f)
else
    vescen = 0; vesvol = 0;
end

constraint = zeros(size(cts.vol)); %constraints are a big ugly mess right now
switch param.constraint %write constraints to initial starting volume
    case 'none'
    case 'box' %intensity is ^2.3 to better match protein and prevent bad binarizations/overlap
        constraint(1:end,1:end,[1 end]) = pix^2.5; %constraint(1:end,1:end,end) = 1; %z end panes
        constraint(1:end,[1 end],1:end) = pix^2.5; %constraint(1:end,end,1:end) = 1; %y end panes
        constraint([1 end],1:end,1:end) = pix^2.5; %constraint(end,1:end,1:end) = 1; %x end panes
        disp('Warning: with a complete box, some particles may be impossible to place')
        %ts.model.constraintbox = constraint;
    case 'tube'
        constraint(1:end,1:end,[1 end]) = pix^2.5; %constraint(1:end,1:end,end) = 1; %z end panes
        constraint(1:end,[1 end],1:end) = pix^2.5; %constraint(1:end,end,1:end) = 1; %y end panes
        %ts.model.constrainttube = constraint;
    case 'sides'
        constraint(1:end,1:end,[1 end]) = pix^2.5; %constraint(1:end,1:end,end) = 1; %z end panes
        %ts.model.constraintsides = constraint;
end


%generate model and add (in case input vol had stuff in it)
[cts.model.targets, cts.splitmodel] = helper_randomfill(cts.vol+constraint,cts.particles.targets,param.iters,...
    vescen,vesvol,param.density,'type','target','graph',opt.graph); 
cts.vol = max(cts.vol,cts.model.targets); %to avoid overlap intensity between transmem and vesicle
%cts.vol = cts.vol+cts.model.targets; %old sum without overlap fix
cts.model.particles = cts.vol;

%change targets/distractors into a single repeating loop of any number of sets of particles?
%not sure how to implement a single opt to retrieve multiple sets of particles. loop to helper?
%also need to fetch particles before grid/membrane for ease of use

%{
if ~strcmp(param.distract,'none') %DISTRACTORS
[cts.particles.distractors] = helper_input(param.distract,pix); %load distractor particles

%generated distraction filler iterations and add to volume to generate the sample
iters = round( param.iters*sqrt(numel(cts.particles.distractors(1,:))) ); %distractor iters
[cts.model.distractors] = helper_randomfill(cts.vol+constraint,cts.particles.distractors,iters,...
    param.density,'type','distractor','graph',opt.graph);
cts.vol = max(cts.vol,cts.model.distractors); %fix for transmembrane overlaps
cts.model.particles = cts.vol;
end
%}

if param.beads~=0 %bead generation and placement block
    beadstrc = gen_beads(pix,param.beads(2:end)); %external generation of varied beads
    cts.particles.beads = beadstrc;
    [cts.model.beads] = helper_randomfill(cts.vol+constraint,beadstrc,param.beads(1),param.density,...
        'type','bead');
    cts.vol = cts.vol + cts.model.beads; 
    cts.model.particles = cts.vol;
end

if ~param.ice==0 % vitreous ice generator, randomized molecular h2o throughout the volume
    fprintf('Generating vitreous ice')
    [iced, ice] = gen_ice(cts.vol,pix);
    cts.model.ice = ice; cts.vol = iced; fprintf('   done \n')
end

%folder and file generation stuff
time = string(datetime('now','Format','yyyy-MM-dd''t''HH.mm')); %timestamp
ident = char(strjoin(fieldnames(cts.splitmodel),'_')); %combine target names to one string
if length(ident)>60, ident=ident(1:60); end %truncation check to prevent invalidly long filenames
foldername = append('model_',time,'_',ident,'_pixelsize_',string(pix)); %combine info for folder name

%move to output directory in user/tomosim
cd(getenv('HOME')); if ~isfolder('tomosim'), mkdir tomosim; end, cd tomosim
mkdir(foldername); cd(foldername);
WriteMRC(cts.vol,cts.pix(1),append(ident,opt.suffix,'.mrc'))
save(append(ident,opt.suffix,'.mat'),'cts','-v7.3')

%output text file of input informations?

cd(userpath)
end