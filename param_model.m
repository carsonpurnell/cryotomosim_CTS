function param = param_model(pix,param)
%param = param_model(pix,param)
%param manager for cts_model
%
%pix            required
%pixel size of the generated model and output .mrc, in angstroms
%
%param - name-value pair arguments in the format ...'name',value... (or in 2021+ ...name=value...)
%density        default 0.4
%    cutoff for occupancy of model-filling steps. rarely encountered unless particles are highly compact
%constraint     default '  &' - 1x3 char defining edge borders. uses helper_constraints
%beads          default [0 50], format [number radius1 radius2... radiusn]
%    number of beads to place, from a set of beads generated based on the input radii (default 50A)
%grid           default [0 0], author uses [15 2000]
%    [thickness radius] of carbon film and grid hole, in nm.
%ice            default 1
%     change to 0 to not generate vitreous ice in the model
%mem            default 0, otherwise an integer
%    1 to generate a super janktastic not at all realistic fascimile of a cell membrane
%filaments      default 0, otherwise 1
%    1 to generate flexible polymers from separate input structures - example cytoskeletal filaments
%graph          default 0
%    1 to have a plot continuously update with particle placement success/failure numbers
%
%
%particles      required, preferred input = 'gui'
%cell array of filenames to load as the target particles, .pdb or .mrc, processed by helper_input
%'gui' uses a GUI to select input files (uses uipickfiles from the FEX if possible)
%if not using 'gui', need to specify files on path or full file paths as strings in a cell array
%when loading, the class ID of an input is the leading string up to the first '__' or '.'
%different objects with the same ID will get placed in the same class of splitmodel
%
%pix is required: defines the pixel size of generated structure volumes
%
%other arguments are name-value pairs:
%layers - number of different particle layers that will be modeled
%density - vector of maximum allowed densities. if there's more layers it will use the last provided density
%iters - vector of iterations per layer. like density, it will use the last provided value if shot
%constraint - sides/box/tube of walls to contain particles
%grid/mem/bead/ice - same as old cts_model, help still should be there
%
% see also helper_constraints, helper_input, helper_pdb2vol
% gen_carbongrid, gen_ice, and gen_beads
arguments
    %easy way to navigate to stored params? modifiable defaults?
    %which params need more thorough validation?
    pix = 'gui' %(1,1) %required input
    param.pix = 'REQUIRED'
    
    param.layers = 1 %number of different particle layers to load and place
    
    %model run and limitations params
    param.density = 0.4 %if moving to target loop, this needs to be able to be a vector
    param.iters = 1000 %auto calculate if not given, would also need to be vectorable
    param.constraint char = '  &' %use & +- to add borders to different axes
    
    %objects and contents params
    param.grid = [150 1e4]
    param.mem = 0 %need subcomponents
    param.filaments = 0
    param.beads = 0
    param.ice = 1 %need more control. also could do with surface ice contamination, and more roughness
end
if strcmp(pix,'gui')
    prompt = {'Pixel Size, in A',...
        'Particle layers',...
        'Maximum density (up to # layers)',...
        'Placement iterations (up to # layers)',...
        'Border Constraints (3-char string for xyz axes)',...
        'Carbon grid (thickness radius, blank for no carbon)',...
        'Number of vesicles',...
        'Filaments (any ~=0 to use, barely-functional WIP)',...
        'Number of beads',...
        'Ice opacity'};
    ptitle = 'Model Generation Parameters';
    
    default = struct2cell(param); %get all the default values back out into a list
    for i=1:numel(default)
        default{i} = num2str(default{i});
    end
    
    p = inputdlg(prompt,ptitle,[1 40],default); %dialogue input happens
    
    fn = fieldnames(param);
    for i=1:numel(fn) %extract the input values back into the param struct
        param.(fn{i}) = str2double(p{i});
    end
else
    param.pix = pix;
end

param.mem = round(param.mem); 

layers = cell(1,param.layers);
iters = zeros(1,numel(param.layers));
for i=1:param.layers %loop through layers to load particles and assign iterations/density vectors
    fprintf('Loading layer %i structures \n',i)
    %need a check or something for when a layer is already parsed files? or just let input return by itself?
    layers{i} = helper_input('gui',param.pix); %load layer - how to deal with saved list of layers?
    param.density(i) = param.density(min(i,end));
    iters(i) = param.iters( min(i,numel(param.iters)) );
    if isempty(iters(i)) || iters(i)==0
        iters(i) = 2500*param.density(i);
    end
end
param.layers = layers; param.iters = iters;

if param.filaments~=0
    param.filaments = helper_filinput(pix,'gui');
end

end