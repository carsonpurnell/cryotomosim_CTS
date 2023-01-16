function param = param_model(pix,param)
%param = param_model(pix,param)
%param manager for cts_model
%
%
%
%

arguments
    %maybe GUI input dlg
    %easy way to navigate to stored params? modifiable defaults?
    
    pix (1,1) %required input
    
    %required for function params
    %param.vol (:,:,:) = zeros(400,400,50) %must be input?
    %param.pix (1,1) %must be input?
    
    %input particles params
    param.layers = 1 %soft required?
    %param.distract = 'none'
    
    %model run and limitations params
    param.density = 0.4 %if moving to target loop, this needs to be able to be a vector
    param.iters = 0 %auto calculate if not given, would also need to be vectorable
    param.constraint string {mustBeMember(param.constraint,{'none','box','tube','sides'})} = 'sides'
    %change constraint to a more flexible x/y/z for different oriented tube/walls?
    
    %objects and contents params
    param.grid = [15 2000]
    param.mem = 0 %need subcomponents
    param.beads = 0
    param.ice = 1 %need more control. also could do with surface ice contamination, and more roughness
    
end
%{
if ~isfield(param,'pix')
    error('A pixel size is required as input')
end
%}
param.pix = pix;

layers = cell(1,param.layers);
iters = zeros(1,numel(param.layers));
for i=1:param.layers %loop through layers to load particles and assign iterations/density vectors
    fprintf('Loading layer %i structures \n',i)
    %need a check or something for when a layer is already parsed files? or just let input return by itself?
    %more likely need to be able to load a saved list of layers
    layers{i} = helper_input('gui',param.pix); %load layer
    param.density(i) = param.density(min(i,end));
    iters(i) = param.iters( min(i,numel(param.iters)) );
    if isempty(iters(i)) || iters(i)==0
        iters(i) = 2000*param.density(i);
    end
end
param.layers = layers;
param.iters = iters;

%input volume
%pixelsize
%targets, distractors - change to just particle layers?
%density, iterations
%constraints
%beads, grid, membrane
%ice

%figure out the bounds of each and different use cases

%what are non-params? need suffix opt separately, that's not a model parameter.

end