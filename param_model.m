function param = param_model(param)
%param manager for cts_model

arguments
    %maybe GUI input dlg
    %easy way to navigate to stored params? modifiable defaults?
    
    %required for function params
    %param.vol (:,:,:) = zeros(400,400,50) %must be input?
    %param.pix (1,1) %must be input?
    
    %input particles params
    param.layers = 1 %soft required?
    %param.distract = 'none'
    
    %model run and limitations params
    param.density = 0.4 %if moving to target loop, this needs to be able to be a vector
    param.iters = [] %auto calculate if not given, would also need to be vectorable
    param.constraint string {mustBeMember(param.constraint,{'none','box','tube','sides'})} = 'sides'
    %change constraint to a more flexible x/y/z for different oriented tube/walls?
    
    %objects and contents params
    param.grid = [15 2000]
    param.mem = 0 %need subcomponents
    param.beads = 0
    param.ice = 1 %need more control. also could do with surface ice contamination, and more roughness
    
end
%{
if isempty(param.iters) || param.iters==0 %compute iters if not provided
    param.iters = round(param.pix*sqrt(numel(param.vol))/30);
end
%}
%make preliminary list of all the parameters

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