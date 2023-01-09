function param = param_model(param)
%param manager for cts_model

arguments
    %maybe GUI input dlg
    
    param.vol
    param.pix
    
    param.targets = 'gui' %soft required?
    param.distractors = 'none'
    
    param.density
    param.iters
    param.constraints
    
    param.grid
    param.mem %need subcomponents
    param.bead
    
    param.ice
    
end
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