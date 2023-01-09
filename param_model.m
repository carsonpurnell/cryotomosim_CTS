function param = param_model(param)
%param manager for cts_model

arguments
    %maybe GUI input dlg
    %easy way to navigate to stored params? modifiable defaults?
    
    param.vol
    param.pix
    
    param.targets = 'gui' %soft required?
    param.distractors = 'none'
    
    param.density = 0.4 %if moving to target loop, this needs to be able to be a vector
    param.iters = [] %auto calculate if not given, would also need to be vectorable
    param.constraint string {mustBeMember(param.constraint,{'none','box','tube','sides'})} = 'sides'
    %change constraint to a more flexible x/y/z for different oriented tube/walls?
    
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