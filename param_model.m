function param = param_model(in)
%param manager for cts_model

arguments
    in
end
%make preliminary list of all the parameters

%input volume
%pixelsize
%
%density, iterations
%constraints
%beads, grid, membrane
%ice

%figure out the bounds of each and different use cases

%what are non-params? need suffix opt separately, that's not a model parameter.

param = in;

end