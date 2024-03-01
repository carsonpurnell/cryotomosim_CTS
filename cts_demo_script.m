%% get some pdb/cif files
%
%download stuff
%command for downloading pdbs?
% help FUNCTIONNAME provides more details, especially for more input details and options

%% load model parameters
% sets model generation parameters in the 'modparam' variable. uses an input dialog for all parameters
modparam = param_model('gui');
% example of directly declaring parameters, not using the GUI dialog
%modparam = param_model(12,'layers',2,'density',[0.3,0.5]);
%pixel size of 12(required in place of 'gui'), then layers/density are name-value arguments

%% generate model
% generates a model with the specified parameters. param_model can replace the modparam variable, when
% it is convinient - but the variable is usually easier (no need to load inputs multiple times)
modsize = zeros(400,300,50); % determines output model size, in pixels
[cts] = cts_model(modsize,modparam,'suffix','ctsdemo'); %output file names will be ctsdemo.ext

%% load simulation parameters
% sets tomogram simulation parameters, similarly to the model parameter function.
simparam = param_simulate('gui');
% similar non-GUI variable assertions - but param_simulate has no require inputs, all have defaults
%simparam = param_simulate('defocus',-4,'tilt',-50:4:50,'tiltscheme','symmetric');

%% run tilt simulation and reconstruction
% run a simulation, using a GUI to select the model. Either the .mat or .mrc from model generation will
%work for simulation, but only the .mat will be able to generate an atlas for the sim
[detected,conv] = cts_simulate('gui',simparam,'suffix','dosezero');
%example of skipping param_simulate, and providing parameters directly - inside the curly braces
%[detected,conv] = cts_simulate('gui',{'defocus',-5,'dose',0,'raddamage',1},'suffix','dosezero');
