
% use help COMMAND to see a more extensive description of individual commands and more options
% the model/sim functions and parameters have MANY options to cover different imaging and samples


%% get some pdb/cif structure files
% I suggest starting with large proteins with distinct shapes for clarity
% start with a pixel size of 6-12A for fast runtimes while testing. Small pixel sizes get much slower

%% load model parameters
% sets model generation parameters in the 'modparam' variable. uses an input dialog for all parameters
modparam = param_model('gui');
% example of directly declaring parameters, not using the GUI dialog
%modparam = param_model(12,'layers',2,'density',[0.3,0.5]);
% pixel size of 12(required in place of 'gui'), then layers/density are name-value arguments

%% generate model
% generates a model with the specified parameters. param_model can replace the modparam variable, when
% it is convinient - but the variable is usually easier (no need to load inputs multiple times)
modsize = zeros(400,300,50); % determines output model size, in pixels
[cts] = cts_model(modsize,modparam,'suffix','_ctsdemo'); %output file names will end with _ctsdemo.ext

%% load simulation parameters
% sets tomogram simulation parameters, similarly to the model parameter function.
simparam = param_simulate('gui');
% similar non-GUI variable assertions - but param_simulate has no require inputs, all have defaults
%simparam = param_simulate('defocus',-4,'tilt',-50:4:50,'tiltscheme','symmetric');

%% run tilt simulation and reconstruction
% generates s simulation in a subfolder under the model. 'gui' uses a file browser gui to select the model
% the model can be either the .mrc or .mat created earlier but only .mat will be able to generate an atlas.
[detected,conv] = cts_simulate('gui',simparam,'suffix','demosim');
% example of skipping param_simulate, and providing parameters directly - inside the curly braces
%[detected,conv] = cts_simulate('gui',{'defocus',-5,'dose',0,'raddamage',1},'suffix','dosezero');

%% review data
% all generated data is saved in the /tomosim directory, which should be in the home folder.
% the model folder is named according to datetime, the input particles, and pixel size
% simulation folders are generated within the model folder, named with the dose and 
% if using imod, you can navigate to the folder and easily view all outputs with 3dmod *.mrc.

