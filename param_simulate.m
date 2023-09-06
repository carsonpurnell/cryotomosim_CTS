function [param] = param_simulate(guiinput,param)
%[param] = param_simulate(guiinput,param)
%records microscope/image details for CTS simulation
%  Inputs
%
%guiinput        default 0
%if guiinput == 'gui', then a GUI dialogue will be shown for manual inputs. The GUI will use any paramter
%values that were supplied, and otherwise have defaults.
%
%name-value pair inputs (param):
%voltage -  in KeV
%aberration - in mm
%sigma - arbitrary value of the CTF envelope function
%defocus - in nm
%tilt - [min increment max], min:increment:max, or otherwise a list of tilt angles in the simulation
%dose - in e/A^2. one value is distributed, 0 skips dose simulation, else must be a vector matching tilts
%tiltscheme - 'symmetric' or otherwise a tilt angle that splits the half-tilt segments
%      tilt acquisition proceeds from the split angle to positive, then the split towards negative
%pix overrides a pixel size from an input model
%tiltax determines the tilt axis, input 'X' else 'Y'. still some janky rotations
%raddamage - arbitrary scale for radiation damage imparted to tilts (default 1)
%scatter - scales inelastic(lossy) electron scattering (nonlinear, scales distance relative to IMFP)
%ctfoverlap sets the degree of overlap between CTF computation strips (default 2), 0 skips CTF simulation
%tilterr scales tilt angle randomization from 0 to N, with 1 the error range equals the tilt increment

arguments
    guiinput = 0
    
    %microscope parameters
    param.voltage {mustBePositive} = 3e2
    param.aberration {mustBePositive} = 2.7
    param.sigma {mustBePositive} = 0.9
    
    %scan/session parameters
    param.defocus = -5
    param.tilt = [-60 2 60]
    param.dose = 60 %0 dose skips dose sim
    param.tiltscheme = 0 %flip angle for bidirectional, otherwise 'symmetric'
    param.pix = 0 %0 to not override input
    param.tiltax = 'Y' %using X now works but generates super thick tomograms
    
    %simulation/computational parameters
    param.raddamage = 1;
    param.scatter = 1;
    param.ctfoverlap {mustBeNonnegative,mustBeInteger} = 2 %if 0 skip ctf
    % Q factor for CTF convolution?
    %randomization to tilt angles?
    param.tilterr = 0; %0 no randomization, 1 range==tilt increment
end

if strcmp(guiinput,'gui') %basic GUI for manual input of values
    
    prompt = {'Microscope voltage (kV)',...
        'Microscope spherical aberration (mm)',...
        'Microscope Sigma value for CTF: ',...
        'Scan focus (negative for defocus, nm)',...
        'Tilt angles [min increment max] (single projection with same min/max and nonzero increment)',...
        'Scan dose, in e/A^2 (for dose weighting, enter weights in tilt order)',...
        'Tilt scheme ("symmetric" or otherwise the flip angle)',...
        'Pixel size (leave at 0 unless overriding)',...
        'tilt axis, X or Y (X kind of wonky)',...
        'scale of radiation damage (0=none, 1 = standard)',...
        'scale of inelastic (lossy) electron scattering',...
        'Processing: extent of defocus overlap for CTF calculation'};
    ptitle = 'Microscope and imaging parameters';
    
    default = struct2cell(param); %get all the default values back out into a list
    for i=1:numel(default)
        default{i} = num2str(default{i});
    end
    
    p = inputdlg(prompt,ptitle,[1 40],default); %dialogue input happens
    
    fn = fieldnames(param);
    for i=1:numel(fn) %extract the input values back into the param struct
        param.(fn{i}) = str2double(p{i});
    end
    
    param.dose = str2double(split(p{6}))'; %should work regardless of list or single value
    if ~strcmp(p{7},'symmetric'), p{7}=str2double(p{7}); end %if not symmetric, get the number
    param.tiltscheme = p{7};
    param.tiltax = p{9};
    
    param.tilt = str2double(split(p{5}))'; %extract the list of potential tiltangles
end

if numel(param.tilt)==3 %if 3 numbers are input try to resolve them as a vector
    tilt = param.tilt(1):param.tilt(2):param.tilt(3);
    err = param.tilt(2)/2; %calculate the base error size
    if ~isempty(tilt), param.tilt=tilt; end
    if param.tilterr~=0 %randomly generate angle errors and add to tiltangles
        err = (-err*1+rand(1,numel(tilt))*err*2);
        param.tilterr = err*param.tilterr;
        %param.tilt = param.tilt+param.tilterr;
    end
end

if numel(param.dose)~=numel(param.tilt) && numel(param.dose)~=1 %check that dose works with the tilts
    warning('dose and tilt are incompatible sizes - they must be equal length or there must be a single dose') 
end

end