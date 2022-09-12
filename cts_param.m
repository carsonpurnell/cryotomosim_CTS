function [param] = cts_param(guiinput,param)
%WIP help
%if guiinput =='gui' an input GUI is used for manual input, it will load any values given over defaults
%voltage in KeV
%aberration in mm
%sigma arbitrary
%defocus in nm
%tilt as [min increment max] in command line or just the numbers in GUI, 3 numbers are resolved as
%   1:2:3 in matlab convention, otherwise the vector is the list of angles to acquire
%dose in e/A^2. one value is distributed across all tilts, 0 does not simulate electron detection, and a
%   vector with length==number of tilts is used for the dose for those tilts in order
%tiltscheme is 'symmetric' for dose-symmetric or a number that sets the start point for two-phase imaging
%   single-phase requires setting the value to either the min or max tilt angle, depending on direction
%(order is start point to positive/max, then start point towards negative)
%pix is used to override the real pixel size of a model/mrc, otherwise leave as 0
%tiltax sets tilt axis, 'X' or anything else becomes 'Y' - still janky
%raddamage scales radiation damage (arbitrary)
%scatter scales inelastic/lossy electron scattering (nonlinear, scales distance)
%ctfoverlap sets how much overlap between CTF strips (default 2), 0 skips CTF convolution entirely

arguments
    guiinput = 0
    
    %microscope parameters
    param.voltage {mustBePositive} = 3e2
    param.aberration {mustBePositive} = 2.7
    param.sigma {mustBePositive} = 1
    
    %scan/session parameters
    param.defocus = -5
    param.tilt = [-60 2 60]
    param.dose = 60 %do a check for dose==0 to not do electron detection at all/perfect detector
        %single = full scan, list==tilts in electrons for that scan
    param.tiltscheme = 0 %flip angle for bidirectional, otherwise 'symmetric'
    param.pix = 0 %0 to not override input
    param.tiltax = 'Y' %using X now works but generates super thick tomograms
    
    %simulation/computational parameters
    param.raddamage = 1;
    param.scatter = 1;
    param.ctfoverlap {mustBeNonnegative,mustBeInteger} = 2 %if 0 skip ctf
    %randomization to tilt angles?
    param.tilterr = 0;
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
    err = param.tilt(2)/2;
    if ~isempty(tilt), param.tilt=tilt; end
    if param.tilterr~=0
        err = (-err+rand(1,numel(tilt))*err);
        param.tilt = param.tilt+err*param.tilterr;
    end
end


if numel(param.dose)~=numel(param.tilt) && numel(param.dose)~=1 %check that dose works with the tilts
    warning('dose and tilt are incompatible sizes - they must be equal length or there must be a single dose') 
end

end