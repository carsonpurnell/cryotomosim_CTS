function [param] = cts_param(guiinput,param)
%WIP help
%if guiinput =='gui' an input GUI is used for manual input, it will load any values given over defaults
%voltage in KeV
%aberration in mm
%sigma arbitrary
%defocus in nm
%tilt as [min increment max] in command line or just the numbers in GUI
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
    param.tilt (1,3) = [-60 2 60]
    param.dose = 60 %do a check for dose==0 to not do electron detection at all/perfect detector
        %activate dose weighting with a list of length==tilts?
        %single = full scan, list==tilts in electrons for that scan
    param.tiltscheme = 0 %flip angle for bidirectional, otherwise 'symmetric'
    %img.symmetric = 0
    %param.doseweight = 0
    param.pix = 0 %0 to not override input
    param.tiltax = 'Y' %using X now works but generates super thick tomograms
    
    %simulation/computational parameters
    param.raddamage = 1;
    param.scatter = 1;
    param.ctfoverlap {mustBeNonnegative,mustBeInteger} = 2
    %randomization to tilt angles?
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
        %'Save these parameters as a file for later use? (1 for yes): '};
        %'Scan dose distribution, 1 for weighting --UNUSED--',...
    
    ptitle = 'Microscope and imaging parameters';
    
    default = struct2cell(param); %get all the default values back out into a list
    for i=1:numel(default)
        default{i} = num2str(default{i});
    end
    
    %{
    default = {num2str(microscope.voltage),num2str(microscope.aberration),num2str(microscope.sigma),...
        num2str(img.tilt),num2str(img.defocus),num2str(img.dose),...
        num2str(img.symmetric),num2str(img.doseweight),num2str(img.raddamage),...
        num2str(img.ctfoverlap),img.tiltax,num2str(img.pix)};
    %}
    p = inputdlg(prompt,ptitle,[1 40],default); %dialogue input happens
    
    fn = fieldnames(param);
    for i=1:numel(fn) %extract the input values back into the param struct
        param.(fn{i}) = str2double(p{i});
    end
    
    param.tilt = str2double(split(p{5}))';
    param.dose = str2double(split(p{6}))'; %should work regardless of list or single value
    if ~strcmp(p{7},'symmetric'), p{7}=str2double(p{7}); end %if not symmetric, get the number
    param.tiltscheme = p{7};
    param.tiltax = p{9};
end

end