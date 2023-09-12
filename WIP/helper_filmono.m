function particle = helper_filmono(input,pix,prop)
% particle = helper_filmono(input,pix,prop)
%
%
% need to add an option to save the struct as a .fil, use the same location and probably filename
arguments
    input
    pix
    prop = 'gui'
end

if strcmp(input,'gui')
    [input, path] = uigetfile({'*.pdb;*.pdb1;*.mat;*.cif;*.mmcif'},'Select input file'); %select file
    input = fullfile(path,input); %parse filename for functional off-path use
end
%if ~isfile(input), error('file does not appear to exist - use "gui" for off-path files'); end
if strcmp(prop,'gui') || numel(prop)~=4
    prompt = {'Helical angle','Monomer step size','Filament flexibility','Minimum length'};
    tmp = inputdlg(prompt,'filament properties'); %present GUI dialogue for manual input
    prop = zeros(1,4);
    for i=1:4
       prop(i) = str2double(tmp{i}); %convert properties from char cells
    end
end

%do the file loading
dat = helper_pdb2dat(input,2,0,1,0);
[sum,~,~,split] = helper_atoms2vol(pix,dat.adat,[0,0,0]); 
% 3 is for temporary fix for scatter val in Z number space
particle.name = dat.name;
particle.filprop = prop;
particle.sum = sum;
particle.vol = split;
particle.adat = dat.adat;
particle.modelname = dat.modelname;

% example test filaments?
% script with a single 3d filament curve test and a second module for filling some volume?
% option for generating a preview filament from input properties?
end