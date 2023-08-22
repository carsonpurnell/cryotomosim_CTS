function particle = helper_filmono(input,pix,prop)
% particle = helper_filmono(input,prop)
%
%
%
arguments
    input
    pix
    prop = 'gui'
end
%input = 'actin_mono_fil2.cif'; 
%ang = -166.15; step = 27.3; flex = 12; minL = 20;
%prop = [-166.15,27.3,12,20];
if strcmp(input,'gui')
    [input, path] = uigetfile({'*.pdb;*.pdb1;*.mrc;*.cif;*.mmcif'},'Select input file'); %select file
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
sum = helper_atoms2vol(pix,dat.adat,[0,0,0])*3;
particle.name = dat.name;
particle.filprop = prop;
particle.sum = sum;
particle.adat = dat.adat;
particle.modelname = dat.modelname;

% example test filaments?
% script with a single 3d filament curve test and a second module for filling some volume?

end