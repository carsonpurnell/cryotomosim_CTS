function list = util_loadfiles(filter,prompt,multi)
% utility to automatically use the superior uipickfiles if found, otherwise use matlab's file UI
%
%
%prompt = 'Select input structure files';
%filter = '*.pdb;*.pdb1;*.mrc;*.cif;*.mmcif;*.mat';
arguments
    filter = '*.mat' %default backstop filetype filter
    prompt = 'Select input files' %window display title
    multi = []; %empty vector to multiselect
end

if exist('uipickfiles','file')==2 %&& strcmp(list,'gui') % preferred method of using GUI to find target files
    filter = replace(filter,'*','\'); filter = replace(filter,';','$|'); %replace separators
    filter = append(filter,'$'); %append maybe important last symbol
    %list = uipickfiles('REFilter','\.mrc$|\.pdb$|\.mat$|\.pdb1$|\.cif$|\.mmcif$','Prompt',prompt); 
    list = uipickfiles('REFilter',filter,'Prompt',prompt,'NumFiles',multi); 
    if ~iscell(list) || numel(list)==0, error('No files selected, aborting.'); end
else%if strcmp(list,'gui')
    if isempty(multi), multi='on'; else multi='off'; end %#ok<SEPEX> %parse multiselect flag
    [list, path] = uigetfile({filter},prompt,'MultiSelect',multi);
    if numel(string(list))==1, list={list}; end %fix single selection not being in a cell
    if ~iscell(list) || numel(list)==0, error('No files selected, aborting.'); end
    for i=1:numel(list) %make the list full file paths rather than just names so it works off-path
        list{i} = fullfile(path,list{i}); 
    end
end