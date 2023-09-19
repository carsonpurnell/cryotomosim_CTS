function list = util_loadfiles(filter,prompt,multi)
% utility to automatically use the superior uipickfiles if found, otherwise use matlab's file UI
% inputs
%filter: file extension filter string, in the format '*.mmcif;*.mat'  defaults to '*.mat'
%prompt: window display title as string, default 'Select input structure files'
%multi: exact number of input files, or [min,max]. empty vector [] allows any number, and is default
arguments
    filter = '*.mat' %default backstop filetype filter
    prompt = 'Select input files' %window display title
    multi = []; %empty vector to multiselect
end
if exist('uipickfiles','file')== 2 %&& strcmp(list,'gui') %better FEX browser UI
    filter = append(replace(replace(filter,'*','\'),';','$|'),'$'); %parse filter
    list = uipickfiles('REFilter',filter,'Prompt',prompt,'NumFiles',multi); 
    if ~iscell(list) || numel(list)==0, error('No files selected, aborting.'); end
else%if strcmp(list,'gui')
    if isempty(multi), multi='on'; else multi='off'; end %#ok<SEPEX> %parse multiselect flag
    [list, path] = uigetfile({filter},prompt,'MultiSelect',multi);
    if numel(string(list))==1, list={list}; end %fix single selection not being in a cell
    if ~iscell(list) || numel(list)==0, error('No files selected, aborting.'); end
    for i=1:numel(list) %append paths to allow reliable off-path use of files
        list{i} = fullfile(path,list{i}); 
    end
end