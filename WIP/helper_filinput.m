function particles = helper_filinput(pix,input)

arguments
    pix
    input = 'gui'
end

elseif strcmp(list,'gui')
    [list, path] = uigetfile({'*.pdb;*.pdb1;*.mrc;*.cif;*.mmcif'},'Select input files','MultiSelect','on');
    if numel(string(list))==1, list={list}; end
    if ~iscell(list) || numel(list)==0, error('No files selected, aborting.'); end
    for i=1:numel(list) %make the list full file paths rather than just names so it works off-path
        list{i} = fullfile(path,list{i}); 
    end
end


end