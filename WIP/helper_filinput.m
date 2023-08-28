function particle = helper_filinput(pix,list)

arguments
    pix
    list = 'gui'
end

if strcmp(list,'gui') %select the files themselves
    [list, path] = uigetfile({'*.pdb;*.pdb1;*.mrc;*.cif;*.mmcif;*.fil'},'Select input files','MultiSelect','on');
    if numel(string(list))==1, list={list}; end
    if ~iscell(list) || numel(list)==0, error('No files selected, aborting.'); end
    for i=1:numel(list) %make the list full file paths rather than just names so it works off-path
        list{i} = fullfile(path,list{i}); 
    end
end

for i=1:numel(list) %loop through and parse inputs through filmono for structs and directly for .fil
    [~,~,ext] = fileparts(list{i});
    switch ext
        case '.fil'
            qq = load(list{i},'-mat'); particle(i) = qq.monomer;  %#ok<AGROW>
            [sum,~,~,split] = helper_atoms2vol(pix,particle(i).adat,[0,0,0]); 
            %if ~iscell(split)
            convfac = 3; %approximate conversion factor from scatter value to Z numbers
            for j=1:numel(split)
                split{j} = split{j}*convfac;
            end
            particle(i).vol = split;  %#ok<AGROW>
            particle(i).sum = sum*convfac; %#ok<AGROW>
        case {'.cif','.mmcif','.pdb','.pdb1'}
            particle(i) = helper_filmono(list{i},pix); %#ok<AGROW>
    end
end

end