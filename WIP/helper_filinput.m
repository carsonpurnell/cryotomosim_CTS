function particle = helper_filinput(pix,list)
%particle = helper_filinput(pix,list)
%

arguments
    pix
    list = 'gui'
end

if strcmp(list,'gui') %select the files themselves
    filter = '*.pdb;*.pdb1;*.mrc;*.cif;*.mmcif;*.fil';
    list = util_loadfiles(filter); %auto parser for both uigetfile and uipickfiles
    %{
    [list, path] = uigetfile({filter},'Select input files','MultiSelect','on');
    if numel(string(list))==1, list={list}; end
    if ~iscell(list) || numel(list)==0, error('No files selected, aborting.'); end
    for i=1:numel(list) %make the list full file paths rather than just names so it works off-path
        list{i} = fullfile(path,list{i}); 
    end
    %}
end

for i=1:numel(list) %loop through and parse inputs through filmono for structs and directly for .fil
    [~,~,ext] = fileparts(list{i});
    switch ext
        case '.fil'
            qq = load(list{i},'-mat'); qq.monomer.perim = 0;
            particle(i) = qq.monomer;  %#ok<AGROW>
            [sum,~,~,split] = helper_atoms2vol(pix,particle(i).adat,[0,0,0]); 
            convfac = 3; %approximate conversion factor from scatter value to Z numbers
            for j=1:numel(split)
                split{j} = split{j}*convfac;
            end
            particle(i).vol = split;  %#ok<AGROW>
            particle(i).sum = sum*convfac; %#ok<AGROW>
        case {'.cif','.mmcif','.pdb','.pdb1','.mat'}
            particle(i) = helper_filmono(list{i},pix,'gui',1); %#ok<AGROW> %disimilar structures errors
            %not properly parsing this somehow, was there a change in data? perim borking it?
            %monomer = particle(i); 
    end
    if particle(i).perim == 0
        tmp = vertcat(particle(i).adat{:}); %cat all partitions for perim testing
        alphat = alphaShape(double(tmp(:,1:3)),12);
        [~,p] = boundaryFacets(alphat);
        particle(i).perim = p; %#ok<AGROW>
    end
    %}
end

end