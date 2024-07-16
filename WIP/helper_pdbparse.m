function data = helper_pdbparse(file)

[path,filename,ext] = fileparts(file);
switch ext %parse structure files depending on filetype
    case '.mat'
        try q = load(file); data = q.data;
        catch warning('Input is not a pdb2vol-generated .mat file'); end %#ok<SEPEX>
    case {'.cif','.mmcif'} %cif-parsed .mat files seem much larger than .pdb - what's happening?
        data = internal_cifparse(file);
    case {'.pdb','.pdb1'}
        data = internal_pdbparse(file);
end

end

function [data] = internal_pdbparse(pdb)
%fid = fileread(pdb); 
text = textscan(fileread(pdb),'%s','delimiter','\n'); %read in file text data (faster to include comment lines)
text = text{1}; %extract text from being in a 1x1 cell array

ix = strncmp(text,'TER',3); text(ix) = []; %clear terminator lines
ix = contains(text,{'REMARK','ANISOU','CONECT'}); text(ix) = [];  %remove remark, temperature, connect records
%ix = strncmp(text,'HETATM',6); text(ix) = []; %delete heteroatoms for sanity

modstart = find(strncmp(text,'MODEL ',6)); %find start of each model entry
modend = find(strncmp(text,'ENDMDL',6)); %find the end of each model entry

if isempty(modstart) %if single-model, extract start and end of atom lines
    %modelspan = strncmp(text(1:round(end/2)),'ATOM  ',6)+strncmp(text(1:round(end/2)),'HETATM',6);
    modelspan = find(strncmp(text,'ATOM  ',6)+strncmp(text,'HETATM',6)); %index all valid atom records
    modstart = modelspan(1); modend = modelspan(end);
    
    %endhet = find(strncmp(text(modstart:end),'HETATM',6)); %disabled by jj, hetatm currently ignored
    %if ~isempty(endhet), endhet = endhet(end); else endhet = 0; end %#ok<SEPEX>
    %modend = max(atomend,endhet)+modstart-1; %adjust for having searched only part of the list for speed
    
    model{1} = text(modstart:modend); models = 1;
elseif numel(modstart)==numel(modend) %if counts match, populate model counts
    models = numel(modstart); model = cell(models,1);
    for i=1:models %extract lines for each individual model
        model{i} = text(modstart(i)+1:modend(i)-1);
    end
elseif numel(modstart)~=numel(modend) %check if model numbers are the same
    error('failed to determine model numbers')
end

data = cell(numel(model),2);
for i=1:models
    chararray = char(model{i}); %convert to char array to make column operable
    chararray(:,[1:30,55:76]) = []; %delete columns outside coords and atom id, faster than making new array
    
    atomvec = upper(strtrim(string(chararray(:,25:26)))); %process atom ids to only letters
    %atomvec1 = chararray(:,25:26); atomvec1 = upper(strrep(string(atomvec1),' ','')); %slightly slower
    data{i,1} = atomvec; %store atoms
    
    coords = [str2num(chararray(:,1:8)),str2num(chararray(:,9:16)),str2num(chararray(:,17:24))]'; %#ok<ST2NM>
    %using str2num because str2double won't operate on 2d arrays, and can't add spaces while vectorized
    
    data{i,2} = single(coords); %store coordinates
    data{i,3} = 'NA';
end

end

function [data] = internal_cifparse(pdb)
text = textscan(fileread(pdb),'%s','delimiter','\n'); %read in each line of the text file as strings
text = text{1}; %extract text from being in a 1x1 cell array

modnames = text(strncmp(text,'data_',5)); %retrieve lines storing model names
modnames = strrep(modnames,'data_',''); %remove the leading line identifier from names
modnames = regexprep(modnames,'\_\d$',''); %remove trailing numbers so same-name submodels don't get split
modnames = erase(modnames,'.'); %remove periods that would break use as fieldnames

headstart = find(strncmp(text,'_atom_site.group_PDB',20)); %header id start
headend = find(strncmp(text,'_atom_site.pdbx_PDB_model_num',29)); %header id end
loopend = find(strncmp(text,'#',1)); %all loop ends

data = cell(numel(headstart),2); % pre-initialize storage cell array
for i=1:numel(headstart)
    loopend(loopend<headstart(i)) = []; %remove loop ends before current block
    header = text( headstart(i):headend(i) )'; %pull header lines
    header = replace(header,{'_atom_site.',' '},{'',''}); %clean bad chars from headers
    model = text( headend(i)+1:loopend(1)-1 ); %pull model lines from after header to loop end
    q = strrep(model,'" "','1'); %replace quoted spaces with 1 to fix blankspace errors
    
    q = textscan([q{:}],'%s','Delimiter',' ','MultipleDelimsAsOne',1); %read strings into cells
    q = reshape(q{1},numel(header),[])'; %reshape cells to row per atom
    t = cell2table(q,'VariableNames',header); %generate table from atoms using extracted headers
    
    atoms = t.type_symbol; %re-extract atom ID and coordinates from the temporary table
    x = char(t.Cartn_x); y = char(t.Cartn_y); z = char(t.Cartn_z);
    coord = [str2num(x),str2num(y),str2num(z)]';  %#ok<ST2NM> %merge coordinates into a single array
    
    data{i,1} = atoms; data{i,2} = single(coord); data{i,3} = modnames{i};
end

end