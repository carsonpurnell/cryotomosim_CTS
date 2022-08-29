function [vol,data] = helper_pdb2vol(pdb,pix,trim,savemat)
if nargin<2, error('requires both pdb and pixel size inputs'), end
if nargin<3, trim=0; end %don't trim by default, but does auto-trim singles
if nargin<4, savemat=1; end %make name-val to generate .mat name with prefix/suffix?

%pdb to atoms
[~,~,ext] = fileparts(pdb);
if strcmp(ext,'.mat') %if .mat, load the data from the file
    try q = load(pdb); data = q.data;
    catch warning('Input is not a pdb2vol-generated .mat file'); end %#ok<SEPEX>
elseif strcmp(ext,'.pdb') %if .pdb, parse the file into a data variable
    data = internal_pdbparse(pdb);
end

vol = internal_volbuild(data,pix,trim);

if savemat==1
    %a1 = fullfile(path,append(base,'.mat')); %significantly slower for unknown reason
    a2 = strrep(pdb,'.pdb','.mat');
    if isfile(a2)
        fprintf(' .mat exists, '); 
    else
        fprintf(' saving .mat... ')
        save(a2,'data','-nocompression'); %slightly faster without compression
    end
end

end


function [data] = internal_pdbparse(pdb)
fid = fileread(pdb); 
text = textscan(fid,'%s','delimiter','\n'); %slightly faster to not parse remarks at all
%text = textscan(fid,'%s','delimiter','\n','CommentStyle',{'REMARK'}); %import each line individually
text = text{1}; %fix being inside a 1x1 cell array for no reason

%delete terminator and temperature/ANISOU records that mess with model reading and parsing
ix = strncmp(text,'REMARK',6); text(ix) = []; %clear terminator lines
ix = strncmp(text,'TER',3); text(ix) = []; %clear terminator lines
ix = strncmp(text,'ANISOU',6); text(ix) = []; %delete temp records
ix = strncmp(text,'HETATM',6); text(ix) = []; %delete heteroatoms for sanity

modstart = find(strncmp(text,'MODEL ',6)); %find start of each model entry
modend = find(strncmp(text,'ENDMDL',6)); %find the end of each model entry

if isempty(modstart) %if single-model, extract start and end of atom lines
    modstart = find(strncmp(text(1:round(end/2)),'ATOM  ',6)); modstart = modstart(1);
    endatm = find(strncmp(text(modstart:end),'ATOM  ',6)); endatm = endatm(end);
    endhet = find(strncmp(text(modstart:end),'HETATMjj',6)); %disabled by jj, hetarm currently ignored
    if ~isempty(endhet), endhet = endhet(end); else endhet = 0; end %#ok<SEPEX>
    modend = max(endatm,endhet)+modstart-1; %adjust for having searched only part of the list for speed
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
    
    atomvec = upper(strrep(string(chararray(:,25:26)),' ','')); %process atom ids to only letters
    %atomvec1 = chararray(:,25:26); atomvec1 = upper(strrep(string(atomvec1),' ','')); %slightly slower
    data{i,1} = atomvec; %store atoms
    
    coords = [str2num(chararray(:,1:8)),str2num(chararray(:,9:17)),str2num(chararray(:,18:24))]'; %#ok<ST2NM>
    %using str2num because str2double won't operate on 2d arrays, and can't add spaces while vectorized
    data{i,2} = coords; %store coordinates
end

%defunct old code for looping through each atom record individually
%{
for i=1:models %loop through detected models
    %atomcell = cell(1,numel(model{i}));
    %atomcell = data{i,1};
    %atomst = strings(1,numel(model{i})); %string is unexpectedly slower than cell
    %coordchar = zeros(3,numel(model{i}));
    %vectorize pruning coords and atom ID to avoid doing so much in the loop?
    %veclabel{1:numel(model{i})} = model{i}{:}(77:78); %still doesn't work
    %veccoord = model{i}{:}(31:54); %mixmatch output errors from this
    %t = model{i};
    
    %atom = cellfun(@(x) textscan(x,'%*s%*s%*s%*s%*s%*f%f%f%f%*f%*f%c'),t, 'uni', false);
    %textscan through cell array works, but is >3x slower than the sscanf method inside the loop
    %is cellfun the slow part, or is textscan?
    
    %for j=1:numel(model{i}) %loop through each line of the model
        %line = model{i}{j}; %extract single line from pdb record
        %atomcell{j} = upper(strrep(line(77:78),' ','')); %read atom identifier while pruning spaces
        %this is the fastest known method, append/strjoin slower. needs spaces to avoid runover of long coords
        %fv = [line(31:38),' ',line(39:46),' ',line(47:54)]; %build coords string for atom
        %can sscanf bypass spaces by parsing the numbers directly from the line portion?
        %coords(:,j) = sscanf(fv,'%f',[3 1]); %read all coords for the atom record (>str2double>>>str2num)
        
        %{
        coord = model{i}{j}(31:54); 
        co = [coord(1:8),' ',coord(9:16),' ',coord(17:24)]; 
        coordchar(:,j) = sscanf(co,'%f',[3 1]);
        if strcmp(atomcell{j},''), error('Bad PDB, try resaving with proper software'), end
        %}
        
        %{
        %slightly slower version that parses line components first
        %label = model{i}{j}(77:78); atomcell2{j} = upper(strrep(label,' ','')); 
        coord = model{i}{j}(31:54); 
        co = [coord(1:8),' ',coord(9:16),' ',coord(17:24)]; 
        coords(:,j) = sscanf(co,'%f',[3 1]);
        %}
        
        %tt = textscan(line,'%*s%*s%*s%*s%*s%*f%f%f%f%*f%*f%c'); %slower than above, but only ~10%
        %coords2(1:3,j) = sscanf(line(31:54),'%f',[3 1]); %reads coords, sscanf>str2double>>str2num
        %coords(1,j) = sscanf(line(31:38),'%f'); %x %slightly faster than cell array, thought it would be more
        %coords(2,j) = sscanf(line(39:46),'%f'); %y
        %coords(3,j) = sscanf(line(47:54),'%f'); %z
    %end
    %data{i,1} = atomcell; 
    %data{i,2} = coordchar;
end
%}

end

function vol = internal_volbuild(data,pix,trim)

%initialize atomic magnitude information
%mag = struct('H',0,'C',6+1.3,'N',7+1.1,'O',8+0.2,'P',15,'S',16+0.6);
edat = {'H',0;'C',6+1.3;'N',7+1.1;'O',8+0.2;'P',15;'S',16+0.6};
elements = edat(:,1);
op = cell2mat(edat(:,2));
%elements = {'H','C',    'N',    'O',    'P', 'S'};
%op =       {'0','6+1.3','7+1.1','8+0.2','15','16+0.6'};

%{
%shang/sigworth numbers (HCNOSP): backwards C and N?
25, 108,  130, 97,  100 and 267   V·Å3
 0, 1.55, 1.7, 2.0, 1.8 and 1.8 radii
interaction parameters for voltage 100=.92 200=.73 300=.65 mrad/(V*A) (multiply to get real value?)
% hydrogens per atom of c=1.3, n=1.1, o=0.2, s=0.6 to bypass needing to add hydrogens manually

%messy heteroatom version
%mag = struct('H',0,'C',6,'N',7,'O',8,'P',15,'S',16,...
%    'MG',12,'ZN',30,'MN',25,'F',9,'CL',17,'CA',20);

%hydrogen intensity instead shifted onto average H per organic atom, because PDB inconsistently use H records
%currently using atomic number, should use a more accurate scattering factor measure
%}

%faster, vectorized adjustments and limits to coordinates and bounding box
[a,b] = bounds(horzcat(data{:,2}),2); %bounds of all x/y/z in row order
adj = max(a*-1,0)+pix; %coordinate adjustment to avoid indexing below 1
lim = round( (adj+b)/pix +1); %array size to place into, same initial box for all models

models = numel(data(:,2)); emvol = cell(models,1); %pre-allocate stuff
if models==1, trim=1; end
for i=1:models
    atomid = data{i,1}; %single column, hopefully for speed
    coords = round((data{i,2}+adj)./pix); %vectorized computing rounded atom bins outside the loop
    
    %convert atomic labels into atom opacity information outside the loop for speed
    [~,c] = ismember(atomid,elements); % get index for each atom indicating what reference it is
    atomint = op(c); %logical index the atom data relative to the atomic symbols
    
    em = zeros(lim'); %initialize empty volume for the model
    
    for j=1:numel(atomint) %faster loop, use vectorized converted atomic info faster than struct reference
        x=coords(1,j); y=coords(2,j); z=coords(3,j);
        em(x,y,z) = em(x,y,z)+atomint(j);
    end
    
    %{
    for j=1:numel(atomid) %slower loop, key-value struct is slow (but almost the fastest method)
        opacity = mag.(atomid{j}); %get atom mag from record - this is the slow step (faster than inline though)
        x=coords(1,j); y=coords(2,j); z=coords(3,j); %parse coords manually, no method to split vector
        %tmp = num2cell(coords(:,j)); [x1,y1,z1] = tmp{:}; %works but is much slower
        em(x,y,z) = em(x,y,z)+opacity; %write mag to the model vol
    end
    %}
    
    %{
    if em==em2
        %sliceViewer(em); figure(); sliceViewer(em2);
    else
        warning('not the same')
    end
    %}
    if trim==1 %trim empty planes from the border of the model (for everything except .complex models)
        em = em(:,any(em ~= 0,[1 3]),:); 
        em = em(any(em ~= 0,[2 3]),:,:); 
        em = em(:,:,any(em ~= 0,[1 2]));
    end
    emvol{i} = em;
end

vol = reshape(emvol,1,numel(emvol)); %make list horizontal because specifying it initially doesn't work
end