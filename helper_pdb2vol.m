function [vol,data] = helper_pdb2vol(pdb,pix,trim,savemat)
if nargin<2, error('requires both pdb and pixel size inputs'), end
if nargin<3, trim=0; end %don't trim by default, but does auto-trim singles
if nargin<4, savemat=1; end %make name-val to generate .mat name with prefix/suffix?
%need user input or name/val to set what each model is for class/ID

%pdb to atoms
[~,~,ext] = fileparts(pdb);
if strcmp(ext,'.mat') %if .mat, load the data from the file
    try q = load(pdb); data = q.data;
    catch warning('Input is not a pdb2vol-generated .mat file'); end %#ok<SEPEX>
elseif strcmp(ext,'.pdb') %if .pdb, parse the file into a data variable
    data = internal_pdbparse(pdb);
end

%probably faster to store coords as normal array and atom id as a separate cell array or string array
%save each variable independently to the mat

vol = internal_volbuild(data,pix,trim);

if savemat==1
    %a1 = fullfile(path,append(base,'.mat')); %significantly slower for unknown reason
    a2 = strrep(pdb,'.pdb','.mat');
    save(a2,'data','-nocompression'); %slightly faster without compression
end

end


function [data] = internal_pdbparse(pdb)
fid = fileread(pdb); 
text = textscan(fid,'%s','delimiter','\n','CommentStyle',{'REMARK'}); %import each line individually
text = text{1}; %fix being inside a 1x1 cell array for no reason

%delete terminator and temperature/ANISOU records that mess with model reading and parsing
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
    %stringarray = string(model{i});
    chararray = char(model{i}); %convert to char array to make column operable
    chararray(:,[1:30,55:76]) = []; %delete columns outside coords and atom id
    %atomvec = chararray(:,26:27); atomvec = upper(strrep(string(atomvec),' ',''));
    atomvec = upper(strrep(string(chararray(:,25:26)),' ','')); %process atom ids to only letters
    
    data{i,1} = atomvec;
    %data{i,2} = coords;
    
    %can't add spaces in this vector
    %coordchar = [chararray(:,[1:8]),chararray(:,[9:17]),chararray(:,[18:24])];
    coords = [str2num(chararray(:,1:8)),str2num(chararray(:,9:17)),str2num(chararray(:,18:24))]'; %#ok<ST2NM>
    data{i,2} = coords;
    %coords(1:3,end);
    %coords2 = [str2double(chararray(:,1:8)),str2double(chararray(:,9:17)),str2double(chararray(:,18:24))]';
    %str2double can't op on vectors
    %x = str2num(chararray(:,1:8))';
    %y = str2num(chararray(:,9:16))';
    %z = str2num(chararray(:,17:24))';
    %stcoords = string(coordchar);
    %co2 = convertCharsToStrings(stcoords);
    %coord2 = sscanf(coordchar,'%f%f%f')
    %coord = textscan(coordchar,'%f%f%f',numel(atomvec))
    %coord{1}'
    %co = {coord{1}';coord{2}';coord{3}'}
end

%for i=1:models %loop through detected models
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
%end

end

function vol = internal_volbuild(data,pix,trim)

%initialize atomic magnitude information
mag = struct('H',0,'C',6+1.3,'N',7+1.1,'O',8+0.2,'P',15,'S',16+0.6);
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
lim = round( (adj+b)/pix +1); %array size to place into, same box for all models

models = numel(data(:,2)); emvol = cell(models,1); %pre-allocate stuff
if models==1, trim=1; end
for i=1:models
    atoms = data{i,1}; %single column, hopefully for speed
    %coords = data{i,2}; %column per atom, hopefully faster indexing
    coords = round((data{i,2}+adj)./pix); %vectorized computing rounded atom bins outside the loop
    em = zeros(lim'); %initialize empty array
    for j=1:numel(atoms)
        opacity = mag.(atoms{j}); 
        x=coords(1,j); y=coords(2,j); z=coords(3,j); 
        em(x,y,z) = em(x,y,z)+opacity;
        
        %ever so slightly slower to do inline reference to atom mag
        %x3=coords(1,j); y3=coords(2,j); z3=coords(3,j); em(x3,y3,z3) = em(x3,y3,z3)+mag.(atoms{j});
        %pt = coords(:,j); %get coord for the atom record
        %co = round((pt+adj)./pix); %vectorize coordinate math
        
        %still don't know how to properly split vector into singles
        %x1=co(1); y1=co(2); z1=co(3); em(x1,y1,z1) = em(x1,y1,z1)+opacity;
        
        %x2=idx(1,j); y2=idx(2,j); z2=idx(3,j); em(x2,y2,z2) = em(x2,y2,z2)+opacity; 
        %no speed increase with inline reference to atom
        
        %x3=coords2(1,j); y3=coords2(2,j); z3=coords2(3,j); em(x3,y3,z3) = em(x3,y3,z3)+opacity;
    end
    if trim==1 %for bundles, should probably do singles here or helper_input for efficiency
        em = em(:,any(em ~= 0,[1 3]),:); 
        em = em(any(em ~= 0,[2 3]),:,:); 
        em = em(:,:,any(em ~= 0,[1 2]));
    end
    emvol{i} = em;
end

emvol = reshape(emvol,1,numel(emvol)); %make list horizontal because specifying it initially doesn't work
vol = emvol;
end
