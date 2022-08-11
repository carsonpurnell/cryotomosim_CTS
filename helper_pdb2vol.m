function [vol,data] = helper_pdb2vol(pdb,pix,trim,savemat)
if nargin<2, error('requires both pdb and pixel size inputs'), end
if nargin<3, trim=0; end %don't trim by default, but does auto-trim singles
if nargin<4, savemat=1; end %make name-val to generate .mat name with prefix/suffix?
%need user input or name/val to set what each model is for class/ID

%pdb to atoms
[path,base,ext] = fileparts(pdb);
if strcmp(ext,'.mat') %if .mat, load the data from the file
    try q = load(pdb); data = q.data;
    catch warning('Input is not a pdb2vol-generated .mat file'); end
elseif strcmp(ext,'.pdb') %if .pdb, parse the file into a data variable
    data = internal_pdbparse(pdb);
end

vol = internal_volbuild(data,pix,trim);

if savemat==1
    %a1 = fullfile(path,append(base,'.mat'));
    a2 = strrep(pdb,'.pdb','.mat');
    %save(a1,'data'); % %30% of runtime for some reason?
    save(a2,'data'); %slightly faster
end

end


function data = internal_pdbparse(pdb)
%fid = fopen(pdb); 
fid = fileread(pdb); 
text = textscan(fid,'%s','delimiter','\n','CommentStyle',{'REMARK'}); %import each line individually
text = text{1}; %fix being inside a 1x1 cell array for no reason

%delete terminator and temperature/ANISOU records that mess with model reading and parsing
ix = strncmp(text,'TER',3); text(ix) = [];
ix = strncmp(text,'ANISOU',6); text(ix) = [];
ix = strncmp(text,'HETATM',6); text(ix) = []; %delete heteroatoms for sanity

modstart = find(strncmp(text,'MODEL ',6)); %find start of each model entry
modend = find(strncmp(text,'ENDMDL',6)); %find the end of each model entry

if isempty(modstart) %if single-model, extract start and end of atom lines
    modstart = find(strncmp(text(1:round(end/2)),'ATOM  ',6)); modstart = modstart(1);
    endatm = find(strncmp(text(modstart:end),'ATOM  ',6)); endatm = endatm(end);
    endhet = find(strncmp(text(modstart:end),'HETATMjj',6)); 
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
for i=1:models %loop through detected models
    atom = cell(1,numel(model{i}));
    coords = zeros(3,numel(model{i}));
    for j=1:numel(model{i}) %loop through each line of the model
        line = model{i}{j}; %extract single line for reading
        atom{j} = upper(strrep(line(77:78),' ','')); %read atom identifier string, prune spaces
        if strcmp(atom{j},''), error('Bad PDB, try resaving with proper software'), end
        %current method, same as old but a space-adding step for insurance slows it slightly
        fv = [line(31:38), ' ', line(39:46),' ', line(47:54)]; %faster than append to add spaces
        coords(:,j) = sscanf(fv,'%f',[3 1]); %single-line read net faster than individually
        
        %can be buggy due to overflow of one axis into the next, so must go with slower one
        %coords2(1:3,j) = sscanf(line(31:54),'%f',[3 1]); %reads coords, sscanf>str2double>>str2num
        %coords(1,j) = sscanf(line(31:38),'%f'); %x %slightly faster than cell array, thought it would be more
        %coords(2,j) = sscanf(line(39:46),'%f'); %y
        %coords(3,j) = sscanf(line(47:54),'%f'); %z
        %c1=line(31:38); c2=line(39:46); c3=line(47:54);
        %fx = strjoin({c1,c2,c3}); %slower than append for some reason
        %fz = append(c1,' ',c2,' ',c3); %slower than expanded append for some reason
        %fc = append(line(31:38),' ',line(39:46),' ',line(47:54)); %still way too slow
        %the following method works, but is significantly slower
        %qq = cellstr(line([31:38;39:46;47:54]));
        %t2 = sscanf(sprintf('%s ', qq{:}),'%f',[3 1]); %it works!
    end
    data{i,1} = atom; data{i,2} = coords;
end

end

function vol = internal_volbuild(data,pix,trim)

%initialize atomic magnitude information
mag = struct('H',0,'C',6+1.3,'N',7+1.1,'O',8+0.2,'P',15,'S',16+0.6);
%{
%shang/sigworth numbers (HCNOSP): backwards C and N?
25, 108,  130, 97,  100 and 267   V·Å3
 0, 1.55, 1.7, 2.0, 1.8 and 1.8 radii
interaction parameters for voltage 100=.92 200=.73 300=.65 mrad/(V*A) (multiply to get real value?)

%messy heteroatom version
%mag = struct('H',0,'C',6,'N',7,'O',8,'P',15,'S',16,...
%    'MG',12,'ZN',30,'MN',25,'F',9,'CL',17,'CA',20);

%hydrogen intensity instead shifted onto average H per organic atom, because PDB inconsistently use H records
%currently using atomic number, should use a more accurate scattering factor measure
%have seen h=25, c=130, n=108, o=97, s=100, p=267 for va^3 scattering potentials - need inorganic atoms
% hydrogens per atom of c=1.3, n=1.1, o=0.2, s=0.6 to bypass needing to add hydrogens manually
%}

%faster, vectorized adjustments and limits to coordinates and bounding box
[a,b] = bounds(horzcat(data{:,2}),2); %bounds of all x/y/z in row order
adj = max(a*-1,0)+pix; %coordinate adjustment to avoid indexing below 1
lim = round( (adj+b)/pix +1); %array size to place into, same box for all models

models = numel(data(:,2)); emvol = cell(models,1); %pre-allocate stuff
if models==1, trim=1; end
for i=1:models
    atomid = data{i,1}; %single column, hopefully for speed
    points = data{i,2}; %column per atom, hopefully faster indexing
    em = zeros(lim'); %initialize empty array
    for j=1:numel(atomid)
        atom = atomid{j}; coord = points(:,j); %get id and coord for the atom record
        co = round((coord+adj)./pix); %vectorize coordinate math
        x=co(1); y=co(2); z=co(3); %still can't split vector into singles
        em(x,y,z) = em(x,y,z)+mag.(atom); %add atom mag to the proper coordinate
    end
    if trim==1 %for bundles, should probably do singles here or helper_input for efficiency
        em = em(:,any(em ~= 0,[1 3]),:); 
        em = em(any(em ~= 0,[2 3]),:,:); 
        em = em(:,:,any(em ~= 0,[1 2]));
    end
    emvol{i} = em;
end

emvol = reshape(emvol,1,numel(emvol)); %make horizontal because specifying it initially doesn't work
vol = emvol;
end
