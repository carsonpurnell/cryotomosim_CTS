function particles = helper_pdb2dat(file,pix,trim,centering,savemat)
%particles = helper_pdb2vol(pdb,pix,trim,centering,savemat)
%generates an EM density map(s) from an atomic structure definition file
%
%pdb - atomic structure file. valid types: .pdb/a, .cif/.mmcif, or .mat generated by this function
%pix - pixel size of the output density map in angstroms
%trim - trim empty space around volume (0 none, 1 trim group as one, 2 trim each individually)
%savemat - save a .mat of atomic information, loads much faster than full structure files
%
%vol - cell array of EM density maps
%data - cell array of atomic IDs and coordinates, equivalent to the saved .mat
arguments
    file
    pix
    trim = 1 %0 none, 1 by sum, 2 each individually
    centering = 0;
    savemat = 1 %by default, save a .mat file if possible as a much faster alternative
end
%calculate the scattering potentials in a subfunction? sum them over 2-4 angstroms?
%probably too low-res for specific angstrom distributions of signal to matter
%DO need to get a reasonable value for each atom/voxel, rather than hamfisting by Z
%E ai * exp(bi * s^2) seems to work, but is frequency dependent? on the E freq, or otherwise?
%sum over a spectrum of S values somehow for total signal? values far lower than even z
%data not centered, only happens in the pts2vol subfunct

%make the new general-purpose loader - parse and generate names/flags, and output a struct of all data
%source filename, group name, flags, sumvol, individual pts/vols/names/perim/shape? 
%should directly load .mat as well, and just check the fields (also reproject the vols/sumvol
%use knnsearch to find potential bonds? nearest 5 neighbors and check 2-5 to avoid self distance

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
%[vol,sumvol,names] = internal_volbuild(data,pix,trim,centering);
names = '';
sumvol = 0;
vol = {0};
%data
%disp(filename)
id = strsplit(filename,{'__','.'});

particles.name = id{1};
particles.flags = 'TODO';
%particles.bonds = 'notimplemented';
%particles.pix = pix;
%particles.atomcoords = data(:,2)'; %hopefully vertical now
for i=1:size(data,1)
    com = mean(data{i,2},1); %need radius from geometric, not mass center
    tmpco = data{i,2}-com;
    tmpint = atomdict(data{i,1},'sc')';
    %tmpz = dict_atoms(data{i,1});
    %tmpatomint = dict_intensity(tmpz,'sc');
    %if tmpatomint~=tmpint, disp('nonidentity'); end
    particles.adat{i} = [tmpco,tmpint];
    
    tmpco = double(unique(tmpco,'rows'));
    %instead of alphashape on everything, prune to fraction of points and AS that - much faster
    %then add >10% of all points back in and unique check at the end?
    %is it faster to iteratively find borders from chunks of points then do it again at the end?
    alphat = alphaShape(double(tmpco),12); %surprisingly slow
    [~,p] = boundaryFacets(alphat);
    %p = boundaryiter(tmpco); %faster iterative border finding
    %size(p2), size(p)
    %tr = delaunay(tmpco(:,1),tmpco(:,2),tmpco(:,3)); %don't have freeboundary function?
    %[f,pfree] = freeBoundary(tr);
    
    %k = boundary(tmpco,0.1); p2 = tmpco(k,:); %boundary is slower, just alphashape/facets with more overhead
    
    n = size(tmpco,1);
    ix = randperm(n); ix = ix(1:round(n/100));
    pi = tmpco(ix,1:3);
    p = single([p;pi]); %need to add back 1-3% or so of points to prevent inside placements
    particles.perim{i} = unique(p,'rows');
    mn = data{i,3};
    if strcmp(mn,'NA')
        particles.modelname{i} = id{min(i,end)};
    else
        particles.modelname{i} = data{i,3};
    end
end
%data
%need to center on 0,0,0
%particles.atomid = data(:,1)';
%particles.modelnames = names';
%particles.radius = max(pdist2(particles.atomcoords,single([0,0,0])));
%particles.vol = vol;
%particles.sumvol = sumvol;
%{
for i=1:numel(particles.atomid)
    particles.atomint{i} = atomdict(particles.atomid{i},'sc');
    particles.adat{i} = [particles.atomcoords{i},particles.atomint{i}'];
end
%}

if savemat==1 %.mat saving and check if file already exists
    outsave = fullfile(path,append(filename,'.mat'));
    if isfile(outsave) %don't overwrite an existing file
        fprintf(' .mat exists, '); 
    else %save without compression if no .mat file exists
        fprintf(' saving .mat... '); save(outsave,'data','-nocompression');
    end
end

end

function piter = boundaryiter(pts)
%pts = unique(pts,'rows');
n = size(pts,1);
ix = randperm(n); ix = repmat(ix,[2,1]);
iters = 5;
l = round(n/iters);
alpha = 12;

b = cell(1,iters);
for i=1:iters
    j = ix(1+l*(i-1):i*l);
    tmp = unique(pts(j,:),'rows');
    shape = alphaShape(tmp,alpha*2);
    [~,b{i}] = boundaryFacets(shape);
end
pcat = cat(1,b{:}); %concatenate boundary points
pcat = unique(pcat,'rows');

sh = alphaShape(pcat,alpha);
[~,piter] = boundaryFacets(sh);
end


function [data] = internal_pdbparse(pdb)
fid = fileread(pdb); 
%text = textscan(fid,'%s','delimiter','\n'); %slightly faster to not parse remarks at all
text = textscan(fid,'%s','delimiter','\n');%,'CommentStyle',{'REMARK'}); %import lines, ignoring comments slow
text = text{1}; %fix being inside a 1x1 cell array

%{
ix = strncmp(text,'REMARK',6); text(ix) = []; %clear remark lines
ix = strncmp(text,'ANISOU',6); text(ix) = []; %delete temp records
ix = strncmp(text,'CONECT',6); text(ix) = []; %delete bond information for now
%}
ix = strncmp(text,'TER',3); text(ix) = []; %clear terminator lines
ix = contains(text,{'REMARK','ANISOU','CONECT'}); text(ix) = []; 
%remove remark, temperature, and connectivity records
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
    
    coords = [str2num(chararray(:,1:8)),str2num(chararray(:,9:16)),str2num(chararray(:,17:24))]; %#ok<ST2NM>
    
    %using str2num because str2double won't operate on 2d arrays, and can't add spaces while vectorized
    data{i,2} = single(coords); %store coordinates
    data{i,3} = 'NA';
end

end

function [data] = internal_cifparse(pdb)
fid = fileread(pdb); 
text = textscan(fid,'%s','delimiter','\n'); %read in each line of the text file as strings
text = text{1}; %fix being inside a 1x1 cell array

modnames = text(strncmp(text,'data_',5)); %retrieve lines storing model names
modnames = strrep(modnames,'data_',''); %remove the leading line identifier from names
%need to remove trailing numbers
%disp(modnames)
modnames = regexprep(modnames,'\_\d$',''); %remove trailing numbers so submodels don't get split
modnames = erase(modnames,'.'); %remove periods that would break use as fieldnames
%disp(modnames)
%other special characters that would break fieldnames?
%probably want to modularize this back into multiple functions for better access

%ix = strncmp(text,'HETATM',6); text(ix) = []; %clear hetatm lines to keep CNOPS atoms only

headstart = find(strncmp(text,'_atom_site.group_PDB',20)); %header id start
headend = find(strncmp(text,'_atom_site.pdbx_PDB_model_num',29)); %header id end
loopend = find(strncmp(text,'#',1)); %all loop ends

data = cell(numel(headstart),2);
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
    
    data{i,1} = atoms; data{i,2} = single(coord)'; data{i,3} = modnames{i};
end


end

function [vol,sumvol,names] = internal_volbuild(data,pix,trim,centering)
%initialize atomic magnitude information in arrays - faster than accessing a struct dictionary
%mag = struct('H',0,'C',6+1.3,'N',7+1.1,'O',8+0.2,'P',15,'S',16+0.6); %atomic number+fractional H counts
H=1;
edat = {'H',0;'C',6+1.3;'N',7+1.1;'O',8+0.2;'P',15;'S',16+0.6;...
    'MG',12;'ZN',30;'MN',25;'F',9;'CL',17;'CA',20};
%H = 25;
%edat = {'H',0;'C',108+1.28*H;'N',130+1.13*H;'O',97+0.08*H;'P',267;'S',100+0.41*H};
%need to reformulate these into vectors, too annoying to work with as cell and refactored anyway
elements = edat(:,1); atomdict = cell2mat(edat(:,2));
%[elements,atomdict] = atomdictfn; %computed scattering values, mismatch against carbon/water/membrane
%atomint = atomdict2(edat{1,2});

%{
%shang/sigworth numbers (HCNOSP): backwards C and N?
25, 108,  130, 97,  100 and 267   V��3
 0, 1.55, 1.7, 2.0, 1.8 and 1.8 radii
interaction parameters for voltage 100=.92 200=.73 300=.65 mrad/(V*A) (multiply to get real value?)
% hydrogens per atom of c=1.3, n=1.1, o=0.2, s=0.6 to bypass needing to add hydrogens manually
from creighton 1993 - c1.28, n1.13, o.08, p0, s.41
https://books.google.com/books?hl=en&lr=&id=hu8T_kI1LrkC&oi=fnd&pg=PR11&ots=jLzwwqkY-h&sig=9WhfbpfFd0mj61ErskgnxWSeGHY#v=onepage&q&f=false
%hydrogen intensity instead shifted onto average H per organic atom, because PDB inconsistently use H records
%currently using atomic number, should use a more accurate scattering factor measure
%}

names = data(:,3);
ix = find(contains(names,'origin')); %get the index, if any, of the name origin in the model
%check to clear out other dummy submodels?
%very much need to condense the accumulated goblin code
%make the if switch only set the origin, everything else outside for less clutter
if centering==1 %&& isempty(ix)
    trim=0; %don't trim if the modeled is centered, because that uncenters it
    if isempty(ix) %assert origin at 0,0,0 if no submodel providing one
        origin = [0,0,0];
    else %meaure the centroid of the origin submodel when provided
        origin = mean(data{ix(1),2},2);
        data(ix,:) = []; names(ix) = []; %remove dummy submodels and names
    end
    [a,b] = bounds(vertcat(data{:,2}),1); %bounds of all x/y/z in row order
%     span = max(origin-a,b-origin); %get spans measured from the origin
%     spanpix = ceil(span/pix)+1*0;
%     lim = spanpix*2+1; %get pixel box from span from origin
%     adj = spanpix*pix+pix*1-origin; %adjustment to apply to coordinates to place them into the pixel box
    %{
elseif centering==1 && 5==4 %&& ~isempty(ix) %&& 5==4
    if ~isempty(ix)
        data(ix,:) = []; names(ix) = []; %remove dummy submodels and names
    end
    trim=0; %don't trim if a centroid is imposed, need to revise input options
    origin = mean(data{ix,2},2);
    data(ix,:) = []; names(ix) = []; %remove the origin model for cleanliness
    [a,b] = bounds(horzcat(data{:,2}),2); %bounds of all x/y/z in row order
    span = max(origin-a,b-origin); %get spans measured from the origin
    spanpix = ceil(span/pix)+1*0;
    lim = spanpix*2+1; %get pixel box from span, always off to ensure origin perfect center
    adj = spanpix*pix+pix*1-origin;
    %}
else
    if ~isempty(ix)
        data(ix,:) = []; names(ix) = []; %remove dummy submodels and names
    end
    [a,b] = bounds(vertcat(data{:,2}),1); %bounds of all x/y/z in row order
    origin = (a+b)/2; %get the geometric center of the points
%     span = max(origin-a,b-origin); %get spans measured from the origin
%     spanpix = ceil(span/pix)+1*0;
%     lim = spanpix*2+1; %get pixel box from span, always off to ensure origin perfect center
%     adj = spanpix*pix+pix*1-origin;
end
span = max(origin-a,b-origin); %get spans measured from the origin
spanpix = ceil(span/pix)+1*0;
lim = spanpix*2+1; %get pixel box from span, always off to ensure origin perfect center
adj = spanpix*pix+pix*1-origin;

models = numel(data(:,2)); emvol = cell(models,1); %pre-allocate stuff
for i=1:models
    atomid = data{i,1}; %single column, hopefully for speed
    [~,c] = ismember(atomid,elements); %index each atom entry with the associated e- magnitude dictionary
    
    badentries = find(c<1); %find entries not in the element register
    c(badentries)=[]; data{i,2}(:,badentries) = []; %remove bad entries
    
    coords = round((data{i,2}+adj)./pix); %vectorized computing rounded atom bins outside the loop
    atomint = atomdict(c); %logical index the atom data relative to the atomic symbols
    em = zeros(lim); %initialize empty volume for the model
    
    for j=1:numel(atomint) %faster loop, use vectorized converted atomic info faster than struct reference
        x=coords(j,1); y=coords(j,2); z=coords(j,3); %fetch individual coordinates
        em(x,y,z) = em(x,y,z)+atomint(j);
    end
    emvol{i} = em;
end

if trim==1 %trim empty planes from the border of the model (for everything except .complex models)
    emvol = ctsutil('trim',emvol);
end
sumvol = 0;
if trim==2 %don't sumvol unmatched vols, leave as 0
    for i=1:numel(emvol)
        emvol{i} = ctsutil('trim',emvol{i});
    end
else
    sumvol = sum( cat(4,emvol{:}) ,4); %sum all volumes
end

vol = reshape(emvol,1,numel(emvol)); %make list horizontal because specifying it initially doesn't work
end

function [el,sc] = atomdictfn
el = {'H','C','N','O','P','S','F','Na','Mg','Cl','K','Ca','Mn','Fe'}; %element symbols to use for lookup
sc = [0.5288,2.5088,2.2135,1.9834,5.4876,5.1604,1.8012,4.7758,5.2078,4.8577,8.9834,9.9131,7.5062,7.1637];
%scattering potentials computed as sum of first 5 parameters of atom form factor, holding s=0
end

function atomint = atomdict2(atomid)
el = {'H','C','N','O','P','S','F','Na','MG','Cl','K','Ca','Mn','Fe'}; %element symbols to use for lookup
sc = [0.5288,2.5088,2.2135,1.9834,5.4876,5.1604,1.8012,4.7758,5.2078,4.8577,8.9834,9.9131,7.5062,7.1637,0];
z = [1,6,7,8,15,16,9,11,12,17,19,20,25,26,0]; %15==0 for bad entries
H = z(1);
hp = [0,1.3,1.1,0.2,0,0.6,0,0,0,0,0,0,0,0,0];
sc = sc+H*hp; %add hydrogen contributions
%scattering potentials computed as sum of first 5 parameters of atom form factor, holding s=0

[~,ix] = ismember(atomid,el);
%errs = find(ix<1); ix(errs) = 15;
ix(ix<1 | ix>15) = 15;
atomint = single(z(ix));
end

function atomint = atomdict(atomid,mode)
if nargin<2, mode='z'; end
el = {'H','C','N','O','P','S','F','Na','MG','Cl','K','Ca','Mn','Fe'}; %element symbols to use for lookup
%scattering potentials computed as sum of first 5 parameters of atom form factor, holding s=0
sc = [0.5288,2.5088,2.2135,1.9834,5.4876,5.1604,1.8012,4.7758,5.2078,4.8577,8.9834,9.9131,7.5062,7.1637,0];
z = [1,6,7,8,15,16,9,11,12,17,19,20,25,26,0]; %15==0 for bad entries
hp = [0,1.3,1.1,0.2,0,0.6,0,0,0,0,0,0,0,0,0]; %average hydrogens per atom
switch mode
    case 'z', atomint = z+z(1)*hp;
    case 'sc', atomint = sc+sc(1)*hp;
end

[~,ix] = ismember(atomid,el);
%errs = find(ix<1); ix(errs) = 15;
ix(ix<1 | ix>15) = 15;
atomint = single(atomint(ix));
end
function atomZ = dict_atoms(atomid)
symbol = {'H','C','N','O','P','S','F','Na','MG','Cl','K','Ca','Mn','Fe'};
z = [1,6,7,8,15,16,9,11,12,17,19,20,25,26,0]; %15==0 for bad entries
[~,ix] = ismember(atomid,symbol); ix(ix<1 | ix>15) = 15; %match symbols to atomic number (or other)
atomZ = single(z(ix));
end
function atomint = dict_intensity(atomZ,mode)
sc = [0.5288,2.5088,2.2135,1.9834,5.4876,5.1604,1.8012,4.7758,5.2078,4.8577,8.9834,9.9131,7.5062,7.1637,0];
z = [1,6,7,8,15,16,9,11,12,17,19,20,25,26,0]; %15==0 for bad entries
hp = [0,1.3,1.1,0.2,0,0.6,0,0,0,0,0,0,0,0,0]; %average hydrogens per atom
switch mode
    case 'z'
        atomint = z+z(1)*hp;
    case 'sc'
        atomint = sc+sc(1)*hp;
end
atomint = single(atomint(atomZ));
end