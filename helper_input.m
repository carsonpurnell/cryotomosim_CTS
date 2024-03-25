function [particleset] = helper_input(list,pixelsize,sv)
%[particleset] = helper_input(list,pixelsize,sv)
%outputs a cell array of 3d volumes ready to input into other tomosim functions
%list is a cell array of input files(pdb or mrc) and workspace variables, which can be mixed. 
%list=='gui' opens a broswer for selecting inputs(files only) or for each time it is used in the cell array
%pixelsize is required if any files are input. 

%any other file formats that are important to have supported?
arguments
    list
    pixelsize double
    sv = 1 %save generated .mat intermediates by default
end

if isstruct(list) && isfield(list,'type') %if the input is a formatted particle list, record and end
    particleset = list; return
end
%list = internal_load(list); %internal call to either uipickfiles or uigetfiles
filter = '*.pdb;*.pdb1;*.mrc;*.cif;*.mmcif;*.mat';
list = util_loadfiles(filter);

%types = {'single','bundle','complex','cluster','group','assembly','memplex','membrane','inmem','outmem'};
%complex one flag to make it place everything separately, rather than needing a complex for each type
%location flag (any/default, membrane, inside/outside vesicles) control locmaps
%grouping/class flag (complex, assembly, random pick from group, sum of group) control split placement
%class? flag for complex, centering complex, randomizing group, or set to sum
modelext = {'.pdb','.pdb1','.cif','.mmcif','.mat'};
flaglist = ["membrane" "vesicle" "cytosol" "complex" "assembly" "cluster" "bundle"];
%loc flags: membrane embedded, vesicle inside, cytosol outside - otherwise anywhere
%subpart flags: complex to place each subpart separately, assembly the same but not all subparts
%grouping flags? randomly pick from group, either use that submodel ID or always use the group ID?
%clustering flags: cluster for clumps and bundle for linear bundles, both separate class methods
%retry flag that changes how many placement attempts are made to increase prevalence? or add to other flags?

%instead of type, use label. first entry? use for generating filenames so complex are not bloated?
%does assembly imply complex, or should it add complex to the flags?

for i=1:numel(list)
    fprintf('Loading input %i ',i)
    [~,filename,ext] = fileparts(list{i}); %get file name and extension
    
    id = strsplit(filename,{'__','.'}); %extract class IDs from filename, delimited by . or __
    %{
    %do flag checks first
    %flagcheck = find(contains(id,'membrane'));
    %flagcheck = find(contains(id,'vesicle'));
    %flagcheck = find(contains(id,'cytosol'));
    %else any
    %would be a big sprawling mess of fallthroughs, need something better
    %}
    %or do a group-level find-contains for location flags etc?
    
    %do a single find-contains against all valid flags, collect them, remove from ids
    %remove duplicates and keep them all as a string array in a flag field?
    
    flagix = matches(id,flaglist); %not using find due to being slower
    %flags = []; %will be empty if no flags detected, usually 1x0 empty
    tmp.flags = id(flagix); %id(flagix) = []; %remove flags from list when detected
    tmp.flags = unique(tmp.flags); %remove duplicate flag entries for cleanliness
    %kind of a mess in randomfill, checking partial flags. split them up or make a sorter funct?
    
    %figure out the relevant trimming/centering for vol loading
    %
    trim = 2; centering = 0;
    if any(ismember(tmp.flags,{'complex','assembly'}))
        trim = 1;
    end
    if any(ismember(tmp.flags,{'membrane'}))
        trim = 0; centering = 1;
    end
    %'sum' flag to use the sumvol to make it easy to use a complex as a single class?
    %if plex trim=1? ismember{x,{'complex','assembly'})
    %if mem trim =0 and center = 1?
    %}
    
    %convert to vols and get names
    
    %assign names from file/filename
    %if lumper, use first id string for all names
    %if splitting, run through and replace NA names with same ID name
    
    % old janky class parser and loader thing
    %
    tmp.type = id{end}; %type is the last item in the parsed name, if at all
    if contains(filename,'distract','IgnoreCase',true); tmp.type = 'distract'; end
    %{
    if ismember(tmp.type,types)==0, tmp.type='single'; end %default to single with no type ID in name
    trim=1;
    if ismember(tmp.type,{'memplex','membrane'}) %don't trim complexes and membrane
        %redundant for memplex/mem, pdb2vol currently sets trim=0 for centering==1
        trim=0;
    elseif ismember(tmp.type,{'complex','assembly'}) %could make the else/otherwise or leave blank if preset
        trim=1;
    elseif ismember(tmp.type,{'single','group','cluster','bundle','inmem','outmem'}) %fully trim each subvol
        trim=2;
    end %trim anything except complex/assem
    %rework this as a switch case? if is messy
    
    if ismember(tmp.type,{'memplex','membrane'}) %don't center membrane embed stuff
        centering = 1;
    else
        centering = 0;
    end
    %}
    
    if iscellstr(list(i)) && ismember(ext,modelext)
        fprintf('read: %s ',filename)
        [tmp.vol,tmp.sumvol,names] = helper_pdb2vol(list{i},pixelsize,trim,centering,sv); 
        %read pdb and construct as volume at pixel size
        %pdb names read in as 'NA', cif are in column cell array of strings
        %fprintf('generating at %g A ',pixelsize)
    elseif iscellstr(list(i)) && strcmp(ext,'.mrc')
        fprintf('loading %s  ',filename)
        [tmpmrcvol, head] = ReadMRC(list{i});
        %fprintf('resizing from %g to %g pixel size',head.pixA,pixelsize)
        tmpmrcvol = imresize3(tmpmrcvol,head.pixA/pixelsize);
        tmp.vol{1} = tmpmrcvol;
        names = {'NA'};
    elseif iscellstr(list(i)) %#ok<*ISCLSTR>
        error('Error: item %i (%s) in the input list is a string, but not a valid file type',i,list{i})
    end
    fprintf('generating at %g A ',pixelsize)
    
    %disp(names); disp(id);
    
    %new name parser, import from pdb2vol and replace when needed
    %need to intelligently get group/single
    %disp(names)
    %names = mat2cell(names,numel(names));
    %disp(names)
    %size(names)
    %class(names)
    %names = string(names); %does not appear to be necessary
    for j=1:numel(names)
        %disp(names{j}); disp(id{j});
        if strcmp(names{j},'NA') %replace empty names with something parsed from the filename
            names{j} = id{min(j,end)};
        end
        if ~isempty(sscanf(names{j},'%f')) %detect id that do not start with a letter
            names{j} = strcat('fix_',names{j}); %append a letter when necessary
        end
        names{j} = strrep(names{j},'-','_');
    end
    names = reshape(names,1,[]); %reshape to single row to unbreak concatenation in randomfill
    tmp.id = names;
    %disp(names); disp(id);
    %disp(names)
    
    %{
    %old id/name parser
    id = strrep(id,'-','_'); %change dashes to underscore, field names can't have dashes
    % parse names block, might go after loading files
    for j=1:numel(id) %loop through ID parts to make them functional for field names
        id{j} = string(id{j}); %convert to string for consistency with other functions
        if ~isempty(sscanf(id{j},'%f')) %detect id that do not start with a letter
            id{j} = strcat('fix_',id{j}); %append a letter when necessary
        end
    end
    %}
    
    tmp.file = {filename}; %tmp.id = id; %store filename and classification id of object
    
    %{
    %id specification from filename
    if numel(tmp.vol)==1 || numel(tmp.vol)==numel(id)-2
        tmp.id = tmp.id(1:numel(tmp.vol));
    else
        postnum = {1:numel(tmp.vol)}; %because string doesn't work on cell arrays that are not variables
        tmp.id = append(tmp.id{1},'_',string(postnum{:}));
    end
    %}
    
    %tmp.sumvol = sumvol;
    
    %tmp.vol = helper_preproc(tmp.vol,proc);
    %need to filter mrc to make density maps clean, pdb are already good to go
    particleset(i) = tmp; %#ok<AGROW> %store in multidim struct for ease of use
    fprintf('  done\n')
end

end

function list = internal_load(list)
prompt = 'Select input structure files';
filter = '*.pdb;*.pdb1;*.mrc;*.cif;*.mmcif;*.mat';
multi = [];
if strcmp(list,'gui') && exist('uipickfiles','file')==2 &&6==8%preferred method of using GUI to find target files
    filter = replace(filter,'*','\'); filter = replace(filter,';','$|'); %replace separators
    filter = append(filter,'$'); %append maybe important last symbol
    %list = uipickfiles('REFilter','\.mrc$|\.pdb$|\.mat$|\.pdb1$|\.cif$|\.mmcif$','Prompt',prompt); 
    list = uipickfiles('REFilter',filter,'Prompt',prompt,'NumFiles',multi); 
    if ~iscell(list) || numel(list)==0, error('No files selected, aborting.'); end
elseif strcmp(list,'gui')
    if isempty(multi), multi='on'; else multi='off'; end %parse multiselect flag
    [list, path] = uigetfile({filter},prompt,'MultiSelect',multi);
    if numel(string(list))==1, list={list}; end
    if ~iscell(list) || numel(list)==0, error('No files selected, aborting.'); end
    for i=1:numel(list) %make the list full file paths rather than just names so it works off-path
        list{i} = fullfile(path,list{i}); 
    end
end
end