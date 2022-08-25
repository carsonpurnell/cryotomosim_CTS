function [output] = helper_input(list,pixelsize,sv)
%outputs a cell array of 3d volumes ready to input into other tomosim functions
%list is a cell array of input files(pdb or mrc) and workspace variables, which can be mixed. 
%list=='user' opens a broswer for selecting inputs(files only) or for each time it is used in the cell array
%accepts a list of mixed pdb, mrc, and variables. alternatively, list=='user' allows browsing for only files
%pixelsize is required if any files are input. [pixelsize resolution] projects pdb files as that resolution.
%proc is passed to helper_preproc method. 0 = none, 1 = crop, 2+ = masking/thresholding+cropping
%files are processed with masking and cropping empty planes by default. proc==0 only performs cropping.
%the list of outputs is sorted largest to smallest by default. sorting==0 to skip sorting.
%LIMITATION: you can't preset the number of manual browser inputs, you only get one attempt

%no variables? require using separate group/variable constructor to make .mat files?

%uipickfiles optional (much better) GUI for selecting files

%any other file formats that are important to have supported?

%can this be reformulated to use argument validation?

%need a desc/filename field for structs too so that things can be traced
%systemize names, break name at first underscore?

if nargin<2, error('first 2 inputs are required'), end
%by default, preproc everything
if nargin<3, sv=1; end 
%resolution projection is not used anymore, resizing also breaks intensity scaling
if numel(pixelsize)==1, pixelsize(2)=pixelsize(1); end


if strcmp(list,'gui') %preferred method of using GUI to find target files
    try %use uipickfiles if available, it's significantly better
        list = uipickfiles('REFilter','\.mrc$|\.pdb$|\.mat$'); %need to add .mat once generator is made
        if ~iscell(list) || numel(list)==0, error('No files selected, aborting.'); end
    catch %if not uipickfiles, use inferior matlab version
        [list, path] = uigetfile({'*.pdb;*.mrc'},'Select input files','MultiSelect','on');
        %similar checks to make sure files were selected, and if no bail immediately with error
        if numel(string(list))==1, list={list}; end
        if ~iscell(list) || numel(list)==0, error('No files selected, aborting.'); end
        for i=1:numel(list) %make the list full file paths rather than just names so it works off-path
            list{i} = fullfile(path,list{i}); 
        end
    end
end

%output = cell(1,numel(list));
%fprintf('Input list: '), disp(list) %too much command window spam
types = {'single','bundle','complex','cluster','group'};

for i=1:numel(list)
    fprintf('Loading input %i  ',i)
    [~,filename,ext] = fileparts(list{i}); %get file name and extension
    id = strsplit(filename,{'__','.'}); %extract class ID from filename (leading string up to . or __)
    type = id{end}; %last item in the ID, need to match to a list somehow
    if sum(strcmp(type,types))==0, type='single'; end %should fix empty types
    tmp.type=type;
    trim=0;
    if strcmp(type,'bundle') || strcmp(type,'single'), trim=1; end
    %need to do trimming based on type
    %typecheck = strcmp(type,types)
    %check2 = types(typecheck); check2=check2{1}; %string of the actual type
    %can use .bundle to get the type, along with .complex and .cluster and so on
    %currently the last item before extension
    
    id = strrep(id,'-','_'); %change dashes to underscore, field names can't have dashes
    for j=1:numel(id)
        id{j} = string(id{j});
        if ~isempty(sscanf(id{1},'%f')) %fix names to start with a letter
            id{1} = strcat('a_',id{1});
        end
    end
    %store filename and classification id of object
    tmp.file = {filename}; 
    tmp.id = id;
    
    if iscellstr(list(i)) && (strcmp(ext,'.pdb') || strcmp(ext,'.mat'))  
        fprintf('reading %s ',filename)
        tmp.vol = helper_pdb2vol(list{i},pixelsize(1),trim,sv); %read pdb and construct as volume at pixel size
        fprintf('generating at %g pixelsize ',pixelsize(1))
    elseif iscellstr(list(i)) && strcmp(ext,'.mrc')
        fprintf('loading %s  ',filename)
        [tmp, head] = ReadMRC(list{i});
        fprintf('resizing from %g to %g pixel size',head.pixA,pixelsize(1))
        tmp.vol = imresize3(tmp,head.pixA/pixelsize(1));
        %need to filter mrc to make density maps clean, pdb are already good to go
        
%     elseif ~iscellstr(list(i)) && numel(size(list{i}))==3
%         fprintf('using %ix%ix%i input array  ',size(list{i}))
%         tmp.vol = list{i}; tmp.file = {'var'};
    elseif iscellstr(list(i)) %#ok<*ISCLSTR>
        error('Error: item %i in the input list is a string, but not a .pdb, .mat, or .mrc',i)
    end
    
    %id specification stuff
    if numel(tmp.vol)==1 %for single models, use the id (first split string)
        tmp.id = (tmp.id(1));
    elseif numel(tmp.vol)==numel(id)-2 %for specifically named models, distribute the ID parts
        for j=1:numel(tmp.vol)
            tid{j} = tmp.id{j}; %#ok<AGROW>
        end
        tmp.id = tid;
    else
        postnum = {1:numel(tmp.vol)}; %because string doesn't work on cell arrays that are not variables
        tmp.id = append(tmp.id{1},'_',string(postnum{:}));
    end
    
    %tmp.vol = helper_preproc(tmp.vol,proc);
    
    output(i) = tmp; %#ok<AGROW> %store in multidim struct for ease of use
    
    fprintf('  done\n')
end

%sorting not relevant anymore with a struct per class. holding temporarily in case of within-class sorting
% if sorting~=0 %sort largest to smallest
%     counts = cellfun(@numel,output(1,:));
%     [~, ix] = sort(counts,'descend');
%     output = output(:,ix); %output = output(2,ix);
% end

end