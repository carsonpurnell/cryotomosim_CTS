function [atlas,roinames,indvol] = helper_particleatlas(cts,opt)

arguments
    cts struct %might make this work other than needing the cts struct
    opt.individual = 0 %by default don't save individual splits, only global atlas image
    opt.dynamotable = 0 %switch to generate dynamo table
    %opt.path %hoping that this inherits the path from higher-level functions
end

%move other important particles to splitmodel if they exist
if isfield(cts.model,'beads'), cts.splitmodel.beads = cts.model.beads; end
if isfield(cts.model,'mem'), cts.splitmodel.mem = cts.model.mem; end
if isfield(cts.model,'grid'), cts.splitmodel.grid = cts.model.grid; end

roinames = fieldnames(cts.splitmodel); %retrieve component names
indvol = cell(1,numel(roinames));
atlas = zeros(size(cts.splitmodel.(roinames{1})));

for i=1:numel(roinames)
    %filename = append('ind',string(i),'_',roinames{i},'.mrc');
    tmp = cts.splitmodel.(roinames{i});
    %if opt.bin==1; indvol=imbinarize(rescale(indvol{i})); end
    
    bin = imbinarize(rescale(tmp));
    atlas = atlas+bin*i;
    
    filename = append('ind',string(i),'_',roinames{i});%,'.mrc');
    if opt.individual==1
        WriteMRC(bin,cts.pix,append(filename,'.mrc'))
    end
    if opt.dynamotable==1
        dynamotable(tmp,filename);
    end
    
end

ident = char(strjoin(roinames,'_'));
WriteMRC(atlas,cts.pix,append('atlas_',ident,'.mrc'))

end


function dynamotable(split,filename)

cl = imclose(split,strel('sphere',2)); %close small holes to avoid chopping most particles up
%er = imerode(cl,strel('sphere',1));

%binarizing doesn't matter with this scheme, so operating directly on image and keeping watershed label
d = -bwdist(~cl); %calculate distances for watershed
mask = imextendedmin(d,2);

%generate area mask, don't actually know why i'm using 2
d2 = imimposemin(d,mask); %local minima to make segmentations less chopped
w = watershed(d2); w(~cl) = 0; %label{i} = w; %do watershed on inverse distances
if 1==2 %diagnostic output of the watershed region map
    WriteMRC(w,2,append(filename,'_wshed_','.mrc'))
end
regions = regionprops3(~w==0,"Centroid"); %identify particle centroids

dtable = zeros(size(regions.Centroid,1),35); %pregenerate the full table with zeros
dtable(:,1) = 1:size(regions.Centroid,1); %particle number
dtable(:,24) = regions.Centroid(:,2); %x coords xyflip because matlab backwards points coordinates
dtable(:,25) = regions.Centroid(:,1); %y coords
dtable(:,26) = regions.Centroid(:,3); %z coords

%[path,out,ext] = fileparts(file{i}); %out = append(string(out),'.tbl');
fid = fopen(append(filename,'.tbl'),'wt'); %writing to file per column after transposing
fprintf(fid,[repmat('%g ', 1, size(dtable,2)) '\n'],transpose(dtable));
fclose(fid);

end