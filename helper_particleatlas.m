function [atlas,roinames,indvol] = helper_particleatlas(cts,individual,dynamotable)
%generates an atlas of the target (and mem/grid/beads) particles 
%
%cts is a cts struct from cts_model
%individual (default 0) if 1 also save individual binary images for particles (space-consuming)
%dynamotable (default 0) if 1 generates dynamo .tbl files for each type of particle

arguments
    cts struct %might make this work other than needing the cts struct
    individual = 0 %by default don't save individual splits, only global atlas image
    dynamotable = 0 %switch to generate dynamo table
end

%move other important particles to splitmodel if they exist
if isfield(cts.model,'beads'), cts.splitmodel.beads = cts.model.beads; end
if isfield(cts.model,'mem'), cts.splitmodel.mem = cts.model.mem; end
if isfield(cts.model,'grid'), cts.splitmodel.grid = cts.model.grid; end

roinames = fieldnames(cts.splitmodel); %retrieve component names
indvol = cell(1,numel(roinames));
atlas = zeros(size(cts.splitmodel.(roinames{1})));

for i=1:numel(roinames)
    indvol{i} = cts.splitmodel.(roinames{i}); %add model to stack
    
    bin = imbinarize(rescale(indvol{i})); %binarize model to add to atlas
    atlas = atlas+bin*i; %generate label image atlas based on model order for label intensity
    
    filename = append('ind',string(i),'_',roinames{i});
    if individual==1
        WriteMRC(bin,cts.pix,append(filename,'.mrc'))
    end
    if dynamotable==1
        generatetable(indvol{i},filename);
    end
    
end

ident = char(strjoin(roinames,'_'));
if length(ident)>60, ident=ident(1:60); end %truncation check to prevent invalidly long filenames
WriteMRC(atlas,cts.pix,append('atlas_',ident,'.mrc'))
end

function generatetable(split,filename)
cl = imclose(split,strel('sphere',2)); %close small holes to avoid chopping most particles up

d = -bwdist(~cl); %calculate distances for watershed
mask = imextendedmin(d,2); %generate extended minima mask from distance map (not sure why 2)
d2 = imimposemin(d,mask); %merge local minima to make segmentations less chopped
w = watershed(d2); w(~cl) = 0; %watershed the local minima and mask out non-particles
if 1==2 %diagnostic output of the watershed region map
    WriteMRC(w,2,append(filename,'_wshed_','.mrc'))
end
regions = regionprops3(~w==0,"Centroid"); %identify particle centroids from watershed mask

dtable = zeros(size(regions.Centroid,1),35); %pregenerate the full table with zeros
dtable(:,1) = 1:size(regions.Centroid,1); %particle number
dtable(:,24) = regions.Centroid(:,2); %x coords xyflip because matlab backwards points coordinates
dtable(:,25) = regions.Centroid(:,1); %y coords
dtable(:,26) = regions.Centroid(:,3); %z coords

fid = fopen(append(filename,'.tbl'),'wt'); %writing to file per column after transposing
fprintf(fid,[repmat('%g ', 1, size(dtable,2)) '\n'],transpose(dtable));
fclose(fid);
end