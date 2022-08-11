function [proc, mask] = helper_preproc(input,method,verbose)
%generates a processed volume and mask of that volume from an input 3d volume
%method 2 masks statistically
%method 3 also masks with thresholding
%method=1 after all other methods, rescale to [0,1] --default
%method=0 and regardless of method crop empty planes from the output
%verbose=1 to see statistical measures over several iterations
arguments
    input (:,:,:) double
    method = 1 %don't know if this should be string-based or a numeric scale of extent
    verbose (1,1) = 0
end
[a, b] = bounds(input,'all');
%need to rework to either output back to the original scale of the inputs or work at the original scale of the
%input volume. scaling everything has intensity problems

vol = internal_crop(input);

vol = rescale(vol);
%how to organize layers of processing? 1-n that goes high to low work?
%still don't know how to do name-value pairs or optional string flags properly



if method==4
    vol = internal_densitymask(vol,1);
end

if method>1 %&& method<4
    for i=1:5
        %if method>3, vol = internal_densitymask(vol,i); end %heavy processing for density maps with background
        if method>2, vol = internal_thresh(vol,i); end %thresholding if method 3
        vol = internal_masking(vol,verbose); %statistical masking if method>1
    end
end

%rescale if method~=0
if method~=0
    vol = rescale(vol,0,b); %scaling from a creates super weirdness from negative values, can they be prevented?
end

proc = internal_crop(vol);
%generate the final mask - problem in not corresponding to the original input though
mask = single(proc>0);

end

function [crop] = internal_crop(in) %removes array planes that contain only zero to shrink volume
crop = in(:,any(in ~= 0,[1 3]),:); crop = crop(any(crop ~= 0,[2 3]),:,:); crop = crop(:,:,any(crop ~= 0,[1 2]));
end

function [out] = internal_masking(in,verbose)
out = rescale(in,0,1); out(out==0)=NaN;
avg = mean(out,'all','omitnan'); stdev = std(out,0,'all','omitnan'); variance = var(out,0,'all','omitnan');
if verbose==1, fprintf('Mean: %d   Stdev: %d   Variance: %d\n', avg, stdev, variance), end
mask = ((out-avg/4)>0).*((out-stdev/4)>0).*(out>(avg-stdev*3)).*(out>variance/2);
%out = fillmissing(out,'constant',0);
out = mask.*in;
end

function [out] = internal_thresh(in,i)
out = rescale(in);
bin = imbinarize(out,graythresh(out)*(i/(i+5)));
bin = imdilate(bin,strel('sphere',max(6-i,1)));
out = out.*bin;
end

function [out] = internal_densitymask(in,i)
out = rescale(in);
bin = imbinarize(out,graythresh(out)*1.1);
bin = bwareaopen(bin,round( mean(size(bin))*(20/i) ),6);
dil = imclose(bin,strel('sphere',10-i));
dil = imdilate(dil,strel('sphere',4));
out = in.*dil;
out(out==0)=NaN; 
out = rescale(out); 
out = fillmissing(out,'constant',0);
end