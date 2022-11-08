function [output, overlap] = helper_arrayinsert(dest,source,coord,method)
%inserts source into dest array at [coords] by the algorithm of "method".
%method can be sum(default),min,max,mean (self explanatory) replace, overlaptest, or nonoverlap
%replace replaces all values in the target region of dest
%overlaptest outputs overlap=1 if the sum would cause an overlap of objects, output is unmodified dest
%nonoverlap works as overlaptest, but if overlap!=1 also adds source to dest as output
arguments
    dest double
    source double
    coord (1,:) double
    method string {mustBeMember(method,{'sum','replace','nonoverlap','overlaptest','min','max','mean'})} = 'sum'
end

[d1, d2, d3]=size(dest);
[s1, s2, s3]=size(source);

%compute indexes for the target overlap region of the destination
dx=max(1,coord(1)):min(d1,coord(1)+s1-1);
dy=max(1,coord(2)):min(d2,coord(2)+s2-1);
dz=max(1,coord(3)):min(d3,coord(3)+s3-1);

%{
%linear index is slower than direct indexing, probably because of reshaping and indexing problems itself
s = [d1,d2,d3]; %size of array to index
index = dx.' + s(1)*(dy-1) + s(1)*s(2)*reshape(dz-1,1,1,numel(dz));
%if ~dest(index)==dest(dx,dy,dz), fprintf('xx'); end %test identity
%}
%{
dxi = max(1,  coord(1));
dxf = min(d1, coord(1)+s1-1);
dyi = max(1,  coord(2));
dyf = min(d2, coord(2)+s2-1);
dzi = max(1,  coord(3));
dzf = min(d3, coord(3)+s3-1);
%}

%compute indexes for the target overlap region of the source
sx=max(-coord(1)+2,1):min(d1-coord(1)+1,s1);
sy=max(-coord(2)+2,1):min(d2-coord(2)+1,s2);
sz=max(-coord(3)+2,1):min(d3-coord(3)+1,s3);
%{
try
sxi=sx(1); sxf=sx(end); %causing index errors if no overlap? inconsistently generated, rng(1) working
syi=sy(1); syf=sy(end);
szi=sz(1); szf=sz(end);
catch
    disp([sx,sy,sz]) %doesn't fix it
end
%}

%logical indexing is sometimes faster, but due to generating log/linear indexes is always slower
%using top/bottom values instead of ranges can cause errors, cause unknown (impossibly empty)

%is a loop faster by avoiding temporary array nonsense?
%avoid indexing on the right side somehow?

switch method
    case 'sum'
        %dest(index) = source(sx,sy,sz) + dest(index); %definitely slower
        %dest(dx,dy,dz) = source(sx,sy,sz) + dest(dx,dy,dz); %slightly slower
        tmp1 = source(sx,sy,sz) + dest(dx,dy,dz); dest(dx,dy,dz) = tmp1; %slightly faster somehow
    case 'nonoverlap' %first test if there would be overlap to save time
        %dl = logical(dest(dx,dy,dz)); sl = logical(source(sx,sy,sz)); %faster but too inclusive
        dbin = imbinarize(rescale(dest(dx,dy,dz))); sbin = imbinarize(rescale(source(sx,sy,sz)));
        
        overlap = dbin+sbin; overlap = max(overlap(:)); %fastest method to find potential overlaps?
        %ow = dbin+sbin-1; ow = any(ow(:)); %is any() faster than max()?
        if overlap>1 %if overlap, record and output original
            overlap = 1;
        else %if no overlap, add the source to the destination
            dest(dx,dy,dz) = source(sx,sy,sz) + dest(dx,dy,dz);
            overlap = 0;
        end
    case 'overlaptest' %faster than nonoverlap by only testing for overlap, will not do operations
        dbin = imbinarize(rescale(dest(dx,dy,dz))); sbin = imbinarize(rescale(source(sx,sy,sz)));
        
        overlap = dbin+sbin; overlap = max(overlap(:)); %fastest method to find potential overlaps
        if overlap>1 %if overlap, record and output original
            overlap = 1; 
        else %if no overlap, record result
            overlap = 0;
        end
    case 'replace'
        dest(dx,dy,dz) = source(sx,sy,sz);
    case 'min'
        dest(dx,dy,dz) = min( source(sx,sy,sz),dest(dx,dy,dz) );
    case 'max'
        dest(dx,dy,dz) = max( source(sx,sy,sz),dest(dx,dy,dz) );
    case 'mean'
        dest(dx,dy,dz) = (source(sx,sy,sz)+dest(dx,dy,dz))./2;
end

output = dest;
end