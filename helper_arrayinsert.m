function [output, overlap] = helper_arrayinsert(dest,source,coord,method)
%inserts source into dest at [coords]. by default, adds to destination values - method 'sum'.
%'replace' replaces dest values with source instead, min/max/mean should be obvious.
%'nonoverlap' adds if it does not overlap extant values, otherwise returns the original and overlap=1
%'overlaptest' is nonoverlap, but doesn't operate on dest, only checks and returns an overlap value
arguments
    dest double
    source double
    coord (1,:) double
    method (1,1) string {mustBeMember(method,{'sum','replace','nonoverlap','overlaptest','min','max','mean'})} = 'sum'
end

[d1, d2, d3]=size(dest);
[s1, s2, s3]=size(source);

%dest
dx=max(1,coord(1)):min(d1,coord(1)+s1-1);
dy=max(1,coord(2)):min(d2,coord(2)+s2-1);
dz=max(1,coord(3)):min(d3,coord(3)+s3-1);

%source
sx=max(-coord(1)+2,1):min(d1-coord(1)+1,s1);
sy=max(-coord(2)+2,1):min(d2-coord(2)+1,s2);
sz=max(-coord(3)+2,1):min(d3-coord(3)+1,s3);

switch method
    case 'sum'
        dest(dx,dy,dz) = source(sx,sy,sz) + dest(dx,dy,dz);
    case 'replace'
        dest(dx,dy,dz) = source(sx,sy,sz);
    case 'nonoverlap'
        dbin = imbinarize(rescale(dest(dx,dy,dz))); sbin = imbinarize(rescale(source(sx,sy,sz)));
        if max(max(max(dbin+sbin)))>1 %if overlap, record and output original
            overlap = 1; 
        else %if no overlap, add the source to the destination
            dest(dx,dy,dz) = source(sx,sy,sz) + dest(dx,dy,dz);
            overlap = 0;
        end
    case 'overlaptest' %much faster, test only variation
        dbin = imbinarize(rescale(dest(dx,dy,dz))); sbin = imbinarize(rescale(source(sx,sy,sz)));
        if max(max(max(dbin+sbin)))>1 %if overlap, record and output original
            overlap = 1; 
        else %if no overlap, record result
            overlap = 0;
        end
    case 'min'
        dest(dx,dy,dz) = min( source(sx,sy,sz),dest(dx,dy,dz) );
    case 'max'
        dest(dx,dy,dz) = max( source(sx,sy,sz),dest(dx,dy,dz) );
    case 'mean'
        dest(dx,dy,dz) = (source(sx,sy,sz)+dest(dx,dy,dz))./2;
end

output = dest;

end