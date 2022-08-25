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
        %dl = logical(dest(dx,dy,dz)); sl = logical(source(sx,sy,sz)); %faster but too inclusive
        dbin = imbinarize(rescale(dest(dx,dy,dz))); sbin = imbinarize(rescale(source(sx,sy,sz)));
        
        olog = dbin+sbin; olog = max(olog(:)); %fastest method, sum areas and find max to test if there was overlap
        if olog>1 %if overlap, record and output original
            overlap = 1; 
        else %if no overlap, add the source to the destination
            dest(dx,dy,dz) = source(sx,sy,sz) + dest(dx,dy,dz);
            overlap = 0;
        end
    case 'overlaptest' %much faster, test only variation
        %dl = logical(dest(dx,dy,dz)); sl = logical(source(sx,sy,sz)); %faster but too inclusive
        dbin = imbinarize(rescale(dest(dx,dy,dz))); sbin = imbinarize(rescale(source(sx,sy,sz)));
        
        olog = dbin+sbin; olog = max(olog(:)); %fastest method, sum areas and find max to test if there was overlap
        if olog>1 %if overlap, record and output original
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