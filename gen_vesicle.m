function [memvol,count,ves,vescen] = gen_vesicle(vol,num,pix,tries)
%randomly generates and places spherical vesicles into a volume without overlapping contents
%
%inputs:
%vol - 3d volume to place vesicles. does not need to be empty.
%num - number of different vesicles to generate
%pix - pixelsize of generated vesicles
%tries - number of placement attempts for each vesicle. default 2
%
%outputs:
%memvol - volume with only membranes. does not contain any prior contents of vol
%count - counts of successes (s) and failures (f) in attempting to place vesicles
%ves - cell array of each generated vesicle
arguments
    vol
    num
    pix
    tries = 2
end
%clipping out of the Z also conviniently how tomos actually look, but is maybe too random

count.s = 0; count.f = 0;
memvol = vol*0;
vescen = [];
%need to also store maps for each individual vesicle in a cell array corresponding to center
for i=1:num
    radi = (rand*300+150)/pix; %randomly generate inner radius of vesicle (need better range)
    rado = radi+(20+randi(12))/pix; %get outer radius from inner, should be constant something (7-9nm-ish?)
    %reduced outer radius distance for pearson, skew makes it wider
    offset = round(rado+20); %centroid offset to prevent negative values
    
    w = (rado-radi)/1.5; %deviation of the membrane distribution
    sf = [(rado^2)/(radi^2),(radi^2)/(rado^2)]/2; %factor to correct for excess inner density
    
    %fill space between radii with tons of points
    %ptnum = round(radi*5*(pix^3)*pi^2); %need to actually calculate volume of shell
    shellvol = 4/8*pi*(rado^3-radi^3); %volume of shell in pixels
    ptnum = round( 0.2*shellvol*pix^3 )*2; %convert to angstroms, scale to some arbitrary working density
    frac = [ptnum,ptnum*sf(2),ptnum*sf(1)]; %get fractions of the total to distribute between inner and outer
    rti = round(ptnum*sf(2)); rto = ptnum-rti; %partition density between inner and outer radii
    
    %ptrad = rand(ptnum,1)*(rado-radi)+radi; %uniform - flat monolayer
    switch 1
        case 1 %mirrored pearson - relatively hard inner and outer edges
            ptrad = [pearsrnd(radi,w,0.7,3,rti,1);pearsrnd(rado,w,-0.7,3,rto,1)];
        case 2 %mirrored gamma - a bit narrower, more edge smoothing
            ptrad = radi+[betarnd(3.0,6,rti,1);betarnd(6,3.0,rto,1)]*(rado-radi)*3.5;
    end
    %pearson is very slow, calls beta to call gamma which takes most of the time
    %figure(); histogram(ptrad);
    
    ptaz = rand(ptnum,1)*pi*2; %random circular azimuth angles
    %ptel = rand(1,ptnum)*pi*2; %causes asymmetry, polar density accumulation
    ptel = asin(2*rand(ptnum,1)-1); %random elevation angles, corrected for polar density accumulation
    
    [x,y,z] = sph2cart(ptaz,ptel,ptrad); %convert spherical coords to cartesian coords
    
    %generate empty array and round points to positive coords
    tmp = zeros(offset*2,offset*2,offset*2);
    x = round(x+offset); y = round(y+offset); z = round(z+offset);
    lipid = 6; %need to find the typical density of lipid membrane
    for j=1:numel(x) %loop through and add points as density to the shell
        tmp(x(j),y(j),z(j)) = tmp(x(j),y(j),z(j)) + lipid;
    end
    tmp = ctsutil('trim',tmp);
    ves{i} = tmp; %#ok<AGROW> %store trimmed vesicle into output cell array
    
    for q=1:tries %try to place each vesicle N times, allows for duplicates
        loc = round( rand(1,3).*size(vol)-size(tmp)/2 ); %randomly generate test position
        [vol,err] = helper_arrayinsert(vol,tmp,loc,'nonoverlap');
        count.f = count.f + err;
        if err==0
            [memvol] = helper_arrayinsert(memvol,tmp,loc); %to avoid weirdness with carbon grid doubling
            count.s = count.s+1; 
            vescen(end+1,:) = loc+round(size(tmp)/2); %#ok<AGROW>
        end
    end
    
end

end