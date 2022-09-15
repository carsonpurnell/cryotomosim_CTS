function [memvol,count,ves] = gen_vesicle(vol,num,pix)

arguments
    vol = zeros(300,300,100)
    num = 5
    pix = 13
end
%generate vesicles directly inside the volume?
%generate a struct of different vesicles probably too cumbersome and inflexible
%needs to be done prior to constraints, vesicles too large to stay inside border
%clipping out of the Z also conviniently how tomos actually look

count.s = 0; count.f = 0;
memvol = vol*0;
for i=1:num
    
    %generate random inner rad, compute outer rad from inner with pixelsize
    radi = (rand*300+100)/pix; %randomly generate inner radius of vesicle (need better range)
    rado = radi+50/pix; %get outer radius from inner, should be constant something
    offset = round(rado+20); %centroid offset to prevent negative values
    
    %fill space between radii with tons of points
    %generate a large number of random sph2cart radii,azimuth,elevation to convert into shell coordinates
    %does the shell data need to be pruned or cleaned in some way?
    %might need to add extra layer of points closer to inner and outer radii to get bilayer
    %ptnum = round(radi*5*(pix^3)*pi^2); %need to actually calculate volume of shell
    shellvol = 4/8*pi*(rado^3-radi^3); %in pixels
    ptnum = round( 0.5*shellvol*pix^3 ); %convert to angstroms, scale to some density
    ptrad = rand(1,ptnum)*(rado-radi)+radi;
    ptaz = rand(1,ptnum)*pi*2;
    %ptel = rand(1,ptnum)*pi*2; %cause asymmetry, polar density accumulation
    ptel = asin(2*rand(1,ptnum)-1);
    %convert spherical data to cartesian
    [x,y,z] = sph2cart(ptaz,ptel,ptrad);
    %[a,b] = bounds(x)
    %[a,b] = bounds(y)
    %[a,b] = bounds(z)
    %plot3(x,y,z,'.'); axis equal
    
    
    %generate empty array and round points to positive coords
    tmp = zeros(offset*2,offset*2,offset*2);
    x = round(x+offset);
    y = round(y+offset);
    z = round(z+offset);
    lipid = 3; %need to find the typical density of lipid membrane
    for j=1:numel(x) %loop through and add points as density to the shell
        tmp(x(j),y(j),z(j)) = tmp(x(j),y(j),z(j)) + lipid;
    end
    %trim and store the vesicle into the output array
    tmp = tmp(:,any(tmp ~= 0,[1 3]),:); 
    tmp = tmp(any(tmp ~= 0,[2 3]),:,:); 
    tmp = tmp(:,:,any(tmp ~= 0,[1 2]));
    ves{i} = tmp; %#ok<AGROW>
    
    
    %use internal function similar to randomfill internal testplace to try to place vesicle a few times
    %if placed, probably bail (maybe keep trying? avoids needing the check)
    for q=1:2
        %currently not doing any rotation with imwarp, no point with a sphere
        
        loc = round( rand(1,3).*size(vol)-size(tmp)/2 ); %randomly generate test position
        [vol,err] = helper_arrayinsert(vol,tmp,loc,'nonoverlap');
        count.f = count.f + err;
        if err==0
            [memvol] = helper_arrayinsert(memvol,tmp,loc); %to avoid weirdness with carbon grid doubling
            count.s = count.s+1; 
        end
        
    end
    
end

%disp(count)

%memvol=vol;
end