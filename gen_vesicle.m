function [memvol,ves] = gen_vesicle(vol, num, pix)

arguments
    vol = zeros(300,300,100)
    num = 5
    pix = 5
end
%generate vesicles directly inside the volume?
%generate a struct of different vesicles probably too cumbersome and inflexible
%needs to be done prior to constraints, vesicles too large to stay inside border


for i=1:num
    
    %generate random inner rad, compute outer rad from inner with pixelsize
    radi = (rand*150+80)/pix; %randomly generate inner radius of vesicle (need better range)
    rado = radi+50/pix; %get outer radius from inner, should be constant something
    offset = round(rado+20); %centroid offset to prevent negative values
    
    %fill space between radii with tons of points
    %generate a large number of random sph2cart radii,azimuth,elevation to convert into shell coordinates
    %does the shell data need to be pruned or cleaned in some way?
    %might need to add extra layer of points closer to inner and outer radii to get bilayer
    ptnum = round(radi*5*(pix^3)*pi^2); %need to actually calculate volume of shell
    ptrad = rand(1,ptnum)*(rado-radi)+radi;
    ptaz = rand(1,ptnum)*pi*2;
    %ptel = rand(1,ptnum)*pi*2;
    ptel = asin(2*rand(1,ptnum)-1);
    %convert spherical data to cartesian
    [x,y,z] = sph2cart(ptaz,ptel,ptrad); %assymmetry: poles have higher density!
    %[a,b] = bounds(x)
    %[a,b] = bounds(y)
    %[a,b] = bounds(z)
    %plot3(x,y,z,'.'); axis equal
    
    
    %loop through each point and add as some scaled intensity into an array
    tmp = zeros(offset*2,offset*2,offset*2);
    x = round(x+offset);
    y = round(y+offset);
    z = round(z+offset);
    lipid = 7; %need to find the typical density of lipid membrane
    for j=1:numel(x)
        tmp(x(j),y(j),z(j)) = tmp(x(j),y(j),z(j)) + lipid;
    end
    %trim and store the vesicle into the output array
    tmp = tmp(:,any(tmp ~= 0,[1 3]),:); 
    tmp = tmp(any(tmp ~= 0,[2 3]),:,:); 
    tmp = tmp(:,:,any(tmp ~= 0,[1 2]));
    ves{i} = tmp; %#ok<AGROW>
    
    %use internal function similar to randomfill internal testplace to try to place vesicle a few times
    %if placed, probably bail (maybe keep trying? avoids needing the check)
    for q=1:3
        %currently not doing any rotation with imwarp, no point with a sphere
        
        
    end
    
    
end


%ves = 0;
memvol=vol;
end