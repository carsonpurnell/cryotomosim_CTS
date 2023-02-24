function [memvol,skel,nvecs,vesvol,count,ves] = gen_mem(vol,pix,num,tries,vecpts,memthick,memsize)
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
%vescen - list of vesicle centers of mass
%vesvol - cell array for volumes with each vesicle present individually (a giant memory sink that must change)
arguments
    vol
    pix
    num
    tries = 2
    vecpts = 9
    memthick = [60 24]
    memsize = 3
end
%clipping out of the Z also conviniently how tomos actually look, but is maybe too random

%pixel size <3 seems to infinite loop due to creating only empty blob vols - too much smoothing/dilating?
% future feature: nested membranes with shrinking/smoothing

%need more control options over thickness/radius and variability of both
%cell array of inputs for each? or just vector?

count.s = 0; count.f = 0;
memvol = vol*0;
%vescen = []; 
ves = cell(1,num);
vesvol = memvol; skel = vesvol;
nvecs = zeros(size(memvol,1),size(memvol,2),size(memvol,3),3);
label = 1;
for i=1:num
    switch 2%randi(2)
        case 1
            tmp = vesgen_sphere(pix); %generate spherical vesicles
            ves{i} = tmp; %store trimmed vesicle into output cell array
            tmpskel = vesskeletonize(tmp);
        case 2
            %lower pixel size can create empty blobs regularly
            tmpskel=0; %rs = 1;
            %thick = [28,12];%-rs;
            while ~any(tmpskel==1,'all')
                l = round(300/pix+20);
                sz = [l+randi(l*memsize),l+randi(l*memsize),l+randi(l*memsize)];
                [tmp,tmpskel] = vesgen_blob(sz,memthick,pix,6);
                %figure(); sliceViewer(tmp);
                %rs = rs+1;
            end
            %disp(rs)
    end
    
    for q=1:tries %try to place each vesicle N times, allows for duplicates
        loc = round( rand(1,3).*size(vol)-size(tmp)/2 ); %randomly generate test position
        [vol,err] = helper_arrayinsert(vol,tmp,loc,'nonoverlap');
        
        count.f = count.f + err;
        if err==0
            memvol = helper_arrayinsert(memvol,tmp,loc); %to avoid weirdness with carbon grid doubling
            skel = helper_arrayinsert(skel,tmpskel,loc); %write skeletons to the volume
            vesvol = helper_arrayinsert(vesvol,imbinarize(tmp)*label,loc); %label image of binary membranes
            
            %close all; sliceViewer(tmp); figure(); sliceViewer(tmpskel);
            
            norm4d = helper_volsurfnorm(tmpskel,vecpts);
            for j=1:3 
                nvecs(:,:,:,j) = helper_arrayinsert(nvecs(:,:,:,j),norm4d(:,:,:,j),loc); 
            end
            
            %vescen(label,:) = loc+round(size(tmp)/2); %#ok<AGROW>
            count.s = count.s+1; label = label+1;
            break
        end
    end
    
end
%sliceViewer(skel); %check skels in whole vol
%figure(); sliceViewer(memvol);
%skel = vesskeletonize(memvol); %generate skeleton of the final membrane volume
%possibly reimplement this per individual vesicle for speed, and to enable blobby ones

%also compute normals here based on the skelmap or a modified working version?
%put normals in dim 2-4? maybe use a sparse array or linear array for normal vecs?
[x,y,z] = ind2sub(size(skel),find(skel>0));
pts = [x,y,z];
%size(pts)
%v = vertexNormal(pts);
%surfnorm(x,y,z);
%disp(vescen)
%disp(count.s)
end


function skel = vesskeletonize(memvol)
bw = bwdist(~memvol); %calculate distances inside the shape
mask = rescale(imgradient3(bw))>0.5; %generate an inverse mask that approximates the border, minus the mid
skel = (bw.*~mask)>max(bw,[],'all')/2-1; %apply the mask to the distance map and threshold edge noise
skel = ctsutil('edgeblank',skel,2); %clear edge ples to depth 2, nothing will be placed there in any case
skel = bwareaopen(skel,20);
end

%need a lot more control options for spheres
function ves = vesgen_sphere(pix)

radi = (rand*700+150)/pix; %randomly generate inner radius of vesicle (need better range)
rado = radi+(14+randi(14))/pix; %get outer radius from inner, should be constant something (7-9nm-ish?)
%reduced outer radius distance for pearson, skew makes it wider
offset = round(rado+20); %centroid offset to prevent negative values
%still not sure how to do the radius and what the radial density curve should look like

w = (rado-radi)/1.5; %deviation of the membrane distribution
sf = [(rado^2)/(radi^2),(radi^2)/(rado^2)]/2; %factor to correct for excess inner density

%fill space between radii with tons of points
%ptnum = round(radi*5*(pix^3)*pi^2); %need to actually calculate volume of shell
shellvol = 4/8*pi*(rado^3-radi^3); %volume of shell in pixels
ptnum = round( 0.4*shellvol*pix^3 )*1; %convert to angstroms, scale to some arbitrary working density
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
%need to reformulate the math so that density is hard-bound between ri/ro in angstroms
%figure(); histogram(ptrad);

ptaz = rand(ptnum,1)*pi*2; %random circular azimuth angles
%ptel = rand(1,ptnum)*pi*2; %causes asymmetry, polar density accumulation
ptel = asin(2*rand(ptnum,1)-1); %random elevation angles, corrected for polar density accumulation

[x,y,z] = sph2cart(ptaz,ptel,ptrad); %convert spherical coords to cartesian coords

%generate empty array and round points to positive coords
tmp = zeros(offset*2,offset*2,offset*2);
x = round(x+offset); y = round(y+offset); z = round(z+offset);
lipid = 5.5; %need to find the typical density of lipid membrane
for j=1:numel(x) %loop through and add points as density to the shell
    tmp(x(j),y(j),z(j)) = tmp(x(j),y(j),z(j)) + lipid;
end
ves = ctsutil('trim',tmp);

end


function [blob,skel] = vesgen_blob(sz,thick,pix,beta)
ptvol = zeros(sz); 
if numel(thick)==1, thick(end+1)=10; end
%thickness as 1x2 vector of min,max membrane thickness?
thickness = thick(1)+rand*thick(2); r = (thickness/2)/pix; %actually computing the radius, not the diam

n = round(sqrt(sum(sz))); 
%might need to change to using points and generating an interpoland for better scaling
%keeping things on a fixed grid makes sizing difficult

%beta 3 is pretty good for quite round membrane bodies
%beta = 5; %1 = uniform, lower is U, higher is increasingly sharp peak
%can decrease for wilder shapes, increase for nearer sphere
minl = sz*0.22; maxl = sz*0.75; %min and max of points for use as blob origins
pts = minl+(maxl-minl).*betarnd(beta,beta,n,3); %generate beta dist points to allow center/siding
pts = round(pts);
for i=1:size(pts,1)
    x = pts(i,1); y = pts(i,2); z = pts(i,3);
    ptvol(x,y,z) = 1;
end
blobvol = bwdist(ptvol)<min(sz)/r/2;
for i=1:3
    smoothvol = imgaussfilt3(single(blobvol),r*2-i*1); %smooth out towards a rounder overall shape
    blobvol = smoothvol>0.2;
end
memvol = ctsutil('trim',smoothvol>0.1); %trim vol to save space
p = round(r*2);
memvol = padarray(memvol,[p p p]); %pad to prevent anything from clipping into edges
%trim and pad the vol to center and make sure nothing is touching the edge

skel = bwperim(memvol);

distmap = bwdist(skel); 
mask = distmap<r;
leafperim = bwperim(mask); %sliceViewer(leafperim+mask);

smmask = imgaussfilt3(single(mask*0.3+leafperim*1.5),2);
sm2 = smmask.*mask;

memnoise = rand(size(skel))*0.4.*sm2;
%sliceViewer(memnoise+sm2);
%need to compute density more precisely (variable component?) - also need control variables
dens = 0.6*pix^3;
blob = (memnoise+sm2)*dens;

end