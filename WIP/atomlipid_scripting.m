% atomistic lipid generation
sz = 300;
sp = 0.99;

[sh,pts,pts1,pts2] = potato(sz,sp);
plot(sh)

%{
plot3(pts(:,1),pts(:,2),pts(:,3),'.'); axis equal; hold on
plot3(pts1(:,1),pts1(:,2),pts1(:,3),'.'); axis equal
plot3(pts2(:,1),pts2(:,2),pts2(:,3),'.'); axis equal
%}

%{
%% generate potato by connecting points
%s = rand(20,3);
%size input and sphericity input
%size scales number of points and base radius, sphericity scales radius between fixed and variable 1/sph
%should give good spectrum of control with few needed parameters
%sz = 300; sp = 0.4; %sphericity scale - needs more impact, by magnitude is awkward
%interesting bugfeature: sp~.8 usually makes double membranes
%>~.85 is double-thick and not a good membrane model unfortunately, need to separate layers
%nesting bugfeature gone as the cost of (mostly) fixing the double layer/delamination bug
n = round(8+sz^(0.2+sp));
rad = sz*(0+sp); var = sz*(1-sp)*2; %probably change to 1/sp-1
iters = round(1+(1+1/sp)^0.5);
az = rand(n,1)*180; el = rand(n,1)*180; r = rand(n,1)*var+rad;
[x,y,z] = sph2cart(az,el,r);
R = makehgtform('xrotate',pi/2); R = R(1:3,1:3); %get rotation matrix (3x3 of full matrix)
pts = ([x,y,z])*R; %rotate about axis so points aren't clustered in Z (stays 0-centered though)
R = makehgtform('zrotate',rand*180); R = R(1:3,1:3);
pts = pts*R; %spin about Z randomly so blobs are isotropically disordered in-plane
% surround points with randn hulls?
qq = repmat(pts,round(10*(sp^0.5)),1); %replicate points by 10
d = randn(size(qq))*30/sp; %not unitized for variability?
%d = d./vecnorm(d,2,2); %vector directions, unitized

% vec = randn(size(spts));
% spd = rand(size(vec,1),1)*8+4;
% vec = vec./vecnorm(vec,2,2).*spd;
pts = qq+d;
pts = smiter(pts,1,9);
sh = alphaShape(pts); sh.Alpha = criticalAlpha(sh,'one-region')+sz;
%plot(sh); hold on
%plot3(pts(:,1),pts(:,2),pts(:,3),'.'); axis equal
%}
%% 
%{
%plot3(x,y,z,'.'); axis equal;
%scales badly, need an interpolator function to get good sizes without handholding
for i=1:iters
    sh = alphaShape(pts);
    sh.Alpha = criticalAlpha(sh,'one-region')*1+0*(1.5+i/3);
    pts2 = randtess(.01*i,sh,'s'); %no points?
    pts = [pts;pts2]*1; 
    v = randn(size(pts))*10; %this right here makes no sense
    pts = pts+v;
    pts = smiter(pts,2,19);
    [~,pts] = boundaryFacets(alphaShape(pts)); %getting duplicates here 
end
sh = alphaShape(pts); sh.Alpha = criticalAlpha(sh,'one-region')*(2);
%[~,pts] = boundaryFacets(sh);
[~,ptso] = boundaryFacets(alphaShape(pts,100));
%smoothing thing
%ptsa = smiter(ptso,1,9);
%{
[~,ix] = pdist2(ptso,ptso,'euclidean','Smallest',9);
ix = ix(2:end,:); ptsa = zeros(size(ptso));
for i=1:size(ix,2)
    mm = mean(ptso(ix(:,i),:));
    ptsa(i,:) = mm;
end
%}

ptsa = unique(ptso,'rows');
%the following should remove inner surface while closing envelope gaps
sh = alphaShape(ptsa); sh.Alpha = criticalAlpha(sh,'one-region')*(5);
plot(sh);
%prune sh to boundary points only to speed randtess?
%}
%{
%% potato, but interpolation
n = 100; sz = 100;
rad = sz/5; var = sz;
az = rand(n,1)*180; el = rand(n,1)*180; r = rand(n,1)*var+rad;
[x,y,z] = sph2cart(az,el,r);
R = makehgtform('xrotate',pi/2); R = R(1:3,1:3); %get rotation matrix (3x3 of full matrix)
pts = ([x,y,z])*R; %rotate about axis so points aren't clustered in Z (stays 0-centered though)
R = makehgtform('zrotate',rand*180); R = R(1:3,1:3);
pts = pts*R; %spin about Z randomly so blobs are isotropically disordered in-plane
x=pts(:,1);y=pts(:,2);z=pts(:,3);

t = [0;cumsum(sqrt(diff(x).^2+diff(y).^2+diff(z).^2))];
t = t/t(end);
ti = linspace(0,1,5000);
xx = spline(t,x,ti)';
yy = spline(t,y,ti)';
zz = spline(t,z,ti)';
sh = alphaShape(xx,yy,zz);sh.Alpha = criticalAlpha(sh,'one-region')*2; 

plot3(xx,yy,zz,'.'); axis equal; hold on; plot(sh)
%}
%{
%% spherical generation
pts = vesgen_sphere(1);
x=pts(:,1);y=pts(:,2);z=pts(:,3);
sh = alphaShape(pts,12); %extremely slow with so many points - and also out of memory!
plot(sh)
%plot3(x,y,z,'.'); axis equal
%}

%% generate from point expansion
%{
minl = sz*0.25; maxl = sz*0.75; %min and max of points for use as blob origins
beta = 40; n = 50;
%pts = -sz/2+(sz).*betarnd(beta,beta,n,3); %random points
pts = sz/2*randn(n,3);
plot3(pts(:,1),pts(:,2),pts(:,3),'.'); axis equal
%sh = alphaShape(pts,120);
%plot(sh)
%}

%% functionalized surface shape to a shell shape
thick = 30; %shape = sh;
[shell] = shape2shell(sh,thick);
%{
%% project potato as volume after shelling
thick = 30;
vpts = randtess(thick/1.2,sh,'s');
vec = randn(size(vpts)); vec = thick*vec./vecnorm(vec,2,2);
vpts = vpts+vec;
ai = ones(size(vpts,1),1);
bx = [200,200,200]*5;

shell = alphaShape(vpts,12); %slow, main bottleneck
%shell.Alpha = criticalAlpha(sh,'one-region')*0.0+12; %weirdly slow, secondary bottleneck
%}
plot(shell)

%% shell to point distributions
vpts = randtess(0.3,shell,'v');
%{
dens = .01;
[mi,ma] = bounds(vpts,1);
box = (ma-mi)+10;
vnum = round(dens*prod( box ));
vpts = rand(vnum,3).*box+mi-5;
%invol = inShape(shell,vpts); %now the main bottleneck - alphashape slowness
%vpts = vpts(invol,:);
%}
spts = randtess(10,shell,'s');
vec = randn(size(spts));
spd = rand(size(vec,1),1)*8+4;
vec = vec./vecnorm(vec,2,2).*spd;
spts=spts+vec;
%plot3(vpts(:,1),vpts(:,2),vpts(:,3),'.'); axis equal; hold on
%plot3(spts(:,1),spts(:,2),spts(:,3),'.'); axis equal


%% project surface and interior points as a volume map
fpts = [spts;vpts];
%ai = ones(size(vpts,1),1);
bx = [200,200,200]*4;
%vol = fnpt2vol(8,fpts,ones(size(fpts,1),1)',bx*2,-bx);
[vv,solv,atlas,split] = helper_atoms2vol(8,fpts,bx,-bx/2);
sliceViewer(vv);

%{
%% spherical vesicle test
[pts] = vesgen_sphere(600,30);
ma = max(pts,[],1); mi = min(pts,[],1);
ai = ones(size(pts,1),1)*5.5;
vol = fnpt2vol(10,pts,ai',ma-mi,mi);
sliceViewer(vol);

%% projection block
spts2=spts+randn(size(spts))*1;
apts = [vpts;spts2]; ai = ones(size(apts,1),1)*15;
ma = max(apts,[],1); mi = min(apts,[],1);
vol = fnpt2vol(10,apts,ai',ma-mi,mi);
sliceViewer(vol);
%}
%% internal functs

function [sh,pts,pts1,pts2] = potato(sz,sp)
n = round(8+sz^(0.2+sp));
rad = sz*(0+sp); var = sz*(1-sp)*2; %probably change to 1/sp-1
iters = round(1+(1+1/sp)^0.5);
az = rand(n,1)*180; el = rand(n,1)*180; r = rand(n,1)*var+rad;
[x,y,z] = sph2cart(az,el,r);
R = makehgtform('xrotate',pi/2); R = R(1:3,1:3); %get rotation matrix (3x3 of full matrix)
pts = ([x,y,z])*R; %rotate about axis so points aren't clustered in Z (stays 0-centered though)
R = makehgtform('zrotate',rand*180); R = R(1:3,1:3);
pts = pts*R; %spin about Z randomly so blobs are isotropically disordered in-plane
% surround points with randn hulls?
qq = repmat(pts,round(10*(sp^0.5)),1); %replicate points by 10
d = randn(size(qq))*30/sp; %not unitized for variability?
%d = d./vecnorm(d,2,2); %vector directions, unitized

% vec = randn(size(spts));
% spd = rand(size(vec,1),1)*8+4;
% vec = vec./vecnorm(vec,2,2).*spd;
pts = qq+d;
%pts = smiter(pts,1,9);
sh = alphaShape(pts); sh.Alpha = criticalAlpha(sh,'one-region')+sz;

tp = randtess(0.1,sh,'s');
sh = alphaShape(tp); sh.Alpha = criticalAlpha(sh,'one-region')+sz/8;
[~,pts] = boundaryFacets(sh);

pts1 = smiter(pts,1,9); %smiter not great, can average between faces. need dist cutoff at least.
pts2 = smiter(pts,1,30);
%pts = (pts+pts1)/2;

sh = alphaShape(pts1); sh.Alpha = criticalAlpha(sh,'one-region')+sz/5;

end

function ptsa = smiter(ptso,iter,nb)
ptsa = zeros(size(ptso));
for j=1:iter
    [~,ix] = pdist2(ptso,ptso,'euclidean','Smallest',nb);
    ix = ix(2:end,:);
    for i=1:size(ix,2)
        mm = mean(ptso(ix(:,i),:));
        ptsa(i,:) = mm;
    end
    ptso = ptsa;
end
end

function [shell] = shape2shell(shape,thick)
vpts = randtess(thick/1.2,shape,'s');
vec = randn(size(vpts)); vec = thick*vec./vecnorm(vec,2,2);
vpts = vpts+vec;
shell = alphaShape(vpts,24);
end

function [pts] = vesgen_sphere(r,thick)
radi = r;
rado=radi+thick;
%radi = (rand*700+150)/pix; %randomly generate inner radius of vesicle (need better range)
%rado = radi+(14+randi(14))/pix; %get outer radius from inner, should be constant something (7-9nm-ish?)
%reduced outer radius distance for pearson, skew makes it wider
%offset = round(rado+20); %centroid offset to prevent negative values
%still not sure how to do the radius and what the radial density curve should look like

w = (rado-radi)/1.5; %deviation of the membrane distribution
sf = [(rado^2)/(radi^2),(radi^2)/(rado^2)]/2; %factor to correct for excess inner density
%correction factor seems a bit off. little too much inner density still?

%fill space between radii with tons of points
%ptnum = round(radi*5*(pix^3)*pi^2); %need to actually calculate volume of shell
shellvol = pi*(rado^3-radi^3)*0.2; %volume of shell in pseudoatoms
%ptnum = round( 0.2*shellvol*1^3 )*2; %convert to angstroms, scale to some arbitrary working density
%frac = [ptnum,ptnum*sf(2),ptnum*sf(1)]; %get fractions of the total to distribute between inner and outer
rti = round(shellvol*sf(2)); %rto = ptnum-rti; %partition density between inner and outer radii
rto = round(shellvol*sf(1)); %slightly more points to balance out?
ptnum = rti+rto;
%ptrad = rand(ptnum,1)*(rado-radi)+radi; %uniform - flat monolayer
switch 1
    case 1 %mirrored pearson - relatively hard inner and outer edges
        ptrad = [pearsrnd(radi,w,0.7,3,rti,1);pearsrnd(rado,w,-0.7,3,rto,1)];
    case 2 %mirrored gamma - a bit narrower, more edge smoothing
        ptrad = radi+[betarnd(3.0,6,rti,1);betarnd(6,3.0,rto,1)]*(rado-radi)*3.0;
end
%pearson is very slow, calls beta to call gamma which takes most of the time
%need to reformulate the math so that density is hard-bound between ri/ro in angstroms
%figure(); histogram(ptrad);

ptaz = rand(ptnum,1)*pi*2; %random circular azimuth angles
%ptel = rand(1,ptnum)*pi*2; %causes asymmetry, polar density accumulation
ptel = asin(2*rand(ptnum,1)-1); %random elevation angles, corrected for polar density accumulation

[x,y,z] = sph2cart(ptaz,ptel,ptrad); %convert spherical coords to cartesian coords
pts = [x,y,z];
ves = 0;
%{
%generate empty array and round points to positive coords
tmp = zeros(offset*2,offset*2,offset*2);
x = round(x+offset); y = round(y+offset); z = round(z+offset);
lipid = 5.5; %need to find the typical density of lipid membrane
for j=1:numel(x) %loop through and add points as density to the shell
    tmp(x(j),y(j),z(j)) = tmp(x(j),y(j),z(j)) + lipid;
end
ves = ctsutil('trim',tmp);
%}

end

function vol = fnpt2vol(pix,pts,atomint,sz,offset)
if nargin<5, offset=[0,0,0]; end
volpts = round((pts-offset)/pix+0.5); %shift so things round well
emsz = floor(sz/pix); vol = zeros(emsz);
for i=1:3
    ix = find(volpts(:,i) > emsz(i) | volpts(:,i) < 1);
    volpts(ix,:) = [];
    atomint(:,ix) = [];
end
for i=1:numel(atomint)
    %x=volpts(1,i); y=volpts(2,i); z=volpts(3,i); %fetch individual coordinates
    x=volpts(i,1); y=volpts(i,2); z=volpts(i,3);
    %c = volpts(:,i); x=c(1); y=c(2); z=c(3); %simultaneous pull is slower for some reason
    try
    vol(x,y,z) = vol(x,y,z)+atomint(i);
    catch
        [x,y,z], atomint(i)  %#ok<NOPRT>
    end
end
end