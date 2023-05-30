%% atomistic carbon grid generator
pix = 12;
%boxsize = pix*[200,300,50];
pad = [pix,pix,0];
radius = 1e4; %hole radius
thick = 150; %carbon thickness
hcen = [boxsize(1)/2+randi(600)-300,radius+30+randi(200)]; %offsets for the hole center
filmsize = boxsize+pad*2; filmsize(3) = thick;


%% direct flat pt gen
%carbon density 2-2.3 g/cm?
density = 2.0/12*(1e8^-3)*6.022e23; %carbons per A^3, approx 0.1
atomfrac = 4;
ccount = round(density*prod(filmsize)/atomfrac); %
carbons = rand(ccount,3).*filmsize-pad;
carbons(:,3) = carbons(:,3)+(boxsize(3)-thick)/2;

h = sqrt( (carbons(:,1)-hcen(1)).^2 + (carbons(:,2)-hcen(2)).^2 );% +rand(ccount,1)*0;
carbons = carbons(h>radius,:); %filter out points inside the grid hole
carbons(:,4) = 2.5088*atomfrac/2; %scatter val for carbon - /2 because intensity too high weirdly

c1 = helper_atoms2vol(pix,carbons,boxsize,[0,0,0]);

%% comparison fill alphashape density
tic

psn = round(prod(filmsize)/5000); %
ps = rand(psn,3).*filmsize-pad;
h = sqrt( (ps(:,1)-hcen(1)).^2 + (ps(:,2)-hcen(2)).^2 );% +rand(ccount,1)*0;
ps = ps(h>radius,:);
ps(:,3) = ps(:,3)+(boxsize(3)-thick)/2;

%sh = alphaShape(carbons(:,1:3),pix); %way too many points
%[~,ps] = boundaryFacets(sh);
vec = randn(size(ps)); mag = rand(size(ps,1),1)*60;
vec = mag.*vec./vecnorm(vec,2,2);

sh = alphaShape(ps+vec,40);
density = 2.0/12*(1e8^-3)*6.022e23; %carbons per A^3, approx 0.1
atomfrac = 4;
vp = randtess(density/atomfrac,sh,'v'); %fewer pts with more density for speed during modelgen
vp(:,4) = 2.5088*atomfrac/2;
plot(sh)

toc
c2 = helper_atoms2vol(pix,vp,boxsize,[0,0,0]);
%dif = c1-c2;
figure(); sliceViewer(c2);
figure(); histogram(c1(c1>0)); hold on; histogram(c2(c2>0));
%figure(); histogram(dif(dif>0));
carbons = vp;

%x = 1:100; y = 1.3.^x; plot(x,y);
%plot3(carbons(:,1),carbons(:,2),carbons(:,3),'.'); axis equal
%[vol] = helper_atoms2vol(pix,carbons,boxsize); sliceViewer(vol);
%[vol,solv] = helper_atoms2vol(pix,carbons,boxsize);
%sliceViewer(vol+solv);