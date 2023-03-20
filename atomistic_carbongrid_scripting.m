%% atomistic carbon grid generator
pix = 12;
boxsize = pix*[300,500,50];
pad = [pix,pix,0];
radius = 1e4; %hole radius
thick = 150; %carbon thickness
hcen = [boxsize(1)/2+randi(120),radius+20+randi(120)]; %offsets for the hole center
filmsize = boxsize+pad*2; filmsize(3) = thick;

density = 2.5/(1e8)^3/12*6.022e23; %carbons per A^3, approx 0.1
ccount = round(density*prod(filmsize)/100);
carbons = rand(ccount,3).*filmsize-pad;
carbons(:,3) = carbons(:,3)+(boxsize(3)-thick)/2;

h = sqrt( (carbons(:,1)-hcen(1)).^2 + (carbons(:,2)-hcen(2)).^2 );% +rand(ccount,1)*0;
carbons = carbons(h>radius,:); %filter out points inside the grid hole

%x = 1:100; y = 1.3.^x; plot(x,y);
plot3(carbons(:,1),carbons(:,2),carbons(:,3),'.'); axis equal
