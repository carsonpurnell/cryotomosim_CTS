function [carbon] = gen_carbon(vol,pix,opt)

arguments
    vol (1,3)
    pix = [] %empty or 0 returns atoms, real value returns the volume
    opt.thick = 150
    opt.radius = 1e4
end


%% atomistic carbon grid generator
%pix = 12;
%boxsize = pix*[200,300,50];
%boxsize = vol;

if pix>0
    vol = vol*pix; %get atomistic scale when working from pixels
end

%{
pad = [20,20,0]; %padding to avoid edge effects
radius = opt.radius; %hole radius
thick = opt.thick; %carbon thickness

hcen = [vol(1)/2+randi(600)-300,opt.radius+30+randi(200)]; %offsets for the hole center
filmsize = vol+pad*2; filmsize(3) = opt.thick;
%}

%{
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
carbonsdeprec = carbons;
c1 = helper_atoms2vol(pix,carbons,boxsize,[0,0,0]);
%}

%% comparison fill alphashape density

%{
psn = round(prod(filmsize)/18000); % 18000 just looks nice and is fast, not evaluated or hypothesis-driven
ps = rand(psn,3).*filmsize-pad;
h = sqrt( (ps(:,1)-hcen(1)).^2 + (ps(:,2)-hcen(2)).^2 ); %find points inside hole
ps = ps(h>opt.radius,:); %remove points in hole
ps(:,3) = ps(:,3)+(vol(3)-opt.thick)/2; %adjust grid to the center of the volume

%sh = alphaShape(carbons(:,1:3),pix); %way too many points
%[~,ps] = boundaryFacets(sh);
vec = randn(size(ps)); mag = rand(size(ps,1),1)*60;
vec = mag.*vec./vecnorm(vec,2,2);

sh = alphaShape(ps+vec,40);
%}
sh = carbonshape(vol,opt);


d2 = 2.0/12*(1e8^-3)*6.022e23; %carbons per A^3, approx 0.1
atomfrac = 4;
vp = randtess(d2/atomfrac*3,sh,'v'); %fewer pts with more density for speed during modelgen
vp(:,4) = 2.5088*1*atomfrac;

carbon = vp;
if pix>0
    carbon = helper_atoms2vol(pix,vp,vol,[0,0,0]);
end



end

function edge = carbonshape(vol,opt)

pad = [20,20,0]; %padding to avoid edge effects
%radius = opt.radius; %hole radius
%thick = opt.thick; %carbon thickness

hcen = [vol(1)/2+randi(600)-300,opt.radius+30+randi(200)]; %offsets for the hole center
filmsize = vol+pad*2; filmsize(3) = opt.thick;

psn = round(prod(filmsize)/18000); % 18000 just looks nice and is fast, not evaluated or hypothesis-driven
ps = rand(psn,3).*filmsize-pad;
h = sqrt( (ps(:,1)-hcen(1)).^2 + (ps(:,2)-hcen(2)).^2 ); %find points inside hole
ps = ps(h>opt.radius,:); %remove points in hole
ps(:,3) = ps(:,3)+(vol(3)-opt.thick)/2; %adjust grid to the center of the volume

%sh = alphaShape(carbons(:,1:3),pix); %way too many points
%[~,ps] = boundaryFacets(sh);
vec = randn(size(ps)); mag = rand(size(ps,1),1)*60;
vec = mag.*vec./vecnorm(vec,2,2);

edge = alphaShape(ps+vec,40);

end