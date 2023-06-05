function [carbon] = gen_carbon(vol,pix,opt)

arguments
    vol (1,3)
    pix = [] %empty or 0 returns atoms, real value returns the volume
    opt.thick = 140+randi(20)
    opt.radius = 1e4+randi(400)-200
end

if pix>0
    vol = vol*pix; %get atomistic scale when working from pixels
end

%{
pad = [20,20,0]; %padding to avoid edge effects
radius = opt.radius; %hole radius
thick = opt.thick; %carbon thickness

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

sh = alphaShape(ps+vec,40);
%}
edge = carbonshape(vol,opt);

density = 2.0/12*(1e8^-3)*6.022e23; %carbons per A^3, approx 0.1
atomfrac = 1; %pseudoatomic factor for speed
vp = randtess(density/atomfrac*3,edge,'v'); %
vp(:,4) = 2.5088*1*atomfrac;

carbon = vp;
if pix>0
    carbon = helper_atoms2vol(pix,vp,vol,[0,0,0]);
end

end

function edge = carbonshape(vol,opt)
pad = [20,20,0]; %padding to avoid edge effects

%centering etc needs more control - at least transparency
hcen = [vol(1)/2+randi(600)-300,opt.radius+30+randi(200)]; %offsets for the hole center
filmsize = vol+pad*2; filmsize(3) = opt.thick;

psn = round(prod(filmsize)/18000); % 18000 just looks nice and is fast, not evaluated or hypothesis-driven
ps = rand(psn,3).*filmsize-pad;
h = sqrt( (ps(:,1)-hcen(1)).^2 + (ps(:,2)-hcen(2)).^2 ); %find points inside hole
ps = ps(h>opt.radius,:); %remove points in hole
ps(:,3) = ps(:,3)+(vol(3)-opt.thick)/2; %adjust grid to the center of the volume

vec = randn(size(ps)); mag = rand(size(ps,1),1)*60;
vec = mag.*vec./vecnorm(vec,2,2);

edge = alphaShape(ps+vec,40); %get the shape itself
end