function [carbon,perim] = gen_carbon(vol,pix,opt)
%[carbon,perim] = gen_carbon(vol,pix,opt)
%

arguments
    vol (1,3)
    pix = [] %empty or 0 returns atoms, real value returns the volume
    opt.thick = 150+randi(20)
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
[edge] = carbonshape(vol,opt);

density = 2.0/12*(1e8^-3)*6.022e23; %carbons per A^3, approx 0.1
atomfrac = 1; %pseudoatomic factor for speed
carbon = randtess(density/atomfrac*0.3,edge,'v'); %
cperim = randtess(.32,edge,'s');

perim = edge.Points; %perimeter pts of shape
n = size(carbon,1);
ix = randperm(n); ix = ix(1:round(n/10));
pi = carbon(ix,1:3);
perim = single([perim;pi;cperim]); perim = unique(perim,'rows');

carbon(:,4) = 2.5088*1.5*atomfrac;

if pix>0
    carbon = helper_atoms2vol(pix,carbon,vol,[0,0,0]);
    vol = round(vol/pix);
    for i=1:3
        while size(carbon,i)<vol(i)
            pv = zeros(1,3); pv(i)=1;
            carbon = padarray(carbon,pv,'post');
        end
    end
end


end

function [edge] = carbonshape(vol,opt)
pad = [50,50,0]; %padding to avoid edge effects

% any directional circular offset mostly working, needs better central point rng
theta = rand*2*pi;
t1 = sin(theta); t2 = cos(theta);
q = abs([t1*vol(1),t2*vol(2)]);
r = opt.radius-max(q)/2 -min(q)/4+rand*max(vol)/10;
hcen = [t1,t2]*r+vol(1:2)/2+(rand(1,2)-rand(1,2)).*vol(1:2)/20;

%hcen = [vol(1)/2+(rand-rand)*(1000),opt.radius+30+randi(240)]; %old offsets for the hole center
filmsize = vol+pad*2; filmsize(3) = opt.thick;

psn = round(prod(filmsize)/18000); % 18000 just looks nice and is fast, not evaluated or hypothesis-driven
ps = rand(psn,3).*filmsize-pad;
h = sqrt( (ps(:,1)-hcen(1)).^2 + (ps(:,2)-hcen(2)).^2 ); %find points inside hole
ps = ps(h>opt.radius,:); %remove points in hole
ps(:,3) = ps(:,3)+(vol(3)-opt.thick)/2; %adjust grid to the center of the volume

vec = randn(size(ps)); mag = rand(size(ps,1),1)*60;
vec = mag.*vec./vecnorm(vec,2,2);
%pts = ps+vec;

edge = alphaShape(ps+vec,40); %get the shape itself
end