function out = gen_carbongrid(vol,pix,grid)
arguments
    vol (:,:,:) double
    pix (1,1) double
    grid = [15 2000] %offsets are by default fallbacks
end
c = 6; %atomic number of carbon atoms
r = round(grid(2)*5/pix); %nm diameter to radius in A
thick = round(grid(1)*10/pix); %nm thickness to A thickness
l = round(size(vol,3)/2-thick/2);

if numel(grid)<3, grid(3) = size(vol,1)/2+randi(20); end
if numel(grid)<4, grid(4) = r+20+randi(10)-pix*1; end

density = 2.5/(1e8)^3/12*6.022e23; %convert from 2-3.5g/cm^3 to atom/ang^3, ==~0.1 atoms per a^3
densepix = pix^2.8*density; %convert to average atoms per pixel

carbon = zeros(size(vol,1),size(vol,2),thick);
atoms = round(densepix*numel(carbon)*0.4); %40% random

pts = rand(3,atoms).*size(carbon)'; %pregenerate points
dist = sqrt( (pts(1,:)-grid(3)).^2 + (pts(2,:)-grid(4)).^2 ) +rand(1,atoms)*2;
pts = round(pts(:,dist>r)+0.5); %filter out points inside the grid hole

for i=1:size(pts,2)
    x = pts(1,i); y = pts(2,i);  z = pts(3,i);
    carbon(x,y,z) = carbon(x,y,z) + c;
end

carbon = carbon+round(densepix*c*0.6).*imbinarize(carbon); %60% flat background for speed
vol(:,:,l:l+thick-1) = carbon;
out = vol;
end