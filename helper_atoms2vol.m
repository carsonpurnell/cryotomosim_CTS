function [vol,solv] = helper_atoms2vol(pix,pts,sz,offset)
if nargin<4, offset=[0,0,0]; end
if nargin<3, sz = max(pts,[],1)+pix; end
%if size(pts,2)<4, pts(:,end+1)=1; end %intensity==1 if not given by 4th column
%need rough estimate of average volume for organic atoms
%very approximately 1.8a radii
%eventually might do individual vdw radii individually
avol = 4/3*pi*(1.4^3); %eyeballed volume of the average organic atom
h20 = 2.9; %overrounded number for water magnitude
pts(:,1:3) = round((pts(:,1:3)-offset)/pix+0.5);
emsz = floor(sz/pix); vol = zeros(emsz);
solv = (rand(emsz)-0.5)*1*pix^2+(pix^3);
for i=1:3
    ix = pts(:,i) < emsz(i) & pts(:,i) > 1; %get points inside the box
    pts = pts(ix,:); %drop points outside the box
end
for i=1:size(pts,1)
    x=pts(i,1); y=pts(i,2); z=pts(i,3); mag = pts(i,4); %fetch data per atom
    vol(x,y,z) = vol(x,y,z)+mag;
    solv(x,y,z) = solv(x,y,z)-avol*(rand*.4+0.8);
end
solv = max(solv,0)/32*h20; %compute waters in pixels from remaining volume
end