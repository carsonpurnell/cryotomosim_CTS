function [vol,solv,atlas,split] = helper_atoms2vol(pix,pts,sz,offset)
%[vol,solv,atlas,split] = helper_atoms2vol(pix,pts,sz,offset)
%
if isstruct(pts)
    names = fieldnames(pts); pts = struct2cell(pts);
else
    names = 0;
end
if iscell(pts)
    s = numel(pts); t=1;
    if nargin<3, sz = max(vertcat(pts{:}(:,1:3)),[],1)+pix; end
else
    s = 1; t=0;
    if nargin<3, sz = max(pts(:,1:3),[],1)+pix; end
end
if nargin<4, offset=[0,0,0]; end
%if nargin<3, sz = max(catpts(:,1:3),[],1)+pix; end
%if size(pts,2)<4, pts(:,end+1)=1; end %intensity==1 if not given by 4th column
%need rough estimate of average volume for organic atoms
%very approximately 1.8a radii
%eventually might do individual vdw radii individually
avol = 4/3*pi*(1.4^3); %eyeballed volume of the average organic atom
h20 = 3.041; %computed scatter factor for H2O

emsz = floor(sz/pix); 
solv = (rand(emsz)-0.5)*1*pix^2+(pix^3);
split = zeros([emsz,s]);
for j=1:s
    if t==1
        p = pts{j}; %split{j} = zeros(emsz);
    else
        p = pts;
    end
    if size(p,2)<4, p(:,4)=1; end %intensity==1 if not provided in 4th column
    p(:,1:3) = round((p(:,1:3)-offset)/pix+0.5);
    for i=1:3
        ix = p(:,i) < emsz(i) & p(:,i) > 1; %get points inside the box
        p = p(ix,:); %drop points outside the box
    end
    for i=1:size(p,1)
        x=p(i,1); y=p(i,2); z=p(i,3); mag = p(i,4); %fetch data per atom
        split(x,y,z,j) = split(x,y,z,j)+mag;
        solv(x,y,z) = solv(x,y,z)-avol*(rand*.4+0.8);
    end
end
solv = max(solv,0)/32*h20; %compute waters in pixels from remaining volume
tmp = cat(4,zeros(emsz),split);
[~,atlas] = max(tmp,[],4); atlas = atlas-1;
vol = sum(split,4);
if iscell(names)
    t = split;
    for i=1:s
        split.(names{i}) = t{i};
    end
end
end