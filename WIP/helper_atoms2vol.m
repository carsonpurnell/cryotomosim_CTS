function [vol,solv,atlas,split] = helper_atoms2vol(pix,pts,sz,offset)
%[vol,solv,atlas,split] = helper_atoms2vol(pix,pts,sz,offset)
%projects a list of points as a 3d density volume
%4th dimension sets the weight value for each point, otherwise all weights are 1

%break into subfunctions for speed? not everything needs to run through multiple inputs
if isstruct(pts)
    names = fieldnames(pts); pts = struct2cell(pts);
else
    names = 0;
end
if iscell(pts) %this is a mess, subfunct/streamline
    s = numel(pts); t=1;
    if nargin<3
        offset = min(vertcat(pts{:}(:,1:3)),[],1)-pix;
        sz = max(vertcat(pts{:}(:,1:3)),[],1)+pix-offset;
    elseif nargin<4
        offset = [0,0,0];
    end
else
    s = 1; t=0;
    if nargin<3
        offset = min(pts(:,1:3),[],1)-pix;
        sz = max(pts(:,1:3),[],1)+pix-offset;
    elseif nargin<4
        offset = [0,0,0];
    end
end
%if nargin<4, offset=[0,0,0]; end
%if nargin<3, sz = max(catpts(:,1:3),[],1)+pix; end
%if size(pts,2)<4, pts(:,end+1)=1; end %intensity==1 if not given by 4th column
%need rough estimate of average volume for organic atoms
%very approximately 1.8a radii
%eventually might do individual vdw radii individually
avol = 4/3*pi*(1.65^3); %eyeballed volume of the average organic atom
h20 = 3.041; %computed scatter factor for H2O

emsz = floor(sz/pix); 
solv = (rand(emsz)-0.5)*0.2*pix^2+(pix^3);
%split = zeros([emsz,s]);
%split = zeros(emsz);
%sl = split;
%vl(1:s) = {split};
split = cell(1,s);
for j=1:s
    
    if t==1
        p = pts{j}; %split{j} = zeros(emsz);
    else
        p = pts;
    end
    
    if size(p,2)<4, p(:,4)=1; end %intensity==1 if not provided in 4th column
    m = p(:,4); p = p(:,1:3); p = round( (p-offset)/pix+0.5 );
    %p(:,1:3) = round((p(:,1:3)-offset)/pix+0.5); %very slow intermediate array assignments
    
    [split{j},solv] = internal_accum(p,m,avol,emsz,solv);
    
    %{
    for i=1:3
        ix = p(:,i) <= emsz(i) & p(:,i) >= 1; %get points inside the box
        p = p(ix,:); %drop points outside the box
    end
    for i=1:size(p,1)
        x=p(i,1); y=p(i,2); z=p(i,3); mag = m(i); %fetch data per atom
        %split(x,y,z,j) = split(x,y,z,j)+mag; %slow, 4d indexing very inefficient
        split(x,y,z) = split(x,y,z)+mag;
        solv(x,y,z) = solv(x,y,z)-avol*1;%(rand*.4+0.8);
    end
    %}
    %sliceViewer(solv-sl)
    %[a,b] = bounds(vl{1}-split,'all')
end
solv = max(solv,0)/35*h20; %compute waters in pixels from remaining volume
tmp = cat(4,zeros(emsz),split{:});
[~,atlas] = max(tmp,[],4); atlas = atlas-1;
vol = sum(tmp,4);
if iscell(names)
    %t = split; 
    clear split
    %t
    for i=1:s
        split.(names{i}) = tmp(:,:,:,i+1);
    end
end
end


function [vl,solv] = internal_accum(p,mag,avol,emsz,solv)
vl = zeros(emsz);

for i=1:3
    ix = p(:,i) <= emsz(i) & p(:,i) >= 1; %get points inside the box
    p = p(ix,:); %drop points outside the box
end
for i=1:size(p,1)
    x=p(i,1); y=p(i,2); z=p(i,3); m = mag(i); %fetch data per atom
    %split(x,y,z,j) = split(x,y,z,j)+mag; %slow, 4d indexing very inefficient
    vl(x,y,z) = vl(x,y,z)+m;
    solv(x,y,z) = solv(x,y,z)-avol*1;%(rand*.4+0.8);
end

end