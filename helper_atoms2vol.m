function [vol,solv,atlas,split,acount] = helper_atoms2vol(pix,pts,sz,offset)
%[vol,solv,atlas,split,acount] = helper_atoms2vol(pix,pts,sz,offset)
%projects a list of points as a 3d density volume
%4th dimension sets the weight value for each point, otherwise all weights are 1

%scalar for thresholding ice vs content in the atlas to dilate/erode atlas labels

%round the values and combine outside the loop for speed?
%combining appears much slower, unique requires sorting

%additional arguments or reuse sz/offset for a 0-centered unbounded output?
%collapse sz/offset into a corner-to-counter bound, with a single 2x3 input
%if only 1x3 input, then first corner is assumped 0,0,0
%if input is only 0 or [] empty, center and trim?
%if input is [0,0,0] keep center at 0?

%corners = [offset;sz]; 
%if all(size(corners)==[2,3]) use for start and end corners of box
%elseif all(corners==[0,0,0]) keep 0-centered
%elseif all(size(corners)==[1,3]) use for start corner of box, find end corner automatically
%else (bound not given or given empty) auto box, tight on all sides

%break into subfunctions for speed? not everything needs to run through multiple inputs
if isstruct(pts)
    names = fieldnames(pts); pts = struct2cell(pts); %convert to cell for easy looping
else
    names = 0;
end
if ~iscell(pts) && numel(size(pts))==2
    pts = {pts}; %convert single array to cell for ease of use
end

%t = 1; 
s = numel(pts); %temp patch for old loop definition code and cell/array switcher

if nargin<3 %no box inputs, output tight bounds
    dd = vertcat(pts{:});
    offset = min(dd(:,1:3),[],1)-pix;
    sz = max(dd(:,1:3),[],1)+pix-offset;
elseif all(sz==[0,0,0]) %if box 0, keep the object centered in the volume
    dd = vertcat(pts{:});
    offset = -max(abs(dd(:,1:3)))-pix/2;
    sz = offset*-2;
elseif nargin<4 %box limit only, output corner starting at 0
    offset = [0,0,0];
end
%if size(pts,2)<4, pts(:,end+1)=1; end %intensity==1 if not given by 4th column
% rough constants - need improved values, per-atom vol especially
avol = 4/3*pi*(1.9^3); %eyeballed volume of the average organic atom (radii approx 1.8A)- get per-atom measure?
h20 = 3.041/2; %computed scatter factor for H2O - /2 for similarity to vol and simulate defaults
wd = 6.022e23/18/(1e8)^3; %molecules of water per a^3 - ~1/30 for liquid water
wvol = 32; %eyeballed volume of amorphous ice molecules in angstroms

emsz = floor(sz/pix); 

%solv = (rand(emsz)-0.6)*1.5*pix^2+(pix^3); %set initial solvent density
%solv = imgaussfilt3(solv,0.5); % smoother solvation test
solv = 0;
acount = zeros(emsz,'single');
sptmp = cell(1,s);
for j=1:s
    
    p = single(pts{j});
    %{
    if t==1
        p = pts{j}; %split{j} = zeros(emsz);
    else
        p = pts;
    end
    %}
    
    if size(p,2)<4, p(:,4)=1; end %intensity==1 if not provided in 4th column
    mag = p(:,4); p = p(:,1:3); p = round( (p-offset)/pix+0.5 );
    %p(:,1:3) = round((p(:,1:3)-offset)/pix+0.5); %very slow intermediate array assignments
    
    [sptmp{j},solv,acount] = test_accumarray(p,mag,emsz,solv,avol,acount);
    %[sptmp{j},solv2] = internal_accum(p,mag,avol,emsz,solv);
    %all(tmpvol==sptmp{j},'all')
    
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
%need to compute spread via smoothing
ex = pix/3;

%solv = (rand(emsz)-0.6)*1.5*pix^2+(pix^3)*1; %set initial solvent density
solv = (rand(emsz,'single')*0.5+ones(emsz,'single')*0.75)*pix^3; % initial solvent density from mean 1
%sliceViewer(solv)
solv = imgaussfilt3(solv,0.5);

%solv = imgaussfilt3(solv,ex);
acount = imgaussfilt3(acount,ex);
solv = max(solv-acount,0)/wvol*h20; %compute waters in pixels from remaining volume
%{
size(sptmp{1})
size(sptmp)
for i=1:numel(sptmp)
    figure(); sliceViewer(sptmp{i});
end
%}
tmp = cat(4,zeros(emsz),sptmp{:});
%size(tmp)
[~,atlas] = max(tmp,[],4); atlas = single(atlas-1);
vol = sum(tmp,4);
if iscell(names)
    for i=1:s
        split.(names{i}) = tmp(:,:,:,i+1);
    end
elseif s>0
    split = cell(1,s);
    for i=1:s
        split{i} = tmp(:,:,:,i+1);
    end
else
    split = 0;
end
end


function [tmpvol,tmpsolv,acount] = test_accumarray(p,mag,emsz,solv,avol,acount)
if nargin<6, acount = zeros(emsz,'single'); end
ixf = ones(size(mag),'single');
for i=1:3
    ix = p(:,i) <= emsz(i) & p(:,i) >= 1; %index points inside the box
    %p = p(ix,:); mag=mag(ix); %drop points outside the box
    %ixf = ixf.*ix;
    ixf = ixf & ix;
    %fprintf('%i dropped \n',numel(ix)-sum(ix)) %diagnostic to check if any points eliminated
end
ixf = ixf>0;
p = p(ixf,:); 
mag = single(mag);
mag=mag(ixf); %out of memory error - problematic reuse of same array? is mag too high precision?
tmpvol = accumarray(p,mag,emsz);
acount = acount+accumarray(p,1,emsz)*avol;
tmpsolv = solv;%-acount;
end

function [vl,solv] = internal_accum(p,mag,avol,emsz,solv)
vl = zeros(emsz);
%[u,ua,ub] = unique(p,'rows'); %slower than the rest due to needing to sort
%ixf = ones(size(mag));
for i=1:3
    ix = p(:,i) <= emsz(i) & p(:,i) >= 1; %index points inside the box
    p = p(ix,:); mag=mag(ix); %drop points outside the box
    %ixf = ixf.*ix;
    %fprintf('%i dropped \n',numel(ix)-sum(ix)) %diagnostic to check if any points eliminated
end
%p = p(ixf,:); mag=mag(ixf);
for i=1:size(p,1)
    x=p(i,1); y=p(i,2); z=p(i,3); m = mag(i); %fetch data per atom
    %split(x,y,z,j) = split(x,y,z,j)+mag; %slow, 4d indexing very inefficient
    vl(x,y,z) = vl(x,y,z)+m;
    solv(x,y,z) = solv(x,y,z)-avol*1;%(rand*.4+0.8);
end

end