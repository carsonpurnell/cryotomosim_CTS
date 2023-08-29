function [vol,split] = helper_randfill_fil(vol,con,pix,particles,oi)
if nargin<5
    oi = zeros(1,numel(particles));
    for i=1:numel(particles)
        oi(i) = particles(i).filprop(4);
    end
end

%for i=1:numel(mono)
    namelist = [particles(:).modelname];
    for j=1:numel(namelist)
        split.(namelist{j}) = zeros(size(vol));
    end
%end
%n = 100; 
retry = 4; ori = [0,0,1];

for gg=1:numel(particles)
mono = particles(gg); iters = oi(gg); %temp before implementing internal loop
%sel = randi(numel(mono)); sel = 1;
%r = max(size(mono.sum,[1,2]))/3-4;
fpl=0;
%ang = mono.filprop(1);
step = mono.filprop(2);
%flex = mono.filprop(3);
minlength = mono.filprop(4);

for nn=1:iters
ftry=0; l=0;
while l<minlength-ftry/3 && ftry<10
    %tvol = ~(bwdist(vol)<4); %weirdly slow
    tvol = ~(vol==1);
    l = 0; fvol = vol*0; %initialize output vol
    for i=1:iters*3
    %while l<30
        for j=1:retry
            if l==0 %new start vals until initial placement found
                veci = []; rang = rand*360; %pos = rand(1,3).*size(vol); 
                pos = ctsutil('findloc',tvol); %find more reliably empty start loc
            end
            %need better pos values from bwdist by approx filament radius
            vecc = randc(1,3,veci,deg2rad(mono.filprop(3)+(j-1)*2)); 
            %generate new vector in a cone from prior vector, or any if not found
            %flexibility slightly increases with more retries to attempt filament forced bending
            pos = pos+vecc([2,1,3])*step/pix;
            
            rotax=cross(ori,vecc); rotax = rotax/norm(rotax); %compute the normal axis from the rotation angle
            theta = -acosd( dot(ori,vecc) ); %compute angle between initial and final pos (negative for matlab)
            filang = rang+mono.filprop(1)*i; %rotation about filament axis
            
            spin = imrotate3(mono.sum,filang,ori); %rotate about Z for filament twist (might go last)
            rot = imrotate3(spin,theta,[rotax(1),rotax(2),rotax(3)]); %rotate to the final position
            %would it be faster to rotate atoms and project them?
            
            com = round(pos([1,2,3])-size(rot)/2-vecc*step/pix/2);
            [~,err] = helper_arrayinsert(vol+con,rot,com,'overlaptest');
            if err==0 %place if location is good
                veci = vecc; %new initial vector for cone search to avoid high angle/retry overwrite
                l = l+1; ggg=l; %length counting
                [fvol] = helper_arrayinsert(fvol,rot,com);
                %
                for jj=1:numel(mono.modelname)
                    nm = mono.modelname{jj};
                    spin = imrotate3(mono.vol{jj},filang,ori); %rotate about Z for filament twist (go last?)
                    rot = imrotate3(spin,theta,[rotax(1),rotax(2),rotax(3)]); %rotate to the final position
                    [split.(nm)] = helper_arrayinsert(split.(nm),rot,com);
                end
                %}
                break; %early exit if good placement found
            elseif j==retry && l>minlength
                %vol = vol+fvol; fvol = fvol*0; ggg=l; l=0;
            elseif j==retry && l<minlength%-ftry/3
                fvol = fvol*0; ggg=l; l=0;
                %vol = vol+fvol; fvol=fvol*0; l=0; %ftry=ftry+1;
                %tvol = (bwdist(vol)<8);
            end
        end
        
        %{
            for j=1:1
                %mn = dat{j};
                %l = l+1; %length counting
                %[fvol] = helper_arrayinsert(fvol,rot,com);
                %if err==1, disp('ERR'); end
        %{
        tmp = dat.adat{j}; name = dat.modelname{j};
        org = [1,2,3]; %or [2,1,3] to invert xy
        tmp(:,org) = tmp(:,org)*rotmat(rotax,theta); %rotate to the filament axis, other order appears identical
        tmp(:,org) = tmp(:,org)*rotmat(vec,filax); %rotate about filament axis
        tmp(:,org) = tmp(:,org)+pos-vec/2; %move rotated unit to the target location, halfway along step
        fil.(dat.modelname{j}) = [fil.(dat.modelname{j});tmp];
        %}
            end
        %}
        if err~=0
            ftry = ftry+1; %somehow broken in certain cases, missing filament links
            if l<minlength, fvol=fvol*0; l=0; end %clear filvol if partial filament not long enough
            %if l>minlength, vol = fvol+vol; end
            %fvol=fvol*0; ggg=l; l=0;
            %fvol = zeros(size(fvol));
            %l=0; %try to re-randomize start location to avoid sequence break
            %st = strel('sphere',7); st = st.Neighborhood*1e2;
            %com = round(pos([2,1,3])-size(st)/2-vecc*step/pix/2);
            %fvol = helper_arrayinsert(fvol,st,com);
            break;  %stop filament placements if chain broken
        end
    end
    %l, %rec
end
%fprintf('%i, ',ggg);
%mask = bwlabeln(fvol>0);

CC = bwconncomp(imbinarize(fvol));
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,idx] = max(numPixels);
if ~isempty(idx)
mask = false(size(fvol));
mask(CC.PixelIdxList{idx}) = true;
fvol = fvol.*(mask>0); 
fpl=fpl+1;
else
    fvol = fvol*0; l=0;
end

vol = fvol+vol; 
%split.(mono.modelname{1}) = fvol;
fn = fieldnames(split);
for i=1:numel(fn) %loop through splits and mask out bad placements
    split.(fn{i}) = split.(fn{i}).*(vol>0);
end
%split = split.*(fvol>0);
fvol = vol*0;
end
fprintf('Placed %i filaments over %i iterations \n',fpl,iters);
end
end

%% internal functions
function vec = randc(row,col,ax,ang)
if isempty(ax), ax = randv(1,3); end %if no axis given, randomize one
nrep = row/size(ax,1); %number of replicates needed to match matrix size for cross
ax = ax/norm(ax); %unitize target vector to avoid miscalculation
rax = randv(row,col); %random axes to cross with the center axis
rotax = cross(repmat(ax,nrep,1),rax); %compute orthogonal axes to rotate 
rotax = (rotax'./vecnorm(rotax'))'; %unitize orthogonal axes
vec = zeros(row,col);
for i=1:row
    R = rotmat(rotax(i,:),rand*ang); %rotation vector
    vec(i,:) = ax*R;
end
end

function [vec] = randv(row,col)
vec = randn(row,col); %random normal numbers for evenly-distributed vector directions
vec = (vec'./vecnorm(vec'))'; %unitize vectors to length 1 for sphere vectors
end

function t = rotmat(ax,rad)
ax = ax/norm(ax);
x = ax(1,1); y = ax(1,2); z = ax(1,3);
c = cos(rad); s = sin(rad);

t1 = c + x^2*(1-c);
t2 = x*y*(1-c) - z*s;
t3 = x*z*(1-c) + y*s;
t4 = y*x*(1-c) + z*s;
t5 = c + y^2*(1-c);
t6 = y*z*(1-c)-x*s;
t7 = z*x*(1-c)-y*s;
t8 = z*y*(1-c)+x*s;
t9 = c+z^2*(1-c);

t = [t1 t2 t3
    t4 t5 t6
    t7 t8 t9];
end