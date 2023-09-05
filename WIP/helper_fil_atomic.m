function [pts,dyn] = helper_fil_atomic(box,particles,con)



ori = [0,0,1]; tol = 2;
if isempty(con)
    dyn{1} = zeros(0,3);
else
    dyn{1} = con;
end
dyn{2} = 0; retry = 5;
mn = [particles.modelname]; %round up all names for models
for i=1:numel(mn)
    pts.(mn{i}) = zeros(0,4);
end

for ol=1:numel(particles)
    iters = 0.2*particles(ol).filprop(4)^2;
    ct = 0;
for i=1:iters
    mono = particles(ol);
    %step = mono.filprop(2);
    %ml = mono.filprop(4);
    l=0; fail=0; fil = struct; END=0;
    for j=1:mono.filprop(4)*10
        if END==1 || fail==1; fprintf('bail,'); break; end %never bails here for some reason
        if fail==0 && END==0
        for il=1:retry
            if l==0 %new start vals until initial placement found
                veci = randc(1,3,ori,deg2rad([60,120])); %random in horizontal disc for efficiency
                rang = rand*360; pos = rand(1,3).*box; %prob need better in-box randomizing
                for mmm=1:numel(mono.modelname)
                    fil.(mono.modelname{mmm}) = zeros(0,4);
                end
            end
            
            vecc = randc(1,3,veci,deg2rad(mono.filprop(3)+(il-1)*2)); %random deviation vector
            pos = pos+vecc([1,2,3])*mono.filprop(2); %0-centered placement location from vector path
            
            if any(pos+200<0) || any(pos-200>box)
                err = 1; %if pos is too far out of box, retry without pointless proxtesting
            else
                rotax=cross(ori,vecc); rotax = rotax/norm(rotax); %compute the normal axis from the rotation angle
                theta = -acos( dot(ori,vecc) ); %compute angle between initial and final pos (negative for matlab)
                filang = rang+mono.filprop(1)*j; %rotation about filament axis
                
                r1 = rotmat(rotax,theta); r2 = rotmat(vecc,deg2rad(filang));
                com = pos-vecc*mono.filprop(2)/2;
                ovcheck = particles(ol).perim*r1*r2+com;
                err = proxtest(dyn{1},ovcheck,tol);
            end
            
            if err==0, break; end %if good placement found, early exit
            if err~=0 && il==retry, END=1; fail=1; end %is reporting here %, fprintf('c%i,',i)
        end
        else
            fprintf('avoided '); %never reaches this point? fail and end both don't work at all?
        end
        
        if err~=0, END=1; fail=1; end
        %{
        if err==1 || END==1
            fail = 1;
            fprintf('t')
            break; %this break does not work
            %}
        if err==1 || fail==1 || END==1
            %fprintf('a%i,',i) %this is the only break that gets hit
            %break
        elseif err==0 && fail==0 && END==0
            for iix=1:numel(mono.adat) %loop through and cumulate atoms
                tmp = mono.adat{iix}; %fetch atoms, needed to operate on partial dimensions
                org = [1,2,3]; %or [2,1,3] to invert xy
                %tmp(:,org) = tmp(:,org)*r1; %rotate to the filament orientation
                %tmp(:,org) = tmp(:,org)*r2; %rotate about the filament axis
                tmp(:,org) = tmp(:,org)*r1*r2+com; %move to halfway along current vector
                fil.(mono.modelname{iix}) = [fil.(mono.modelname{iix});tmp];
            end
            %fprintf('%i,',fail) %fail = 0;
            veci = vecc; l=l+1; %store current vector as prior, increment length tracker
        else%if il==retry
            %{
            for mmm=1:numel(mono.modelname)
                    fil.(mono.modelname{mmm}) = zeros(0,4);
            end
            %}
            %fail = 1;
            fprintf('b%i,',i) %now never reaches, caught by first check
            break %this break still seems to fail sometimes. pass-through filaments still happen
        end
    end
    
    if l>mono.filprop(4)*1.0 %&& fail==0
    fn = fieldnames(fil); fail = 0; pos = []; l=0;
    for fsl=1:numel(fn)
        pts.(fn{fsl}) = [pts.(fn{fsl});fil.(fn{fsl})]; 
        n = size(fil.(fn{fsl}),1);
        ix = randperm(n); ix = ix(1:round(n/40));
        lim = fil.(fn{fsl})(ix,1:3);
        dyn{1} = [dyn{1};lim];
        %fil.(fn{fsl}) = zeros(0,4);
    end
    end
    %fil = struct; l=0;
    if fail==0; ct = ct+1; end
end
fprintf('Layer %i, placed %i out of %i \n',ol,ct,iters)
end


end

function err = proxtest(c,pts,tol)
l = min(pts,[],1)-tol; h = max(pts,[],1)+tol; %low and high bounds per dimension
ix = c>l & c<h; % compare all points against the prospective box
ix = all(ix,2); % filter to index of pts inside the box
if ~any(ix), ix=[]; end % check for early end if no points in the box
%ix = find(ix>0); %bottleneck - just too many points. mutable octree should be faster overall
err=0; %with n=100 exhaustive is only slightly slower than kdtree search, but progressive slowdown
if ~isempty(ix) %this thing is taking SO VERY LONG, need more pre-optimization
    buck = 100;%round( size(c,1)/7650 ); %very rough, is probably not linear scale
    % probably needs some sort of depth-based metric, not a flat one depth = log2 (n/leaf)
    %ot = OcTree(c(ix,:),'binCapacity',buck); %significantly slower than kdt build
    
    %mutree = octcubetree(c(ix,:),'leafmax',500); %slightly faster than kd building
    %err = mutreetest(mutree,pts); %WAYYY slower than knn search
    %err = any(err);
    
    modeltree = KDTreeSearcher(c(ix,:),'Bucketsize',buck); %67 with 1K %32 with 10K, 18 100K
    [~,d] = rangesearch(modeltree,pts,tol,'SortIndices',0); %?? 1K,11.4 10K, 85 100K
    d = [d{:}]; if any(d<tol), err=1; end %test if any points closer than tol
end
end

function vec = randc(row,col,ax,ang)
if isempty(ax), ax = randv(1,3); end %if no axis given, randomize one
%if ang is 1x2, use 1st for min 2nd max?
if numel(ang)==1, ang(2)=ang(1); ang(1)=0; end
ang(2) = ang(2)-ang(1); %store difference from min for simpler following code
%ang is IN RADIANS
nrep = row/size(ax,1); %number of replicates needed to match matrix size for cross
ax = ax/norm(ax); %unitize target vector to avoid miscalculation
rax = randv(row,col); %random axes to cross with the center axis
rotax = cross(repmat(ax,nrep,1),rax); %compute orthogonal axes to rotate 
rotax = (rotax'./vecnorm(rotax'))'; %unitize orthogonal axes
vec = zeros(row,col);
for i=1:row
    R = rotmat(rotax(i,:),rand*ang(2)+ang(1)); %rotation vector
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