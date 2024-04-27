
%% integrated polymer walk - atomistic
% working, faster than vol implementation, looks better, and doesn't have pass-through errors
% problem: won't work with membranes, which are still non-atomic. can do fil instead of mem though
%rng(11) %1 275/277/333   %2 84/109/110
profile on
pix = 10;
%input = {'MT.fil','actin.fil','cofilactin.fil','actin.fil'};
%input = {'actin.fil','actin.fil','cofilactin.fil'};
input = {'MT.fil','cofil_actin_split.fil','actin.fil','cofilactin.fil'};
%input = {'cofil_actin_split.fil','actin.fil','cofilactin.fil'};
particles = helper_filinput(pix,input);
box = [500,400,60]*pix; % box size in A

con = helper_atomcon(box,pix); %con = internal_atomcon(box,pix);
[pts,dyn,fil,mu] = helper_fil_atomic(box,particles,con);

profile viewer
[vol,solv,atlas,split] = helper_atoms2vol(pix,pts,box);
sliceViewer(vol+solv);
%%WriteMRC(vol+solv,pix,'filtest.mrc')

%requires atomistic grid and membrane though, and then projecting as a vol
%need to recheck when membrane normals are generated and if they'd break

%% integrated filament walk - vol-based version
%{
% filaments rotate in the wrong direction - XY inversion for everything?
% is this from pdb2vol? in normal cts models too
pix = 10; ori = [0,0,1];
input = 'actin_mono_fil2.cif'; prop = [-166.15,27.3,12,20];
monomeract = helper_filmono(input,pix,prop); monomeract.modelname{1} = 'actin';
input = 'MTring2.cif'; prop = [0,85,5,10];
monomer = helper_filmono(input,pix,prop); monomer.modelname{1} = 'MT';
input = 'cof_fix3.cif'; prop = [-162,24,10,20];
monomercof = helper_filmono(input,pix,prop); monomercof.modelname{1} = 'cofilactin';

%{
%input = 'actin_mono_fil2.cif'; %ang = -166.15; step = 27.3; flex = 12; minL = 20;
%input = 'MTring2.cif'; %ang = 0; step = 85; flex = 5; minL = 8;
%input = 'cof_fix3.cif'; %ang = -160; step = 24; flex = 10; minL = 15;
%dat = helper_pdb2vol(input,pix,0,1,0); 
%dat = helper_pdb2vol('MTring2.cif',pix,0,1,0); ang = 0; step = 85; flex = 5; minL=8;
%dat = helper_pdb2vol('cof_reZ2.pdb',pix,0,1,0); ang = -160; step = 24*1; flex = 10*1.0; minL = 15;
%can save with arbitrary file extensions - .fil or similar. just need to load with load(fil,'-mat')
%}
%monomer = helper_filmono(input,pix,prop);
%}
%
pix = 8;
input = {'MT.fil','actin.fil','cofilactin.fil'};
%input = 'gui'; %works, but unordered, don't have pickfiles fallback yet
particles = helper_filinput(pix,input);

%{
%dat = helper_pdb2dat(input,pix,0,1,0);
%sumv = sum(cat(4,dat{:}),4);
%need dictionary function to transform between atom Z values and scattering magnitudes
%sumv = helper_atoms2vol(pix,dat.adat,[0,0,0])*3; %lower intensity from scatter v Z
%atoms2vol output centered on 0,0,0?
%part of errors is from non-centering, so wildly wrong Z axis borks everything
%measure center and move z d models # to z-flatten things seems to fix it well enough
%minimum repeat for each filament type, maximum length? or default very overlong loop?
%
%}
%{
%r = max(size(sumv,[1,2]))/3-4; %find approximate maximum radius for bwdist comparison efficiency
%mono.vol = dat; 
mono.sum = sumv;
mono.ang = ang; mono.step = step; mono.flex = flex;%*pi/180; 
mono.minlength = minL;
%an minL be derived from step size and flexibility in some way?
%flex = flex*pi/180; %ang = ang*pi/180; %vol method is degree based
%
%}
%rng(3)
mvol = gen_memvol(zeros(500,400,50),pix,4,5)*1;
iters = [5,30,30];
con = helper_constraints(mvol*0,'  &')*pix^2.5;
vol = mvol;
%{
for nn=1:10
ftry=0; l=0;
while l<minL-ftry/3 && ftry<10
    %tvol = ~(bwdist(vol)<4); %weirdly slow
    tvol = ~(vol==1);
    l = 0; rec = 0; fvol = vol*0; %initialize output vol
    for i=1:n
    %while l<30
        for j=1:retry
            %figure out how to reproject (affine translate subpixel?) vol to right spot?
            %restart initial vals until 1st placement
            if l==0 %new start vals
                %pos = rand(1,3).*size(vol); 
                veci = []; rang = rand*360; 
                pos = ctsutil('findloc',tvol); %find more reliably empty start loc
            end
            %need better pos values from bwdist by approx filament radius
            vecc = randc(1,3,veci,flex+deg2rad((j-1)*2)); %generate new vector in a cone from prior vector, or any if not found
            %flexibility slightly increases with more retries to attempt filament forced bending
            pos = pos+vecc([2,1,3])*step/pix;
            
            rotax=cross(ori,vecc); rotax = rotax/norm(rotax); %compute the normal axis from the rotation angle
            theta = -acosd( dot(ori,vecc) ); %compute angle between initial pos and final pos (negative for matlab)
            filang = rang+ang*i; %rotation about filament axis
            
            spin = imrotate3(sumv,filang,ori); %rotate about Z for filament twist (might go last)
            rot = imrotate3(spin,theta,[rotax(1),rotax(2),rotax(3)]); %rotate to the final position
            
            com = round(pos([1,2,3])-size(rot)/2-vecc*step/pix/2);
            [~,err] = helper_arrayinsert(vol+con,rot,com,'overlaptest');
            if err==0 %place if location is good
                %rec=rec+1; 
                veci = vecc; %new initial vector for cone search to avoid high angle/retry overwrite
                l = l+1; %length counting
                [fvol] = helper_arrayinsert(fvol,rot,com);
                break; %early exit if good placement found
            elseif retry==5 && l>6
                %vol = vol+fvol; fvol=fvol*0; l=0; %ftry=ftry+1;
                %tvol = (bwdist(vol)<8);
            end 
        end
        
        %if err==0
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
            %if l<20, fvol=fvol*0; l=0; end %clear filvol if partial filament not long enough
            if l>minL, vol = fvol+vol; end
            fvol=fvol*0; l=0;
            %l=0; %try to re-randomize start location to avoid sequence break
            %st = strel('sphere',7); st = st.Neighborhood*1e2;
            %com = round(pos([2,1,3])-size(st)/2-vecc*step/pix/2);
            %fvol = helper_arrayinsert(fvol,st,com);
            break;  %stop filament placements if chain broken
        end
    end
    %l, %rec
end
disp(ftry)
vol = fvol+vol; 
fvol = vol*0;
end
%}
profile on
%[vol,split] = vol_fill_fil(vol,con,pix,particles);
[vol,split] = helper_randfill_fil(vol,con,pix,particles);

%{
%for i=1:numel(particles)
    %[vol,split] = vol_fill_fil(vol,con,pix,particles);%,iters);
%end
%
%ovol = vol_fill_fil(mvol,con,pix,monomer,iters); %it works, it's just so slow
%ovol2 = vol_fill_fil(ovol,con,pix,monomeract,iters);
%ovol3 = vol_fill_fil(ovol2,con,pix,monomercof,iters);
%}

profile viewer
sliceViewer(vol); 
%
WriteMRC(ovol3,pix,'filmixbig2.mrc')

%% integrated filament walk - atomistic version
pix = 6; ori = [0,0,1];
dat = helper_pdb2dat('actin_mono_fil.cif',pix,2,'z',0); ang = -166.15; step = 27.3; flex = 12;
%dat = helper_pdb2dat('MTring.cif',pix,2,'z',0); ang = 0; step = 85; flex = 3;
%dat = helper_pdb2dat('cofilactin_lead_samename2.cif',pix,2,'z',0); ang = -161; step = 24; flex = 10;
ang = ang*pi/180; flex = flex*pi/180; %monomer = dat.adat;
n = 50; 
pos = rand(1,3)*100; vec = []; rang = rand*360*pi/180;
%place the particles halfway from prior/next, use the line itself as the tan?
fil = struct; %initialize struct 
for i=1:numel(dat.modelname)
    fil.(dat.modelname{i}) = []; %preinitialize empty split struct
end
for i=1:n
    vec = randc(1,3,vec,flex); %generate new vector in a cone from prior vector, or any if not found
    pos = pos+vec*step;
    %d(i,:) = pos;
    
    %floc = pos;%dd(i,:); %placement centroid
    %fax = vec; fax = fax/norm(fax); %fax = [fax(2),fax(1),fax(3)];
    rotax=cross(ori,vec); rotax = rotax/norm(rotax); %compute the normal axis from the rotation angle
    theta = -acos( dot(ori,vec) ); %compute angle between initial pos and final pos (negative for matlab)
    filang = rang+ang*i; %rotation about filament axis
    
    for j=1:numel(dat.modelname)
        tmp = dat.adat{j}; name = dat.modelname{j};
        org = [1,2,3]; %or [2,1,3] to invert xy
        tmp(:,org) = tmp(:,org)*rotmat(rotax,theta); %rotate to the filament axis, other order appears identical
        tmp(:,org) = tmp(:,org)*rotmat(vec,filang); %rotate about filament axis
        tmp(:,org) = tmp(:,org)+pos-vec*step/2; %move rotated unit to the target location, halfway along step
        fil.(dat.modelname{j}) = [fil.(dat.modelname{j});tmp];
        %{
        if rem(i,2)==1
        %if j<=numel(monomer)/2%m==2
            fil.a = [fil.a;tmp];
        else
            fil.b = [fil.b;tmp];
        end
        %}
    end
    
end
%d = randc(1000,3,[0,1,1],pi/2);
%plot3(d(:,1),d(:,2),d(:,3),'o-'); axis equal
[vol,~,atlas] = helper_atoms2vol(8,fil);
volshow(atlas);

%{
%% spline generation - diagonal
x = [0,1,2,3,4,5,6,6.5];
y = [0,1,2,3,4,5.5,6,7];
z = [0,1,2,3,4,5,6.5,8];
d = [x',y',z']*100;
%% curve in xy, flat Z for inspection
x = [0,1,2,3,4,5,6,6.5];
y = [0,2,3,4,5,6,6,7];
z = x*0;
d = [x',y',z']*100;
%need to figure out how to reliably generate splines of a given flexibility
%and then also fill a region with them quickly and efficiently
%start, end, and mid point clouds to draw spline through?
%more scattered start/end for different structures. tight for bundle, 1 axis diff for flow, 2d/3d for mesh
%% random curve in mostly xy
x = rand(3,1)*20;
y = rand(3,1)*20;
z = rand(3,1)*3;
d = [x,y,z]*100;
d(2,:) = (d(2,:)+sum(d([1,3],:))*1.5)/4; %heavily weighted average to prevent high curvature
plot3(d(:,1),d(:,2),d(:,3)); axis equal
%% random curve of specific length
s = [0,0,0]; dv = rand(3,1); dv=dv/norm(dv); %start pos and vector direction
len = 2000; e = dv*len; %length and vector end
n = 5; seg = linspace(0,len,n); d = dv'.*[seg',seg',seg'];
flex = 0.5;
d = d+rand(size(d)).*(len/n*flex);
plot3(d(:,1),d(:,2),d(:,3),'o-'); axis equal
%too averaged out, doesn't have large-scale curvature. will be even worse with longer splines. 
%try starting large and subdividing segments to preserve larger-scale curves?
%make some sort of constrained walk generator for individual or bundles/meshes of filaments?
%filter new splines with pdist2(all,newtest) to pre-remove clashes
%do random walks all at the same time, or one at a time?
%bwdist to find voxels away from others?
%% random curve through step walk
%initial point, and initial vector in any direction
%then do cone random angles by replicating the prior vector, rotate about 0 on perp axis random angle<flex
%actin computed: 25.6-27.1
%cofil computed: 22.5-23.3
%step = 27.3; flex = 12; %actin
step = 25; flex = 12; %cofilactin - testing
step = 85; flex = 3; %microtubules
orig = [0,0,0]; tmp = randn(1,3); tmp = step*tmp/norm(tmp);
n = 15; d = zeros(n,3); pos = orig;
flex = pi/180*flex; %convert to radians
for i=1:n
    pos = pos+tmp;
    d(i,:) = pos;
    rax = randn(1,3);
    rotax=cross(tmp,rax/norm(rax)); %random rotation axis hopefully
    R = rotmat(rotax/norm(rotax),rand*flex);
    tmp = tmp*R;
end
plot3(d(:,1),d(:,2),d(:,3),'o-'); axis equal

%% cscvn attempt - fails due to unknowable curve distances
cv = cscvn(d'); %can use ppval(cv,t) to get the point from the length
tn = fnder(cv); %tangent derivation of spline
%may need to revert to interp1 version and construct own ppfn struct.
%distance value is unreliable, needs to change with lines of different length
%can't use one on the other - still can't convert between query scales
%opt 1: figure out how to convert into the weird space that cscvn uses
%opt 2: get the tangents from the interp1 output

ang = 166.15; ang = ang*pi/180;
step = 27.3;
pix = 5; init = [0,0,1];
dat = helper_pdb2dat('actin_monomer.pdb',pix,2,'z',0);
monomer = dat.adat{1}; fil = [];

for i=1:50
    dist = i*step/14; %distance along spline
    floc = ppval(cv,dist)'; %placement centroid
    fax = ppval(tn,dist)'; fax = fax/norm(fax); %normalized filament axis at location
    rotax=cross(init,fax); %compute the normal axis from the rotation angle
    theta = acosd( dot(init,fax) ); %compute angle between initial pos and final pos
    theta = theta*pi/180;
    t2 = ang*i;
    tmp = monomer;
    tmp(:,1:3) = tmp(:,1:3)*rotmat(rotax,theta); %rotate to the filament axis
    tmp(:,1:3) = tmp(:,1:3)*rotmat(fax,t2)+floc; %rotate about filament axis and move to filament
    fil = [fil;tmp]; %#ok<AGROW>
end
[vol] = helper_atoms2vol(pix,fil);
sliceViewer(vol);
figure
fnplt(cv); hold on
plot3(fil(:,1),fil(:,2),fil(:,3),'o'); axis equal

%{
ix = 65;
tmp = ppval(cv,ix)';
ts = ppval(tn,ix)'; %gives the vector of the tangent line from the intersect
tl = [tmp-5*ts;tmp;tmp+5*ts];
%fnplt(cv); hold on
%plot3(tl(:,1),tl(:,2),tl(:,3),'o-'); axis equal
%fnplt(cv); hold on, 
%plot3(d(:,1),d(:,2),d(:,3),'o'), hold off
%}

%% interp1 attempt - don't know exact tangent, estimating from -1/+1 positions

%ang = 35; %cofilactin - testing 2-unit - minor groove
pix = 5; init = [0,0,1];
%dat = helper_pdb2dat('actin_mono_fil.cif',pix,2,'z',0); ang = -166.15; %actin
dat = helper_pdb2dat('MTring.cif',pix,2,'z',0); ang = 0;
%monomer = dat.adat{1}; %fil{1} = []; fil{2} = fil{1}; ang = -161; %cofilactin - testing
%dat = helper_pdb2dat('cofilmonolead.pdb',pix,2,'z',0); %leading fixes weird angle errors due to centering?
ang = ang*pi/180; monomer = dat.adat;
%fil = cell(2,1);
fil.a = []; fil.b = [];
rang = rand*pi;

%{
former method of interpolating and reconstructing a spline path from a random janky line
%cs = cat(1,0,cumsum(sqrt(sum(diff(d,[],1).^2,2))));
%n = round( (cs(end)+step*1)/step); %number of subunits for polymer of the whole spline
%l = n*step; %length of spline
%dd = interp1(cs, d, linspace(0,l,n),'pchip');
%}
dd = d; %diverge for unknown reason
%axis vectors MUST be normalized. unnormed break everything
%dd = [dd(:,2),dd(:,1),dd(:,3)]; % doesn't solve plate-stacking, but changes axis it occurs
%plate stacking is happening in the Z direction without x-y inv

for i=2:n-1
    %dist = i*step; %distance along spline
    floc = dd(i,:); %placement centroid
    %fax = ppval(tn,dist)'; fax = fax/norm(fax); %normalized filament axis at location
    fax = dd(i+1,:)-dd(i-1,:); fax = fax/norm(fax); %fax = [fax(2),fax(1),fax(3)];
    rotax=cross(init,fax); rotax = rotax/norm(rotax); %compute the normal axis from the rotation angle
    theta = -acos( dot(init,fax) ); %compute angle between initial pos and final pos (negative for matlab)
    
    t2 = rang+ang*i; %rotation about filament axis
    
    for j=1:numel(monomer)
        tmp = monomer{j};
        org = [1,2,3]; %or [2,1,3] to invert xy
        tmp(:,org) = tmp(:,org)*rotmat(rotax,theta); %rotate to the filament axis, other order appears identical
        tmp(:,org) = tmp(:,org)*rotmat(fax,t2); %rotate about filament axis
        tmp(:,org) = tmp(:,org)+floc; %move rotated unit to the target location
        %m = ; 
        if rem(i,2)==1
        %if j<=numel(monomer)/2%m==2
            %fil{2} = [fil{2};tmp];
            fil.a = [fil.a;tmp];
        else
            %fil{1} = [fil{1};tmp];
            fil.b = [fil.b;tmp];
        end
    end
    %something is going wrong to severely warp the monomers at certain angles
    %was combination of rotation matrix non-unitized, inconsistent axis direction, xy confusion
%     tmp(:,1:3) = tmp(:,1:3)*rotmat(rotax,theta); %rotate to the filament axis, other order appears identical
%     tmp(:,1:3) = tmp(:,1:3)*rotmat(fax,t2)+floc; %rotate about filament axis and move to filament
    %tmp(:,1:3) = tmp(:,1:3)*rotmat(rotax,theta)+floc; %rotate to the filament axis, other order appears identical
%     m = rem(i,2);
%     if m==0
%         %fil{2} = [fil{2};tmp];
%         fil.a = [fil.a;tmp];
%     else
%         %fil{1} = [fil{1};tmp];
%         fil.b = [fil.b;tmp];
%     end
    %fil = [fil;tmp]; %#ok<AGROW>
end
[vol,~,atlas] = helper_atoms2vol(8,fil);
%sliceViewer(atlas);
volshow(atlas);
%plot3(dd(:,1),dd(:,2),dd(:,3),'.r-'); axis equal; hold on
%plot3(fil(:,1),fil(:,2),fil(:,3),'.'); axis equal

%%
%d1 = sum(sqrt(diff(x).^2+diff(y).^2+diff(z).^2));
%t = [1,2,3,4,5,6,7,8]; % Assumed time stamp
%plot3(x,y,z); axis equal

% tt = linspace(0,d1);
% xx = interp1(t,x,tt,'spline');
% yy = interp1(t,y,tt,'spline');
% zz = interp1(t,z,tt,'spline');
% figure
% scatter3(x,y,z)
% hold on
% plot3(xx,yy,zz); axis equal

d = [x',y',z']*100;
cs = cat(1,0,cumsum(sqrt(sum(diff(d,[],1).^2,2))));
n = round( (cs(end)+step*1)/step); %number of subunits for polymer of the whole spline
l = n*step; %length of spline
dd = interp1(cs, d, linspace(0,l,n),'pchip'); %unique([CS(:)' linspace(0,l,n)])
%tangents
% m = gradient(dd); %slope of tan at interp pts?
% ix = 44;
% tmp = dd(ix,:);
% tl = [tmp+m(ix,:)*2;tmp;tmp-m(ix,:)*2];


%figure, hold on
fnplt(cscvn(d')); hold on
%plot3(d(:,1),d(:,2),d(:,3),'.b-')
%plot3(m(:,1),m(:,2),m(:,3),'.b-')
%plot3(tl(:,1),tl(:,2),tl(:,3),'.b-')
plot3(dd(:,1),dd(:,2),dd(:,3),'.r-'); axis equal
%axis image, view(3), legend({'Original','Interp. Spline'})

%% 2 from monomer filament generation
ang = 166.15; ang = ang*pi/180;
step = 27.3;
pix = 2; 

dat = helper_pdb2dat('actin_monomer.pdb',pix,2,'z',0);
monomer = dat.adat{1};

bnd = [];
for i=1:3
    %filament start adjustment
    slide = [90*i,0,0,0];
    %poly = monomer+slide;
    for j=1:randi(20)+5
        tmp = monomer;
        %polymerization loop
        R = rotmat([0,0,1],ang*j); %rotation matrix assembler from imrotate3
        tmp(:,1:3) = tmp(:,1:3)*R+[80*i,0,step*j];
        bnd = [bnd;tmp];
    end
end
[vol,solv,atlas,splitvol] = helper_atoms2vol(6,bnd);
sliceViewer(vol);

%% 1 initial filament polymerizer test script, prior to curved filament generator
% actin: 166.15 pitch, 27.3 rise almost exact
ang = 166.15; ang = ang*pi/180;
step = 27.3;
pix = 2; 
clear particles;
input = {'actinpolyalignedunbork.cif','actin_monomer.pdb'};
for i=1:numel(input)
    particles(i) = helper_pdb2dat(input{i},pix,2,1,0); %load test unit
end
particles(2) = helper_pdb2dat(input{i},pix,2,'z',0);
actin = particles(1);
mono = particles(2).adat{1}-[0,0,step*22,0];
%actin oriented along the Z dimension
%actmono = actin.adat{1}; %doesn't work because it's not a real monomer, but a replimer along the chain
actolig = cat(1,actin.adat{:});
%% 

bundle = actolig;
for i=1:3
    tmp = actolig;
    %{
    R = [cos(ang*i) -sin(ang*i) 0; ...
        sin(ang*i)  cos(ang*i) 0; ...
        0           0  1];
    %}
    R = rotmat([0,0,1],ang*i); %rotation matrix assembler from imrotate3
    tmp(:,1:3) = tmp(:,1:3)*R+[80,0,step]*i;
    %{
    tt = makehgtform('zrotate',ang*i);
    tform = affine3d(tt);
    tmp(:,1:3) = transformPointsForward(tform,tmp(:,1:3))+[80,0,step]*i;
    %}
    
    olig2 = tmp;
    bundle = [bundle;olig2];
end
poly = mono;
for i=1:50
    tmp = mono;
    R = rotmat([0,0,1],ang*i); %rotation matrix assembler from imrotate3
    tmp(:,1:3) = tmp(:,1:3)*R+[0,0,step]*i;
    %tt = makehgtform('zrotate',ang*i);
    %tform = affine3d(tt);
    %tmp(:,1:3) = transformPointsForward(tform,tmp(:,1:3))+[0,0,step]*i;
    poly = [poly;tmp];
end
bundle = [bundle;poly+[0,80,0,0]];


[vol,solv,atlas,splitvol] = helper_atoms2vol(5,bundle);
sliceViewer(vol);
%WriteMRC(vol+solv,pix,'upscaletest_5.mrc');
%}

%% internal functions
function con = internal_atomcon(box,pix,n,sc)
if nargin<3
    n = 4+pix^1.5;
    sc = 2400;
end
sz = [max(box),max(box)]; 
dl = 10;
w = 1;
ptsb = internal_gen_atomborder(sz,n/2,sc*1,4)-[0,0,dl*randi(6)*w];
ptst = internal_gen_atomborder(sz,n/2,sc*1,4)+[0,0,dl*randi(6)*w+box(3)];
pts = zeros(0,3);
for i=1:4
    pts = [pts;ptsb-[0,0,dl*i]]; pts = [pts;ptst+[0,0,dl*i]];
end
con = pts;
end
function pts = internal_gen_atomborder(sz,n,sc,sep)
%n = 2.5; % noise magnitude
%sc = 500; % scale of Z noise
%sep = 3;
pad = 20; %padding - scale by input size maybe? prune afterward?
sd = max(sz)+pad*2;
[X,Y] = ndgrid(1:sep:sd,1:sep:sd);
i = min(X-1,sd-X+1); j = min(Y-1,sd-Y+1);
H = exp(-.5*(i.^2+j.^2)/n^2);
Z = real(ifft2(H.*fft2(randn(size(X))))); % 0-centered, approximately normal

pts = [X(:),Y(:),Z(:)*sc];
n = size(pts,1); ix = randperm(n); ix = ix(1:round(n/10));
pts = pts(ix,:);
pts(:,1:2) = pts(:,1:2)-pad;
%size(pts)
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

function [vol,split] = vol_fill_fil(vol,con,pix,particles,oi)
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
retry = 5; ori = [0,0,1];

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
    for i=1:iters*2
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
            elseif retry==5 && l>minlength
                %vol = vol+fvol; fvol = fvol*0; ggg=l; l=0;
            elseif retry==5 && l<minlength%-ftry/3
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

function vec = randc(row,col,ax,ang)
if isempty(ang) || ang==0
    vec = randv(row,col); return %if no angle given, random whole circle vectors
end
if isempty(ax), ax = randv(1,3); end %if no axis given, randomize one
if numel(ang)==1, ang(2)=ang(1); ang(1)=0; end %if only one angle given use 0 as minimum
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
function t = rot_aff(ax,rad)
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

t = [t1 t2 t3 0
    t4 t5 t6 0
    t7 t8 t9 0
    0  0  0  1];
end