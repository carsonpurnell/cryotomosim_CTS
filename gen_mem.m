function [atoms,perim,vol] = gen_mem(sz,pix,sp,thick)


arguments
    sz
    pix = []
    sp = 0.6+rand*0.4
    thick = 30
end
vol = 0;

switch 1%randi(2)
    case 1 %needs a bit more smoothing and less flat faces
        [sh,~,~,~] = blob(sz,sp); %connect and smooth scattered points
    case 2
        [sh] = bubble(sz,sp); %expand spheres from core pts
end
%alternative shape generators? cylinders, planes, sphere, stacks, double layers?
%alt gen 3: curved spline of core points with varying radii for expansion
[shell] = shape2shell(sh,thick);
[pts,head,tail] = shell2pts(shell); %need to use surface/interior separately for atomic density purposes
%better control over thickness and surface layer density - impacts layered CTF artifact a lot.
%spread of surface density also impacts the apparent thickness of final membrane, need to account

%assign atomic IDs to atoms for proper density
atoms = pts;%[head;tail];
perim = shell.Points; %perimeter from shell shape
ix = randi(size(pts,1),1,round(size(pts,1)/50)); % 2% of pts
perim = [pts(ix,1:3);perim]; perim = unique(perim,'rows');

%atoms(:,4) = 6.4*4; %terrible very bad interim density

if ~isempty(pix)
    vol = helper_atoms2vol(pix,atoms);
end

end

function [sh,pts,pts1,pts2] = blob(sz,sp)
n = round(8+sz^(0.2+sp));
rad = sz*(0+sp); var = sz*(1-sp)*2; %probably change to 1/sp-1
iters = round(1+(1+1/sp)^0.5);
az = rand(n,1)*180; el = rand(n,1)*180; r = rand(n,1)*var+rad;
[x,y,z] = sph2cart(az,el,r);
R = makehgtform('xrotate',pi/2); R = R(1:3,1:3); %get rotation matrix (3x3 of full matrix)
pts = ([x,y,z])*R; %rotate about axis so points aren't clustered in Z (stays 0-centered though)
R = makehgtform('zrotate',rand*180); R = R(1:3,1:3);
pts = pts*R; %spin about Z randomly so blobs are isotropically disordered in-plane

% surround points with randn hulls
%same strategy could generate variability inside membrane, including lipid rafts
qq = repmat(pts,round(10*(sp^0.5)),1); %replicate points by 10
d = randn(size(qq))*10/sp^2*(1-sp); %not unitized for variability?
%d = d./vecnorm(d,2,2); %vector directions, unitized

% vec = randn(size(spts));
% spd = rand(size(vec,1),1)*8+4;
% vec = vec./vecnorm(vec,2,2).*spd;
pts = qq+d;
%pts = smiter(pts,1,9);
pts = unique(pts,'rows');
sh = alphaShape(pts); sh.Alpha = criticalAlpha(sh,'one-region')+sz;

tp = randtess(0.01,sh,'s');
sh = alphaShape(tp); sh.Alpha = criticalAlpha(sh,'one-region')+sz/2;
[~,pts] = boundaryFacets(sh);

pts1 = smiter(pts,1,9); %smiter not great, can average between faces. need dist cutoff at least.
pts2 = smiter(pts,1,30);
%pts = (pts+pts1)/2;
pts1 = unique(pts1,'rows'); %prune duplicates

sh = alphaShape(pts1); sh.Alpha = criticalAlpha(sh,'one-region')+sz/2;
end

function ptsa = smiter(ptso,iter,nb)
ptsa = zeros(size(ptso));
for j=1:iter
    [~,ix] = pdist2(ptso,ptso,'euclidean','Smallest',nb);
    ix = ix(2:end,:);
    for i=1:size(ix,2)
        mm = mean(ptso(ix(:,i),:));
        ptsa(i,:) = mm;
    end
    ptso = ptsa;
end
end

function [shell] = shape2shell(shape,thick)
shellpts = randtess(thick/10.0,shape,'s'); %might be too rough at 10, smaller divisor is smoother and slower
vec = randn(size(shellpts)); vec = thick*vec./vecnorm(vec,2,2);
shellpts = shellpts+vec;
shell = alphaShape(shellpts,24);
end

function [pts,head,tail] = shell2pts(shell)
surfvar = 12;
atomfrac = 1; %make operable? 4 super rough at higher pixel sizes

tail = randtess(0.5/atomfrac,shell,'v'); %/4 for speed and memory limits. need larger hydrophobic dict
head = randtess(20/atomfrac,shell,'s'); %was 20,testing for less bilayer

vec = randn(size(head));
spd = rand(size(vec,1),1)*surfvar+0;
%spd = max(spd,0); spd = min(spd,surfvar*2); %control distant fuzzyness
vec = vec./vecnorm(vec,2,2).*spd;
head=head+vec;
pts = [head;tail];
pts(:,4) = 6.2/4 *atomfrac;
end