% testing thickness computation from surface points
pix = 10;
box = [400,300,40]*pix;
angles = -60:5:60;
%rng(7)
mang = max(abs(angles));
mtilt = tand(mang)*0.6; % angle for initial grid padding to cover tilt area - don't know why 0.6 is enough
% make surfaces - displaced in Z
% trying with identical surf pairs first
res = pix/2; %oversample to prevent holes in data
%need to stretch coverage a lot, fill values are static and not NN or average of near edge
pad = [pix*2,pix*2]; axis = 1; % huge padding only against tilt axis, otherwise small
pad(axis) = round(box(axis)*mtilt); % huge padding to cover for tilting
[x,y] = meshgrid(pix/2-pad(1):res:box(1)+pad(1),pix/2-pad(2):res:box(2)+pad(2));
%{
%surfinit = [x(:),y(:)];
%n = pix; sc = 4800; sep = 0.5;

%pts = border(max(box),n,sc,sep);
%pts = border2(box,pix);
%definition: s{1} is the top/first E interaction, s{2} is the bottom/last interaction
%sampsurf(1:2) = {[x(:),y(:)]}; 
%s1 = surfinit; s2 = surfinit;
%sampsurf{1}(:,3) = box(3)/2+sin(x(:)/90).*sqrt(x(:)+200);
%sampsurf{2}(:,3) = -box(3)/2;%+randn(size(surfinit,1),1)*5;%+sin(x(:)/90).*sqrt(x(:));
%noise = rand(size(x))*9-5; sm = imgaussfilt(noise,27,'FilterSize',171)*40;
%sampsurf{2} = [x(:),y(:),sm(:)];
%sampsurf{1} = border(max(box),n*2,sc*2,pix/2)+[0,0,box(3)/2];
%sampsurf{2} = border(max(box),n,sc,sep)-[0,0,box(3)/2];
%plot3p([s1;s2],'.');
%sampsurf{1} = pts+[0,0,box(3)/2*1.1];

%pts1 = surfice(size(x),pix)*sc(1); %pts = rescale(pts,0,100);
%sampsurf{1} = [x(:),y(:),pts1(:)+box(3)*sc(2)/2];
%w = box(3)/2; scale = 0.5;
%pts2 = surfice(size(x),pix)*sc(1); %pts = rescale(pts,w,w+w*scale);
%sampsurf{2} = [x(:),y(:),-pts2(:)-box(3)*sc(2)/2];
%pts = znoise(x)*sc(1);
%[pts,wn] = colored_noise(max(size(x)),2,-2); pts = pts(1:size(x,1),1:size(x,2))*sc(1);
%sampsurf{1} = [x(:),y(:),pts(:)+box(3)*sc(2)/2];
%both methods are inconsistent noise. larger grids smooth things out, need a way to counteract
%im = perlin_noise(x);
%im = perlin(size(x),'correl',0.1); im = rescale(im,-1,1)*sc(1);
%unfortunately quite slow due to stepwise interpolations
%[im,l] = perlin2D(x,pix); %sliceViewer(l); %histogram(im)
%[field,layers] = helper_perlin(x,pix,2.5,9,8);
%mv = mean(field,'all')*0.5; if mv<0; mv=mv*2-0*sqrt(abs(mv)); end
%field = field-mv; %hold on; %histogram(im)
%im = rescale(im,-1,1)*sc(1);
%sampsurf{1} = [x(:),y(:),field(:)+box(3)*sc(2)/2];

% rotate both surfaces by the same angle about 0
%sampsurfrot = rotsurf(sampsurf,box,pi/15,[0,1,0]);
% s1rot = sampsurf{1}; %s2rot = sampsurf{2}; %initial test: no rotation!
% cen = [box(1)/2,0,-box(3)/2*0];
% R = rotmat([0,1,0],-pi/8);
% %R2 = makehgtform('xrotate',pi/12); R3 = makehgtform('yrotate',pi/12);
% s2rot = (sampsurf{2}-cen)*R+cen; % add z displacement at the end? not correct math though

% re-grid surfaces for pair comparisons - round them to target grid?
%s1rot = prune(sampsurfrot{1},box); s1regrid = regrid(s1rot,box,pix);
%s2rot = prune(sampsurfrot{2},box); s2regrid = regrid(s2rot,box,pix);

%output pruned 0* surfaces for viewing?

%thick = s1regrid-s2regrid;
%pts = surfice(box);
% s1regrid = [round(s1rot(:,1:2)/10),s1rot(:,3)];
% s1regrid = prune(s1regrid,box/pix); %prune first because accum can't handle out-of-bounds
% s1mat = accumarray(s1regrid(:,1:2),s1regrid(:,3),box(1:2)/pix,@mean,mean(s1regrid(:,3)));
%}

sc = [2.5,1.2];
sampsurf{1} = gensurf(x,y,box,pix,sc);
sampsurf{2} = gensurf(x,y,box,pix,-sc);
[tilts,gridt] = thicktilts(sampsurf,box,pix,angles); sliceViewer(tilts);

function gridsurf = gensurf(x,y,box,pix,sc)
[field,layers] = helper_perlin(x,pix,sc(1),8,8);
mv = mean(field,'all')*0.5; 
if mv<0&&sc(2)>0; mv=mv*2-0*sqrt(abs(mv)); end
if mv>0&&sc(2)<0; mv=mv*2-0*sqrt(abs(mv)); end
field = field-mv;
gridsurf = [x(:),y(:),field(:)+box(3)*sc(2)/2];
end

function [s,l] = perlin2D(m,pix)
s = 0;%zeros(size(m));     % Prepare output image (size: m x m)
w = size(m);
i = 8-round(log2(pix));
e = i+12;
l = zeros(0);
j = i:e;%,e,e]
pad = 10;
m = zeros(size(m)+pad);
%i = 4;
%change loop to a resolution-based metric - wavelength or frequency?
for jj=j%while 2 > 1 && i<e
    %i = i + 2;%1;
    i=jj;
    d = randn(size(m));
    for k=1:i-1 %do iterative refinements and shrinking rather than 2^k expansion in one step
        d = interp2(d, 1, 'spline');
        d = d(1:size(m,1), 1:size(m,2));
    end
    %d = rescale(d,-1,1);
    l(:,:,end+1) = (1.35^i) *2* d(1:size(m,1), 1:size(m,2));
    s = s + l(:,:,end);
    %w = w - ceil(w/2 - 1);
end
s = s(pad+1:pad+w(1),pad+1:pad+w(2));
l = l(:,:,2:end);
%s = (s - min(min(s(:,:)))) ./ (max(max(s(:,:))) - min(min(s(:,:))));
end

function [tilts,regridt] = thicktilts(sampsurf,box,pix,angles)
tilts = zeros(box(1)/pix,box(2)/pix,numel(angles));
regridt(1:2) = {zeros(box(1)/pix,box(2)/pix,numel(angles))};%cell(1,2);
for i=1:numel(angles)
    sampsurfrot = rotsurf(sampsurf,box,deg2rad(angles(i)),[0,1,0]);
    %s1rot = prune(sampsurfrot{1},box); 
    regridt{1}(:,:,i) = regrid(sampsurfrot{1},box,pix);
    %s2rot = prune(sampsurfrot{2},box); 
    regridt{2}(:,:,i) = regrid(sampsurfrot{2},box,pix);
    thick = regridt{1}(:,:,i)-regridt{2}(:,:,i); 
    tilts(:,:,i) = thick;
end

end

function surfcell = rotsurf(surfcell,box,theta,ax)
for i=1:2
    %s1rot = surfcell{i}; %s2rot = sampsurf{2}; %initial test: no rotation!
    cen = [box(1)/2,0,0]; %not appropriate for X tilt
    R = rotmat(ax,theta);
    %R2 = makehgtform('xrotate',pi/12); R3 = makehgtform('yrotate',pi/12);
    surfcell{i} = (surfcell{i}-cen)*R+cen;
end
end

function pts = regrid(pts,box,pix)
%tmp = [round(pts(:,1:2)/pix),pts(:,3)];
p = round(pts(:,1:2)/pix); v = pts(:,3);
[p,v] = prune2(p,v,box/pix);
%tmp = prune(tmp,box/pix); %prune first because accum can't handle out-of-bounds
%tmp = sortrows(tmp); %takes more time than it saves from accumarray
%pts = accumarray(tmp(:,1:2),tmp(:,3),box(1:2)/pix,@mean);%,mean(tmp(:,3))); %mean 300x slower
pts = accumarray(p,v,box(1:2)/pix)./accumarray(p,1,box(1:2)/pix);
%pts2 = accumarray(tmp(:,1:2),tmp(:,3),box(1:2)/pix)./accumarray(tmp(:,1:2),1,box(1:2)/pix);
%pts = sm./c;
%mint = mean(pts(pts>0),'all'); pts(pts==0)=mint;
end

function [p,v] = prune2(p,v,box)
p = real(p);
for i=1:2
    ix = p(:,i) <= box(i) & p(:,i) >= 1; %index points inside the box
    p = p(ix,:); v = v(ix); %drop points outside the box
end
end

function pts = prune(pts,box)
pts = real(pts);
for i=1:2
    ix = pts(:,i) <= box(i) & pts(:,i) >= 1; %index points inside the box
    pts = pts(ix,:); %drop points outside the box
end
end

function t = rotmat(ax,rad)
ax = ax/norm(ax);
x = ax(1); y = ax(2); z = ax(3);
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

function pts = surfice(grid,pix)
k = max(grid)*1;%/pix*1; % size of the matrix, must be even and square
%m = randn(k); % k-by-k matrix
mf = fftshift(fft2(randn(k)));
%a = 2; % Parameter 3. 1/f? exponent
d = ((1:k)-(k/2)-1).^2;
dd = sqrt(d(:) + d(:)');
filt = dd .^ -(2);%a; % frequency filter
filt(isinf(filt))=1; % replace +/-inf at DC or zero-frequency component
ff = mf .* filt;
b = ifft2(ifftshift(ff));
pts = rescale(real(b),-1,1);
pts = pts(1:grid(1),1:grid(2));
end

function im = perlin_noise(im)
[n, m] = size(im);
i = 0;
w = sqrt(n*m);

while w > 3
    i = i + 1;
    d = interp2(randn(n, m), i-1, 'spline');
    im = im + i * d(1:n, 1:m);
    w = w - ceil(w/2 - 1);
end
end

function pts = border2(box,pix)
n = pix/2; sc = 4800*5; %sep = 1;
oversample = 1;
%need to stretch coverage a lot, fill values are static and not NN or average of near edge
pad = max(box)/1; %massive pad value to fill holes from high rotations
sd = max(box)+pad;
n = sd/200;
%should only pad for the direction lost to tilting
[x,y] = meshgrid(-pad:pix/oversample:sd,-pad:pix/oversample:sd);

%n = 2.5; % noise magnitude
%sc = 2000; % scale of Z noise
%pad = 20; %padding - scale by input size maybe? prune afterward?
%x = pts(:,1); y = pts(:,2);

i = min(x-1,sd-x+1); j = min(y-1,sd-y+1);
H = exp(-.5*(i.^2+j.^2)/n^2);
z = real(ifft2(H.*fft2(randn(size(x))))); % 0-centered, approximately normal

pts = [x(:),y(:),z(:)*sc];
end
function pts = border(sz,n,sc,sep)
%n = 2.5; % noise magnitude
%sc = 500; % scale of Z noise
%sep = 3;
pad = 20; %padding - scale by input size maybe? prune afterward?
sd = max(sz)*3;
[X,Y] = ndgrid(1:sep:sd,1:sep:sd);
i = min(X-1,sd-X+1); j = min(Y-1,sd-Y+1);
H = exp(-.5*(i.^2+j.^2)/n^2);
Z = real(ifft2(H.*fft2(randn(size(X))))); % 0-centered, approximately normal

pts = [X(:),Y(:),Z(:)*sc];
%n = size(pts,1); ix = randperm(n); ix = ix(1:round(n/10));
%pts = pts(ix,:);
pts(:,1:2) = pts(:,1:2)-max(sz)/1; %size(pts)
end

function [sig,n]=colored_noise(sz,dim,exponent)
n=squeeze(normrnd(0,1,[ones(1,dim)*sz 1]));
dist=ones(sz,1); % create Fourier filter
for i=1:sz
    dist(i)=norm(i-sz/2-1)^2;
end
dist_tot=dist;
for d=2:dim
    dist_add=reshape(dist,[ones(1,d-1) sz]);
    dist_tot=dist_tot+dist_add;
end
dist_tot=sqrt(dist_tot);
filt=dist_tot.^(exponent*dim/2);
filt(isinf(filt))=1; % leave DC unchanged
wnf = fftshift(fftn(n)); % Fourier transform white noise, then fftshift to shift 0-frequency
wnf_filt = wnf.*filt; % multiply with frequency filter
sig = ifftn(ifftshift(wnf_filt)); % ifftshift to first shift back the Fourier transform
sig = rescale(real(sig),-1,1);
end

function outpict = perlin(outsize,varargin)
%   PERLIN(SIZE, {OPTIONS})
%       Generates pseudo-perlin noise fields or volumes. Compared to PERLIN3, 
%       behavior is nonrepeatable, but the execution is much faster.  Instead 
%       of generating a perlin noise volume to create a multichannel image, 
%       PERLIN generates a moving weighted sum of 2-D noise sets.  A potential 
%       benefit of this compromise to avoid the expense of volumetric interpolation 
%       is that page correlation is uniquely controllable as a parameter. 
%   
%   SIZE is a 2 or 3-element vector defining the size of the output image
%       Dim 3 is not limited to the sizes expected for typical image channel formats
%   OPTIONS includes the following key-value pairs:
%       'correl' specifies the page or channel correlation factor (default 1)
%           This is a proxy for the correlation coefficient between pages (image
%           channels).  For 0, pages are uncorrelated.  For the default value, 
%           page correlation is approximately 93%, more closely approximating the 
%           appearance of 3-D noise.  Values above 1 further increase correlation.
%           In the context of an RGB image, CORREL=0 yields a garish rainbow cloud.  
%           As CORREL increases beyond 1, the image approaches grayscale.
%       'interpolation' specifies the interpolation used in scaling noise components
%           supports the methods used by interp2() (default 'spline')
%       'outclass' specifies the output image class (default 'double')
%           supports all standard image class names
%
% Webdocs: http://mimtdocs.rf.gd/manual/html/perlin.html 
% See also: perlin3
correl = 1; 
outclassstrings = {'double','single','uint8','uint16','int16','logical'};
outclass = 'double';
interpmethodstrings = {'nearest','linear','cubic','spline'};
interpmethod = 'spline';
if numel(varargin) > 0
	k = 1;
	while k <= numel(varargin)
		switch lower(varargin{k})
			case 'correl'
				if isnumeric(varargin{k+1})
					correl = varargin{k+1};
				else
					error('PERLIN: expected numeric value for CORREL')
				end
				k = k+2;
			case 'outclass'
				thisarg = lower(varargin{k+1});
				if ismember(thisarg,outclassstrings)
					outclass = thisarg;
				else
					error('PERLIN: unknown output class %s\n',thisarg)
				end
				k = k+2;
			case 'interpolation'
				thisarg = lower(varargin{k+1});
				if ismember(thisarg,interpmethodstrings)
					interpmethod = thisarg;
				else
					error('PERLIN: unknown interpolation method %s\n',thisarg)
				end
				k = k+2;
			otherwise
				error('PERLIN: unknown input parameter name %s',varargin{k})
		end
	end
end
pagesize = outsize(1:2);
s = max(pagesize,2);
if numel(outsize) == 3
    numchan = outsize(3);
else
    numchan = 1;
end
outpict = zeros([pagesize numchan]); 
for c = 1:numchan
    wpict = zeros(s);
    w = max(s);
    k = 0;
    while w > 3
        k = k+1;
        d = interp2(randn(ceil(s/(2^(k-1))+1)), k-1, interpmethod)*k;
        wpict = wpict + k*d(1:s(1),1:s(2));
        w = w-ceil(w/2 - 1);
    end
    if sum(pagesize > 1) == 1
        wpict = wpict(1:pagesize(1),1:pagesize(2));
    end
    
    outpict(:,:,c) = wpict;
end
if numchan > 1 && correl ~= 0
	% this is an attempt to counteract the natural scale-dependence of correlation in these sums
	% for volumes, this may be undesired, but for RGB images and typical use, it's probably for the best
	dv = log10(1E6/prod(pagesize))/(correl*300);
	for c = 2:numchan
		outpict(:,:,c-1) = simnorm(outpict(:,:,c-1),'mean')-0.5;
		outpict(:,:,c) = outpict(:,:,c-1)+outpict(:,:,c)*dv;
	end
	outpict(:,:,c) = simnorm(outpict(:,:,c),'mean')-0.5;
	outpict = outpict+0.5;
else
	for c = 1:numchan
		outpict(:,:,c) = simnorm(outpict(:,:,c),'mean');
	end
end
%outpict = imcast(outpict,outclass);
end
function out = simnorm(inpict,constraint)
%[mn av mx] = imstats(inpict,'min','mean','max');
mn = min(inpict,[],'all'); mx = max(inpict,[],'all'); 
av = mean(inpict,'all');
os = max(abs(av-mn),abs(av-mx));
out = (inpict-av)/(2*os) + 0.5;
end