function con = helper_atomcon(box,pix,n,sc)
if nargin<3
    n = 4+pix^1.5;
    sc = 2400;
end
sz = [max(box),max(box)]; 
dl = 2;
w = 1;
ptsb = internal_gen_atomborder(sz,n/2,sc*1,pix*0.25)-[0,0,dl*randi(6)*w];
ptst = internal_gen_atomborder(sz,n/2,sc*1,pix*0.25)+[0,0,dl*randi(6)*w+box(3)];
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
if n==0; Z = zeros(size(X)); end

pts = [X(:),Y(:),Z(:)*sc];
n = size(pts,1); ix = randperm(n); ix = ix(1:round(n/10));
pts = pts(ix,:);
pts(:,1:2) = pts(:,1:2)-pad; %size(pts)
end