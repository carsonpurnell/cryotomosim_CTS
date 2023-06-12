%% prototyping for nonflat constraints mimicking more natural ice
pix = 8;
binsize = pix*[200,100,50];
%{
%how to trim excess border length when using multiple sidedness?

edges = '  &';
vec = cell(1,numel(edges));
b = zeros(0,3);
vol = zeros(binsize/pix);
for i=1:numel(edges)
    for j=1:3
        vec{j} = 1:size(vol,j); %initialize vector triad
    end
    switch edges(i) %write the subject vector to list the appropriate border indices
        case "+", vec{i} = size(vol,i);%binsize(i);
        case "-"; vec{i} = 1;
        case "&", vec{i} = [1,size(vol,i)];
        case " ", vec{i} = []; %empty vector writes nothing when used so if check not needed
    end
    %tmp = [vec{1};vec{2};vec{3}];
    %b = [b;tmp];
    vol(vec{1},vec{2},vec{3}) = 1;
end

[x,y,z] = ind2sub(size(vol),find(vol>0));
b = [x,y,z]*pix;
d1 = 42;

n = size(b,1); ix = randperm(n); ix = ix(1:round(n/5));
b = b(ix,:);
vec = randn(size(b)); vec = d1*vec./vecnorm(vec,2,2);
b = b+vec;

bshell = repmat(b,[4,1]);
d2 = 16;
vec = randn(size(bshell)); vec = d2*vec./vecnorm(vec,2,2);
bshell = bshell+vec;

%n = size(b,1); ix = randperm(n); ix = ix(1:round(n/10));

sh = alphaShape(bshell,pix*4);
plot(sh)
%plot3(x,y,z,'.')
%plot3(qq(:,1),qq(:,2),qq(:,3),'.')
%sliceViewer(vol);
%}

%% surfgen version
pix = 10;
binsize = pix*[300,200,50]; %curvature is anisotropic, nonsquare grid has uneven noise
sz = [max(binsize),max(binsize)]; n = pix^1.2;

pts = surfgen_scripting(sz,n);

bshell = repmat(pts,[4,1]);
d2 = 12;
vec = randn(size(bshell)); vec = d2*vec./vecnorm(vec,2,2);
bshell = bshell+vec;

%sh = alphaShape(bshell,pix*5);
%plot(sh)



