function [err,ix] = mu_search(mu,test,tol,opt)

arguments
    mu
    test (:,3)
    tol = 2 %proximity tolerance for overlap/error
    opt.short = 1 %default short circuit - ss ends after first error and returns only 1 for error
    %opt.pdist = 0 %default vecnorm per point, flag for using vectorized pdist2 with unique
    %pdist should always be faster than vecnorm if not short-circuit. do binary toggle?
end
%short = strcmp(opt.type,'short');
test = single(test); %convert to single so pdist2 doesn't crap itself
bm = zeros(1,3,8); for i=1:8; bm(1,:,i) = bitget(i,1:3); end %generate bitmask

n = size(test,1);
err = zeros(1,n); %array for overlap results per point
mdepth = size(mu,1); % current depth of mutree

ix = rootsplit(mu,test); %125
%{
%170 inlined - subfunct faster due to cleanup?
ix = ones(2,n); %initialization for index
roots = vertcat(mu{1,1}{2,:}); %all top-level centers
%check for multiple roots to avoid slow pdist 2 check?
%D = distEucSq(roots,test); [~,rix] = min(D,[],1); % slower under most conditions
[~,rix] = pdist2(roots,test,'squaredeuclidean','Smallest',1);
%[D] = mypdist2(roots,test); %would need to filter afterward - and too slow anyway
ix(1,:) = 1; ix(2,:) = rix; %level and layer to start with
%}

%histogram(rix)
%{
if ~isempty(mutree{1,1}{3,1})
cmp = test>mutree{1,1}{2,1}; %get split bin
[~,bin] = max(all(bsxfun(@eq,cmp,bm),2),[],3); %convert bit comparison to numeric branch label
ix(1,:) = 2; ix(2,:) = bin;
end
%}

%loop through binned prelim branches first to take adv of fast pdist2 vectorization?
%bsxfun nav simply not fast when going through so many individual points
%bsxfun apply to all prelim leaf points?
if mdepth>1
    [uix] = fastunique(ix');
    %[uixo] = unique(ix','rows'); %[uixo,iao,ico] = unique(ix','rows');
    %if ~all(all(uixo==uix)) || ~all(all(iao==ia)) || ~all(all(ico==ic)); disp('e'); end
    for i=1:size(uix,1)
        ix = rsplit2(mu,test,ix,uix(i,:));
        %{
        d = uix(i,1); br = uix(i,2); %start bin to work with
        if numel(mu{d,1}{3,br})==8 %if no points, don't nav
            ix2 = all(ix==(uix(i,:)'),1); %logical index of pts in the root to split
            tmp = test(ix2,:); %which pts being worked with
            
            leaves = mu{d,1}{3,br}; %indices of leaves from branch
            bcp = vertcat(mu{d+1,1}{2,leaves}); %centers for child nodes for pdist2
            %mu{d+1,1}{3,leaves}
            [~,c2] = pdist2(bcp,tmp,'squaredeuclidean','Smallest',1); %which branch for the bin pts
            %D = distEucSq(bcp,tmp); [~,cc] = min(D,[],1); %slower
            t2 = zeros(2,numel(c2))+d+1;
            t2(2,:) = leaves(c2);
            ix(:,ix2) = t2;%[d+1,c2]';
        end
        %}
    end
end

%loop through points
for i=1:n
    %ix(:,i)
    %size(mu{ix(1,i),1}{1,ix(2,i)},2); %size(mu{ix(1,i),1}{3,ix(2,i)},1)
    if ix(1,i)==mdepth %if bottom of tree no nav
        %it = ix(:,i);
    elseif size(mu{ix(1,i),1}{1,ix(2,i)},2)==3 %if points in leaf no nav
        %it = ix(:,i);
    %elseif size(mu{ix(1,i),1}{3,ix(2,i)},1)==2 %if corners no nav
        %it = ix(:,i); %appears redundant with above, should always be pts if there are corners
    else
        ix(:,i) = nav(mu,bm,test(i,:),ix(1,i),ix(2,i)); %still too common, bin check not deep?
    end
    if opt.short==1
        te = prox(mu,test(i,:),ix(:,i),tol);
        if opt.short %short circuit loop breaker, return 0
            if te, err=i; %disp('kill');
                return; end
        end
        err(i) = te;
    end
end
if opt.short==0 %pdist 2 faster for searching few bins with many points, short-circuit not as powerful
    [err] = bincheck(mu,test,ix);
    %{
    [branchdet,~,bix] = fastunique(ix');
    %[branchdeto,~,bixo] = unique(ix','rows'); %more bins = better ss, but slower pdist2
    %if ~all(all(bix==bixo)) || ~all(all(branchdet==branchdeto)); disp(1); end
    d2 = zeros(size(ix,2),1)+tol^2;
    for i=1:size(branchdet,1)
        bpts = mu{branchdet(i,1),1}{1,branchdet(i,2)};
        if ~isempty(bpts)%, err=0; end;%break; end
            d2(bix==i,:) = pdist2(bpts,test(bix==i,:),'euclidean','Smallest',1);
            %D = distEucSq(bpts,test(bix==i,:)); [~,dd] = min(D,[],1); %slower
            %
            if any(d2<(tol))
                err=1; break
            end
            %
        end
        %
    end
    err = any(d2<(tol));
    %}
end

end

function [err] = bincheck(mu,pts,ix,bin)
tol = 2;
[branchdet,~,bix] = fastunique(ix'); % more bins= faster SS, slower pdist2. split bins?
%[branchdeto,~,bixo] = unique(ix','rows'); %more bins = better ss, but slower pdist2 via more runs?
%if ~all(all(bix==bixo)) || ~all(all(branchdet==branchdeto)); disp(1); end
d2 = zeros(size(ix,2),1)+tol^2;
for i=1:size(branchdet,1)
    bpts = mu{branchdet(i,1),1}{1,branchdet(i,2)};
    if ~isempty(bpts)%, err=0; end;%break; end
        d2(bix==i,:) = pdist2(bpts,pts(bix==i,:),'euclidean','Smallest',1);
        if any(d2<(tol)); err=1; break; end %early exit if collision
    end
end
err = any(d2<(tol));

end

function [ix,err] = rsplit2(mu,pts,ix,bin) %split a bin of pts across branches
% get branches from mu based on bin - will need to deal with 0/empty tree to streamline later
d = bin(1); br = bin(2); %depth and branch of the bin of test points

err=0;
if numel(mu{d,1}{3,br})==8 % only split when the bin has branches
    ix2 = all(ix==(bin'),1); %logical index of pts in the root to split
    %tmp = pts(ix2,:); %which pts being worked with
    leaves = mu{d,1}{3,br}; %indices of leaves from branch
    bcp = vertcat(mu{d+1,1}{2,leaves}); %centers of child nodes for nearest search
    [~,c2] = pdist2(bcp,pts(ix2,:),'squaredeuclidean','Smallest',1); %which branch for the bin pts
    t2 = zeros(2,numel(c2))+d+1;
    t2(2,:) = leaves(c2); ix(:,ix2) = t2;%[d+1,c2]';
end

end


function ix = rootsplit(mu,test)
n = size(test,1);
ix = ones(2,n);
roots = vertcat(mu{1,1}{2,:}); %all top-level centers
%check for multiple roots to avoid slow pdist 2 check?
%D = distEucSq(roots,test); [~,rix] = min(D,[],1); % slower under most conditions
[~,rix] = pdist2(roots,test,'squaredeuclidean','Smallest',1);
%[D] = mypdist2(roots,test); %would need to filter afterward - and too slow anyway
ix(1,:) = 1; ix(2,:) = rix;
end

function D = distEucSq(X,Y) %somewhat slower than pdist
Yt = Y';
XX = sum(X.*X,2);
YY = sum(Yt.*Yt,1);
D = bsxfun(@plus,XX,YY)-2*X*Yt;
end

function [D] = mypdist2(X,Y) %much slower than pdist2
%PDIST2 computes the squared euclidean distance between rows of X and rows of Y
%input:
% X - n*d n rows in d dimensions
% Y - m*d m rows in d dimensions
%output:
% D - n*m D_ij = norm(X(i,:)-Y(j,:))^2

[n,dx] = size(X); [m,dy] = size(Y);
assert(dx==dy,'X and Y must have same column dimension');
nx = dot(X,X,2); ny = dot(Y,Y,2);
Nx = repmat(nx,1,m); Ny = repmat(ny,1,n);
D = Nx + Ny' - 2*X*Y';
end

function [uq,ia,ic,idx1,idx2,k] = fastunique(A)
%[A_sorted, idx1] = sortrows(A); % 5.09, 7.8/5.06
%[A_sorted, idx1] = sortrows(A,2); % 3.84, 7.1/4.3, %apparently the fastest, is it as stable?
[A_sorted, idx1] = sortrows(A,[2,1]); % 4.77, 7.3,4.6
k    = find([true; any(diff(A_sorted, 1, 1), 2); true]);
idx2 = k(diff(k) >= 1);
ia  = idx1(idx2);
uq    = A(ia, :);

no = nargout(); %number of outputs declared - is there a positional version? don't need it now
if no>1
    numRows = size(A,1); %numCols = size(A,2);
    groupsSortA = A_sorted(1:numRows-1,:) ~= A_sorted(2:numRows,:);
    groupsSortA = any(groupsSortA,2);
    if (numRows ~=0)
        groupsSortA = [true; groupsSortA];       % First row is always a member of unique list.
    end
    
    ic = cumsum(full(groupsSortA));               % Lists position, starting at 1.
    ic(idx1) = ic;                          % Re-reference indC to indexing of sortA.
end
end


function it = nav(mu,bm,pt,depth,br)
if size(mu{depth,1}{1,br},2)==3
    %~isempty(mutree{depth,1}{1,br})%all(size(mutree{depth,1}{3,br})==[2,3])
    it = [depth,br]; %if corner array, nav complete
else
    %depth,br
    bdir = pt>mu{depth,1}{2,br}; %bit comparison for split axes
    [~,child] = max(all(bsxfun(@eq,bdir,bm),2),[],3); %branch index from bit comparison
    leaves = mu{depth,1}{3,br};
    bcp = vertcat(mu{depth+1,1}{2,leaves}); %centers for child nodes for pdist2
    %[~,c2] = pdist2(bcp,pt,'euclidean','Smallest',1);
    child = mu{depth,1}{3,br}(child);
    it = nav(mu,bm,pt,depth+1,child);
end
end
function te = prox(mu,pt,it,tol)
c = mu{it(1),1}{3,it(2)}; %corner array
if ~isempty(c) %if bounds present, check against them before testing against contents
    %c
    %can fail, mutate sometimes leaves points in branches after a non-exhaustive sort
    inside = pt>c(1,:)-tol & pt<c(2,:)+tol; %boolean if pt is potentially inside box
    if ~any(inside)
        te = 1; %fprintf('a,');
    else
        tpts = mu{it(1),1}{1,it(2)};
        d = vecnorm(tpts-pt,2,2); %faster for a single pt test than pdist2
        te = any(d<tol);
    end
else %if empty, can't overlap
    te = 0; %fprintf('b,');
end
end