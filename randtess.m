function X = randtess(N,tess,sampledomain)
% randtess: samples points randomly (uniformly) from the domain or boundary of a simplicial tessellation
% usage: X = randtess(N,tess)
% usage: X = randtess(N,tess,sampledomain)
%
% arguments: (input)
%  N  - scalar, positive integer, the total number of points to sample
%
%  tess - simplicial tesselation object, may be any of a
%         triangulation object or an alpha shape, or a
%         delaunaytriangulation.
%         The code only currently runs only in 2-d or 3-d.
%
%  sampledomain - (optional) character flag, indicates where the
%         sampling will be done, thus either from the boundary
%         surface of the tessellation, or the volume itself.
%         sampledomain may be either 'volume' or 'surface', or any
%         direct contraction of either of those words. 
%
%         Default: sampledomain = 'volume'
%
% arguments: (output)
%  X  - (N,dim) array containing sampled points, where dim is the
%       dimension of the domain of the tesselation
%
%
% See also: rand, triangulation, delaunayTriangulation, alphaShape
% 
% Author: John D'Errico
% e-mail: woodchips@rochester.rr.com
% Release date: 2/26/2023

if (nargin < 3) || isempty(sampledomain)
  % the default is a volume sampling
  sampledomain = 'v';
else
  % There are only two alternatives. validatestring is robust
  % to caps or not. I'll make sure it is character though.
  sampledomain = validatestring(char(sampledomain),{'volume','surface'});
  % just grab the first letter, so 's' or 'v'.
  sampledomain = sampledomain(1);
end

%{
% insure that N is a scalar positive integer
if ~isscalar(N) || (N<=0) || (rem(N,1)~= 0)
  error('N must be a scalar positive integer')
end
%}

% tesselation must be either an alpha shape or a triangulation or a
% delaunayTriangulation. In any case, only 2-d or 3-d tesselations are
% accepted.
switch class(tess)
  case 'alphaShape'
    % extract the vertices
    domain = tess.Points;
    
    % the dimension of the domain space is
    dim = size(domain,2);
    
    % we need to extract the connectivity list for the surface, or for the
    % volume tessellation
    if strcmp(sampledomain,'v')
      % a volume sampling
      Conn = alphaTriangulation(tess);
    else
      % just from the surface
      Conn = boundaryFacets(tess);
    end

  case {'triangulation'  'delaunayTriangulation'}

    % extract the vertices
    domain = tess.Points;
    
    % the dimension of the domain space is
    dim = size(domain,2);
        
    % we need to extract the connectivity list for the surface, or for the
    % volume tessellation
    if strcmp(sampledomain,'v')
      % a volume sampling
      Conn = tess.ConnectivityList;
    else
      % just from the surface
      Conn = freeBoundary(tess);
    end

  otherwise
    error('Tess must be one of these classes: ''alphaShape'', ''triangulation'', or ''delaunayTriangulation''')

end

if dim > 3
  error('I''m sorry. this code only runs in 2-d or 3-d. It could, but I did not extend it beyond 3 dimensions')
end

% the manifold dimension would be dim-1, if this is a boundary surface sampling
% problem.
mdim = dim - strcmp(sampledomain,'s');

% compute volumes of each simplex in the complex
vols = simplexvolume(domain,Conn);


if strcmp(sampledomain,'v')
    %tvol = sum(vols) %total volume of simplexes
    tvol = volume(tess);
    N = round(tvol*N*10); %  points per unit^3
else
    tvol = surfaceArea(tess);
    N = round(tvol*N/100); % don't know how many pts per surface area this is
end
% determine which simplex each point will randomly fall in
cvol = cumsum([0;vols]); %disp(cvol') %size(cvol)
cvol = cvol/cvol(end); %disp(cvol')
q = rand(N,1);
pbins = discretize(q,cvol);

% now sample from the simplexes. First, along one edge
X = domain(Conn(pbins,1),:); %size(X)
Y = domain(Conn(pbins,2),:); %size(Y)
p = rand(N,1);
X = X.*p + Y.*(1-p); % presumes R2016b or later.

for i=2:mdim
  Y = domain(Conn(pbins,i+1),:);
  p = rand(N,1).^(1./i);
  
  X = X.*p + Y.*(1-p); % presumes R2016b or later.
end

function vols = simplexvolume(domain,Conn)
  % compute the volumes or area for each simplex
  
  switch mdim
    case 1
      % these are just edges, so find the length of each edge
      vols = domain(Conn(:,1),:) - domain(Conn(:,2),:);
      vols = sqrt(sum(vols.^2,2));
    case 2
      % these are triangles. do they live in 2 or 3 dimensions?
      % 2-d complex

      % use a 2x2 determinant, in a vectorized form
      v1 = domain(Conn(:,1),:);
      v2 = domain(Conn(:,2),:);
      v3 = domain(Conn(:,3),:);
      
      % translate to the origin
      v1 = v1-v3;
      v2 = v2-v3;
      
      if dim == 2
        % vectorized determinant
        % divide by factorial(2) for the area
        vols = (v1(:,1).*v2(:,2) - v1(:,2).*v2(:,1))/2;
      elseif dim == 3
        % get the area as the magnitude of this cross product
        vols = cross(v1,v2);
        % norm of the vector, divided by factorial(2) for the area
        vols = sqrt(sum(vols.^2,2))/2;
      end

    case 3
      % the volumes of tetrahedra
      v1 = domain(Conn(:,1),:);
      v2 = domain(Conn(:,2),:);
      v3 = domain(Conn(:,3),:);
      v4 = domain(Conn(:,4),:);
      
      % translate to the origin
      v1 = v1-v4;
      v2 = v2-v4;
      v3 = v3-v4;
      
      % vectorized determinant
      vols = v1(:,1).*v2(:,2).*v3(:,3) - ...
            v1(:,1).*v2(:,3).*v3(:,2) - ...
            v1(:,2).*v2(:,1).*v3(:,3) + ...
            v1(:,2).*v2(:,3).*v3(:,1) + ...
            v1(:,3).*v2(:,1).*v3(:,2) - ...
            v1(:,3).*v2(:,2).*v3(:,1);
    
      % divide by factorial(3) for the volume
      vols = vols/6;
      
  end

end

end

