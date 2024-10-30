function vec = randv(row,col,ax,ang) 
% vec = randv(row,col,ax,ang) % ang is an angle in RADIANS
% limitation: only working for 1x3 vectors (relies on exclusively 1x3 cross function)
% generates row vectors of length col on a unit sphere, under constraints in ax/ang (optional)
% if ax/ang are input, vectors are constrained along those axes up to ang offset angle
% vectors are NOT uniformly spaced, but uniformly deflected from the axis (polar density bias)
% ang can be 1x2 for both minimum and maximum deflection (order invariant)
% if ang is input but ax is null, a random axis is used
% ax can accept a list of axes, replicating them to row*col size for ouput (may change, it was a weird choice)
%
% ax must be a 1x3 vector, or empty. if empty, a random axis is used for the entire set
% ang are the constraints in radians. either [min,max] or [max] that defaults to min=0
% this is the deflection angle from ax. [0,pi] for hemisphere, [0.4,0.5] for a narrow belt, etc

arguments
    row
    col
    ax = []
    ang = []
end

%rax = randvec(row,col); %random vectors - finished or to cross with the center axis
rax = randn(col,row); rax = (rax./vecnorm(rax))'; %random initial vectors
if (isempty(ax) && isempty(ang)) || col~=3
    vec=rax; % if ax and ang are empty, return the random unconstrained vectors
else
    if isempty(ax), ax = randn(col,1); ax = (ax./vecnorm(ax))'; end %if no axis given, randomize one
    if numel(ang)==1, ang(2)=ang(1); ang(1)=0; end %if only 1 ang, use as max against min 0
    %nrep = row/size(ax,1); %number of replicates needed to match matrix size for cross
    %might be able to rework ax to be a column of axes. need to be vecnormed.
    
    % now distributes outputs equally among a list of input axes
    ax = (ax'./vecnorm(ax'))'; % unitize target vector
    nrep = ceil(row/size(ax,1)); % needed number of axes for outputs
    axrep = repmat(ax,nrep,1); axrep = axrep(1:row,1:col); % replicate input axes to needed length
    rotax = cross(axrep,rax,2); % compute orthogonal axes to rotate
    rotax = (rotax'./vecnorm(rotax'))'; % unitize orthogonal axes
    vec = zeros(row,col);
    
    for i=1:row %loop because matrix multiplication can't be vectorized?
        R = rotmat(rotax(i,:),rand*(ang(2)-ang(1))+ang(1)); %rot matrix about axis, random angle within range
        vec(i,:) = axrep(i,1:3)*R; %vector of rotated point away form cone center
    end
end
end

function [vec] = randvec(row,col)
vec = randn(col,row); vec = (vec./vecnorm(vec))'; % random uniform vectors on unit sphere
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