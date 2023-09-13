function vec = randv(row,col,ax,ang) %ang IN RADIANS
% vec = randv(row,col)
% generates an array the size of row*col, with one unitized vector per row
% advanced usage for nx3 vectors:  vec = randv(row,col,ax,ang)
% these vectors are incread constrained about the vector ax by the limits ang
% ax must be a 1x3 vector, or empty. if empty, a random axis is used for the entire set
% ang are the constraints in radians. either [min,max] or [max] that defaults to min=0
% this is the deflection angle from ax. [0,pi] for hemisphere, [80,100] for a narrow belt, etc
% the per-vector angle is generated with a uniform distribution, producing a density bias near ax and -ax

rax = randvec(row,col); %random vectors - finished or to cross with the center axis
if (isempty(ax) && isempty(ang)) || col~=3
    vec=rax; % if ax and ang are empty, return the random unconstrained vectors
else
    if isempty(ax), ax = randvec(1,3); end %if no axis given, randomize one
    if numel(ang)==1, ang(2)=ang(1); ang(1)=0; end %if only 1 ang, use as max against min 0
    %ang(2) = ang(2)-ang(1); %store difference from min for simpler following code
    %nrep = row/size(ax,1); %number of replicates needed to match matrix size for cross
    ax = ax/norm(ax); %unitize target vector to avoid miscalculation
    rotax = cross(repmat(ax,row/size(ax,1),1),rax); %compute orthogonal axes to rotate
    rotax = (rotax'./vecnorm(rotax'))'; %unitize orthogonal axes
    vec = zeros(row,col);
    for i=1:row %loop because matrix multiplication can't be vectorized?
        R = rotmat(rotax(i,:),rand*(ang(2)-ang(1))+ang(1)); %rot matrix about axis, random angle within range
        vec(i,:) = ax*R; %vector of rotated point away form cone center
    end
end
end

function [vec] = randvec(row,col)
vec = randn(col,row); %random normal numbers for evenly-distributed vector directions
vec = (vec./vecnorm(vec))'; %unitize vectors to length 1 for sphere vectors
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