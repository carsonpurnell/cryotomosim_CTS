function [tilts,gridt] = helper_thickfromsurf(surfaces,box,pix,angles,ax)
%[tilts,gridt] = helper_thickfromsurf(surfaces,box,pix,angles,axis)

tilts = zeros(box(1)/pix,box(2)/pix,numel(angles));
gridt(1:2) = {zeros(box(1)/pix,box(2)/pix,numel(angles))};
for i=1:numel(angles)
    sampsurfrot = rotsurf(surfaces,box,deg2rad(angles(i)),ax);
    gridt{1}(:,:,i) = regrid(sampsurfrot{1},box,pix);
    %if ~all(isfinite(gridt{1}(:,:,i)),'all'), fprintf('top infinite at angle %i \n',angles(i)); end
    gridt{2}(:,:,i) = regrid(sampsurfrot{2},box,pix);
    %if ~all(isfinite(gridt{2}(:,:,i)),'all'), fprintf('bot infinite at angle %i \n',angles(i)); end
    %thick = gridt{1}(:,:,i)-gridt{2}(:,:,i); 
    tilts(:,:,i) = gridt{1}(:,:,i)-gridt{2}(:,:,i);
end
%[a,b] = bounds(surfaces{1})
%[a,b] = bounds(surfaces{2})
%plot3p([surfaces{1};surfaces{2}],'.'); axis equal
%plot3p([sampsurfrot{1};sampsurfrot{2}],'.'); axis equal
%histogram(surfaces{1}(:,1)); hold on
%histogram(surfaces{2}(:,1));
%histogram(sampsurfrot{1}(:,3)); hold on
%histogram(sampsurfrot{2}(:,3));
end

function surfcell = rotsurf(surfcell,box,theta,ax)
spec = find(ax); spec = 1+rem(spec,2);
cen = zeros(1,3); cen(spec) = box(spec)/2*1; %cen(3) = 0*box(3)/2;
%cen = [box(1)/2,0,0] %not appropriate for X tilt
R = rotmat(ax,theta); %R2 = makehgtform('xrotate',theta); R3 = makehgtform('yrotate',theta);
for i=1:2
    surfcell{i} = (surfcell{i}-cen)*R+cen;
    %if ~all(isfinite(surfcell{i}),'all'),fprintf('infinite at angle %i',rad2deg(theta)); end
end

end

function pts = regrid(pts,box,pix)
%tmp = [round(pts(:,1:2)/pix),pts(:,3)]; %tmp = prune(tmp,box/pix); 
%prune first because accum can't handle out-of-bounds
p = round(pts(:,1:2)/pix); v = pts(:,3);
[p,v] = prune2(p,v,box/pix);
%tmp = sortrows(tmp); %takes more time than it saves from accumarray
%pts = accumarray(tmp(:,1:2),tmp(:,3),box(1:2)/pix,@mean);%,mean(tmp(:,3))); %mean 300x slower
pts = accumarray(p,v,box(1:2)/pix)./accumarray(p,1,box(1:2)/pix);
end

function [p,v] = prune2(p,v,box)
p = real(p);
for i=1:2
    ix = p(:,i) <= box(i) & p(:,i) >= 1; %index points inside the box
    p = p(ix,:); v = v(ix); %drop points outside the box
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