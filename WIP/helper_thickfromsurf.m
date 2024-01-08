function [tilts,gridt] = helper_thickfromsurf(surfaces,box,pix,angles)

[tilts,gridt] = thicktilts(surfaces,box,pix,angles);

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