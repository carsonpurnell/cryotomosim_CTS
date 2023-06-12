%% wavy surface generator
sz = [400,300];
n = 2.5; % noise magnitude

[X,Y] = ndgrid(1:sz(1),1:sz(2));
i = min(X-1,sz(1)-X+1); j = min(Y-1,sz(2)-Y+1);
H = exp(-.5*(i.^2+j.^2)/n^2);
Z = real(ifft2(H.*fft2(randn(sz)))); % 0-centered, approximately normal

pts = [X(:),Y(:),Z(:)*1e3];
n = size(pts,1); ix = randperm(n); ix = ix(1:round(n/10));
pts2 = pts(ix,:);

%figure();
%plot3(pts2(:,1),pts2(:,2),pts2(:,3),'.'); axis equal;