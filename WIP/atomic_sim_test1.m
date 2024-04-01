
% collect atomic model into single set of points
fn = fieldnames(split);
cen = boxsize/2;
atoms = zeros(0,4);
for i=1:numel(fn)
    atoms = [atoms;split.(fn{i})];
end

%%
angles = -60:10:60;
param = param_simulate('pix',8,'tilt',angles,'dose',50);
[tilt,dtilt,cv,cv2,ctf] = atomictiltproj(atoms,param,angles,boxsize,20);
sliceViewer(cv);

%% 
%{
% rotate the model to an angle (eucentric adjustment? rotate about 0?)
angle = 60;
ax = [1,0,0];
atoms(:,1:3) = (atoms(:,1:3)-cen)*rotmat(ax,deg2rad(angle))+cen;

% project a set of slices - higher resolution in Z? start with isotropy
pix = 8;
slabthick = 2;
atoms(:,3) = (atoms(:,3)-min(atoms(:,3),[],'all'))/slabthick;
sz = boxsize; sz(3) = max(atoms(:,3),[],'all');
vol = helper_atoms2vol(pix,atoms,sz);

%sim params
param = param_simulate('pix',pix,'tilt',zeros(1,size(vol,3)));
param.tilt = -40:numel(param.tilt)-41;

% get the transmission wave
d = param.dose*pix^2;
dvol = poissrnd(rescale(vol*-1)*d,size(vol));

% propogate transmission
mid = round(size(dvol,3)/2);
convolved = zeros(size(dvol));
for i=1:size(dvol,3)
    adj = (pix*slabthick*(i-mid))/1e4*3e1;
    param.defocus = -5+adj;
    param.tilt = 0;
    [convolved(:,:,i), ctf, param] = helper_ctf(dvol(:,:,i),param);
end
proj = rescale(sum(convolved,3));

imshow(proj);
%}
%% internal functs

function [convolved,ctf] = flatctf(input,slab,param,pad)

arguments
    input
    slab
    param = {} %direct inputs from cts_param
    pad = 10; %padding added to volume before any computations
end
if iscell(param), param = cts_param(param{:}); end %is this needed anymore?
if param.ctfoverlap==0, convolved=input; ctf=0; return; end %if overlap==0, skip doing CTF

fprintf('CTF parameters: pixels %g angstroms, %i KeV, aberration %g nm, sigma %g, defocus %d um',...
    param.pix, param.voltage, param.aberration, param.sigma, param.defocus)

V = param.voltage*1000; %convert from KeV to eV
cs = param.aberration/1000; %convert from mm to m
pix = param.pix/1e10; %convert from angstroms to m
Dz = param.defocus/1e6; %convert from microns to m

L = relativistic_electrons(V); %compute wavelength from voltage, correcting for relativistic speed

Ny = 1/(2*pix); %nyquist frequency calculation - functionalize?
B = param.sigma*Ny; %envelope factor from nyquist frequency - also incorporates the MTF signal dropoff (approx)
q = 0.07; %amplitude contrast value - 7% is generalization
%envelope/amplitude still needs validation and corroboration to our real data

% crunchy strip math thing - replace with subfunction, extend from 2d to 3d?
%k = 1:size(input,1); %divs = k(rem(size(input,1),k)==0); %find divisible factors from volume size

%binspacing = divs(round(end/2)); %use the middle divisor as the spacing
%bins = size(input,1)/binspacing+1; %determine bins from the spacing and vol size

%binlength = binspacing*param.ctfoverlap; %this is one-sided, not full length
%edge = binlength+binspacing*(param.ctfoverlap-1);
padded = padarray(input,[pad pad],mean(input,'all'));
%bincenter = (linspace(binlength,size(padded,1)-binlength,bins+(param.ctfoverlap-1)*2)); %compute strip centres
% crunchy strip math thing

xl = size(padded,2); %uses dim2 to avoid needing to permute the CTF to the image space
yl = size(padded,1);%binlength*2;
%zl = size(padded,3);%+pad*2; %not using due to 2d implementation
[x,y] = meshgrid(-Ny:2*Ny/xl:Ny-Ny/xl,-Ny:2*Ny/yl:Ny-Ny/yl);%,-Ny:2*Ny/zl:Ny-Ny/zl);
k = sqrt(x.^2+y.^2);%+z.^2); %evaluate inverse distance, identical for all strips
%imshow(rescale(k))
%{
whole-tilt envelope setup stuff - deprecated
xf = size(padded,2); yf = size(padded,1);
[w,u] = meshgrid(-Ny:2*Ny/xf:Ny-Ny/xf,-Ny:2*Ny/yf:Ny-Ny/yf);
kf = sqrt(w.^2+u.^2);
[r,c] = meshgrid(1:xf,1:yf);
circfilt = sqrt((r-xf/2-0.5).^2+(c-yf/2-0.5).^2)<50;
%}

cv = zeros(size(padded)); %pre-initialize output array

%generate weights for overlapping portions of bins
%inc = 0.5/binlength; %make edges not fall to 0/centre not 1
%weight = 1-abs(linspace(-1+inc,1-inc,yl))'; %more simple weight function for strips
%weight = repmat(weight/param.ctfoverlap,[1 xl]); %replicate across the strip length
%[weight] = abs(abs(meshgrid(-1+0.5/binlength:1/(binlength):1-0.5/binlength,1:size(padded,2)))-1);
%weight = permute(weight,[2 1])/overlap; %former method, a bit more convoluted to get the right values

mid = round(size(input,3)/2);
for i=1:size(input,3) %loop through tilts
    adj = (param.pix*slab*(i-mid))/(1e10)*1e0;
    %param.defocus = -5+adj;
    %shift = tand(param.tilt(i)); %proportion of length by tilt to compute vertical displacement
    %for j=1:numel(bincenter)
        %six = round([1+bincenter(j)-binlength, bincenter(j)+binlength]);
        %sdist = pix*(bincenter(j)-size(padded,1)/2)*1; %calculate horizontal distance (dev multiplier)
        %Dzs = Dz + shift*sdist; %average defocus in the strip given tilt angle and horizontal distance
        %in = padded(six(1):six(2),1:end,i); %input slice for convolution
        Dzs = Dz+adj;
        [lg, ctf] = internal_ctf(padded(:,:,i),cs,L,k,Dzs,B,q); %get ctf-convolved subvolume
        %lg = lg.*weight; %scale by weight for gradient overlap of strips
        %cv(six(1):six(2),1:end,i) = cv(six(1):six(2),1:end,i)+lg;
        cv(:,:,i) = lg; %imshow(rescale(padded(:,:,i)));
        %verbose real defoc listing
        %fprintf('tiltangle %g ix %g to %g strip distance %g at defoc %g\n',...
            %param.tilt(i),six(1),six(2),sdist,Dzs)
    %end
    %whole tilt lowpass filter test
    %cv(:,:,i) = real(ifft2(ifftshift(fftshift( fft2(cv(:,:,i)) ).*circfilt )));
end
%sliceViewer(cv);
%imshow(rescale(lg));
convolved = cv(1+pad:end-pad,1+pad:end-pad,:); %extract image area from padded dimensions
ctf = ctf(1+pad:end-pad,1+pad:end-pad);
fprintf('  - modulation done \n')
end

function [out,ctf] = internal_ctf(in,cs,L,k,Dz,B,q)
eq = pi/2*(cs*L^3*k.^4 - 2*Dz*L*k.^2); %main equation for each part of CTF signal wave
env = exp(-(k./(B)).^2); %envelope function of the overall CTF, radial signal falloff
ctf = ( (1-q)*sin(eq) + (1)*q*cos(eq) ) .*env; %evaluate phase and amp components, amplitude reduces halo
out = real(ifft2(ifftshift(fftshift(fft2(in)).*ctf))); %fft stack to translate from ctf fourier to realspace
end

function L = relativistic_electrons(V) %for calculating relativistic wavelengths of electrons
m = 9.1093837e-31; %mass of electron Kg
c = 299792458;  %speed of light m/s
e = 1.60217663e-19; %charge of an electron coulombs
h = 6.62607015e-34; %planck constant m^2 Kg/s
L = h*c/sqrt(e*V*(2*m*c^2+e*V)); %calculation of wavelength L from accelerating voltage and constants
end

function [tilt,dtilt,c2,cv2,ctf2] = atomictiltproj(atoms,param,angles,boxsize,slabthick)
ax = [1,0,0];
cen = boxsize/2;
%angles = param.tilt;
% get the transmission wave
DQE = 0.84;
d = DQE*param.dose/numel(param.tilt)*param.pix^2;
boxsize = param.pix*round(boxsize/param.pix);

tilt = zeros(boxsize(1)/param.pix,boxsize(2)/param.pix,numel(param.tilt));
for t=1:numel(angles)
    disp(angles(t))
    angle = angles(t);
    atomtmp = atoms;
    atomtmp(:,1:3) = (atomtmp(:,1:3)-cen)*rotmat(ax,deg2rad(angle))+cen;
    
    % project a set of slices - higher resolution in Z? start with isotropy
    %pix = 8;
    %slabthick = 10;
    atomtmp(:,3) = (atomtmp(:,3)-min(atomtmp(:,3),[],'all'))/slabthick;
    % fixing boxsize seems to crop out excess slices
    sz = boxsize; sz(3) = max(atomtmp(:,3),[],'all')-min(atomtmp(:,3),[],'all');
    %of = min(atomtmp(:,1:3),[],1);
    of = [0,0,min(atomtmp(:,3),[],'all')];
    %size(of)
    [vol,solv] = helper_atoms2vol(param.pix,atomtmp,sz,of);
    vol = vol+solv;
    %sim params
    %param = param_simulate('pix',param.pix,'tilt',zeros(1,size(vol,3)));
    %param.tilt = -40:numel(param.tilt)-41;
    
    % get the transmission wave
    %d = param.dose*param.pix^2;
    %dvol = poissrnd((vol*1)*d,size(vol)); %extremely slow with many sections - do at the end?
    
    % propogate transmission
    mid = round(size(vol,3)/2);
    convolved = zeros(size(vol));
    tparam = param;
    for i=1:size(vol,3)
        adj = (tparam.pix*slabthick*(i-mid))/1e4*1e2;
        tparam.defocus = param.defocus+adj;
        tparam.tilt = 0;
        %[convolved(:,:,i), ctf, tparam] = helper_ctf(vol(:,:,i),tparam);
    end
    
    [cv2,ctf2] = flatctf(vol,slabthick,param);
    c2(:,:,t) = sum(cv2,3);
    
    tilt(:,:,t) = sum(convolved,3);
end
ctf = 0;
dtilt = poissrnd((d*rescale(tilt*1,0,1))*01,size(tilt));
c2 = poissrnd((d*rescale(c2*1,0,1))*01,size(c2));
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