function [dtilt,atlas,tilt] = cts_simulate_atomic(input,param,opt)
% dtilt = cts_simulate_atomic(input,param,opt)
% simulate a tilteries and reconstruct from an atomic model. more accurate and slower than vol-based method.
% ex
% dtilt = cts_simulate_atomic('gui',{'pix',9,'tilt',60:3:60},'slice',9)
%
arguments
    input = 'gui'
    param = {}
    opt.slice = 0
    opt.suffix = ''
    opt.norm = 1
end
if iscell(param), param = param_simulate(param{:}); end
if param.pix<0.5, error('inappropriate pixel size - check paramters'); end
if opt.slice==0, opt.slice = param.pix*2; end

if strcmp(input,'gui') %load model via GUI or specific filename
    [input, path] = uigetfile({'*atom.mat'},'Select input MRC or generated ts.mat',getenv('HOME')); 
    if input==0, error('At least one file must be selected or input'), end
else
    [path,input,ext] = fileparts(input); input = append(input,ext);
end
q = load(fullfile(path,input));
if ~isfield(q,'dat'), error('not a valid atomic model file'); end
mod = q.dat; split = mod.data;

fn = fieldnames(split);
atoms = zeros(0,4);
for i=1:numel(fn)
    atoms = [atoms;split.(fn{i})];
end

cd(path); %cd to the input file location to prepare session folder
%filename = append(filename,'_',opt.suffix); %generate initial filename
if ~strncmp('_',opt.suffix,1), opt.suffix = append('_',opt.suffix); end
runfolder = append('pix_',string(param.pix),'_dose_',string(sum(param.dose)),opt.suffix);
mkdir(runfolder); cd(runfolder); delete *.mrc; fprintf('Session folder: %s\n',runfolder);

file = fopen('tiltanglesT.txt','w'); fprintf(file,'%i\n',param.tilt+param.tilterr); fclose(file);
file = fopen('tiltanglesR.txt','w'); fprintf(file,'%i\n',param.tilt); fclose(file);

[vol,solv,atlas,splitvol,acount] = helper_atoms2vol(param.pix,split,mod.box);
mod.box = size(atlas)*param.pix;

%param = param_simulate('pix',3,'tilt',angles,'dose',80/4,'tiltax','Y','defocus',-2);
[tilt,dtilt,cv,cv2,ctf] = atomictiltproj(atoms,param,mod.box,opt.slice);
%sliceViewer(dtilt*-1);

WriteMRC(rescale((vol+solv)*-1),param.pix,append('0_model',opt.suffix,'.mrc'));
WriteMRC(atlas,param.pix,append('Atlas',opt.suffix,'.mrc'));

roinames = fieldnames(split); roinames = string(roinames); 
file = fopen(append('Atlas',opt.suffix,'.txt'),'w'); 
fprintf(file,'background\n'); fprintf(file,'%s\n',roinames); fclose(file);
%WriteMRC(atlas2-1,cts.param.pix,append('Atlas',opt.suffix,'.mrc'),1)

prev = append('3_tilt',opt.suffix,'.mrc');
WriteMRC(rescale(dtilt*-1),param.pix,prev);

param.size = round(mod.box/param.pix);
thick = string(round(param.size(3)*1)); %w = string(param.size(1)-50);
%reconstruct and rotate back into the proper space
%radial command for fourier filtering the output, no idea what normal runs use so random numbers
%first number radial cutoff, real tomos ~.35? cutoff slightly smoothes and increases contrast
%lower second number sharper cutoff? or fill value past cutoff?
%-hamminglikefilter should work similarly but only needs one input
%-radial default 0.35 0.035
if strcmp(param.tiltax,'Y'), tmpax = 1; else tmpax = 2; end
w = string(round(param.size(tmpax)*1));
cmd = append('tilt -tiltfile tiltanglesR.txt -RADIAL 0.35,0.035 -width ',w,...
    ' -thickness ',thick,' ',prev,' temp.mrc'); 
disp(cmd); [~] = evalc('system(cmd)'); %run the recon after displaying the command
base = append(opt.suffix,'.mrc');
%cmd = append('trimvol -rx temp.mrc ',append('5_recon_rx',base)); %#ok<NASGU>
%[~] = evalc('system(cmd)'); %run the command and capture outputs from spamming the console
cmd = append('trimvol -mode 1 -yz temp.mrc ',append('5_recon',base)); %#ok<NASGU>
[~] = evalc('system(cmd)'); %run the command and capture outputs from spamming the console
%cmd = append('clip flipz ',append('5_recon',base),' ',append('5_recon_flipz',base)); %#ok<NASGU>
%[~] = evalc('system(cmd)'); %run the command and capture outputs from spamming the console
% rx is supposed to match yz direction, but is inverted
% yz is almost matched to vol, might be ob1 in z direction.
% flipz is identical to yz
if opt.norm ==1
    [vol,head] = ReadMRC(append('5_recon',base));
    delete(append('5_recon',base));
    vol = single(vol);
    normed = (vol-mean(vol,'all'))/std(vol,1,'all');
    WriteMRC(vol,head.pixA,append('5_recon',base),1);
end
delete temp.mrc
cd(userpath)
end


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

%fprintf('CTF parameters: pixels %g angstroms, %i KeV, aberration %g nm, sigma %g, defocus %d um',...
%    param.pix, param.voltage, param.aberration, param.sigma, param.defocus)

V = param.voltage*1000; %convert from KeV to eV
cs = param.aberration/1000; %convert from mm to m
pix = param.pix/1e10; %convert from angstroms to m
Dz = param.defocus/1e6; %convert from microns to m

L = relativistic_electrons(V); %compute wavelength from voltage, correcting for relativistic speed

Ny = 1/(2*pix); %nyquist frequency for later use
B = param.sigma*Ny; %envelope factor from nyquist frequency - also incorporates the MTF signal dropoff (approx)
% make sigma param description more clear as the envelope factor
q = 0.07; %amplitude contrast 7% is approx https://www.sciencedirect.com/science/article/pii/0304399188900034
%phi = param.phase;
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
% functionalize k-space generation for cleaner code?
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

mid = round(size(input,3)/2); %middle slice for given defocus to diverge from
if numel(size(input))>2, iters = size(input,3); else iters=1; end
for i=1:iters %loop through tilts
    adj = (param.pix*slab*(i-mid))/(1e10)*1e0; %adjustment to listed defocus by depth, convert A to m
    Dzs = Dz+adj;
    [cv(:,:,i), ctf] = math_ctf(padded(:,:,i),cs,L,k,Dzs,B,q,param.phase); %get ctf-convolved subvolume
    %cv(:,:,i) = lg; %imshow(rescale(padded(:,:,i)));
    %whole tilt lowpass filter test
    %cv(:,:,i) = real(ifft2(ifftshift(fftshift( fft2(cv(:,:,i)) ).*circfilt )));
end
%sliceViewer(cv);
%imshow(rescale(lg));
convolved = cv(1+pad:end-pad,1+pad:end-pad,:); %extract image area from padded dimensions
ctf = ctf(1+pad:end-pad,1+pad:end-pad);
%fprintf('  - modulation done \n')
end

function [out,ctf] = internal_ctf(in,cs,L,k,Dz,B,A)
phi = 1*pi/3; % phase shift, ideal pi/2 phase imaging, assuming 0 otherwise
eq = pi/2*(cs*L^3*k.^4 - 2*Dz*L*k.^2); % main equation for each part of CTF signal wave
env = exp(-(k./(B)).^2); % envelope function of the overall CTF, radial signal falloff
ctf = sin(phi+eq-A).*env; % evaluate CTF terms (phase, defoc/abb, amplitude) and apply envelope falloff
out = real(ifft2(ifftshift(fftshift(fft2(in)).*ctf))); %fft stack to translate from ctf fourier to realspace
end

function L = relativistic_electrons(V) %for calculating relativistic wavelengths of electrons
m = 9.1093837e-31; %mass of electron Kg
c = 299792458;  %speed of light m/s
e = 1.60217663e-19; %charge of an electron coulombs
h = 6.62607015e-34; %planck constant m^2 Kg/s
L = h*c/sqrt(e*V*(2*m*c^2+e*V)); %calculation of wavelength L from accelerating voltage and constants
end

function [tilt,dtilt,cv,cv2,ctf] = atomictiltproj(atoms,param,boxsize,slabthick)
if param.tiltax=='Y'
    ax = [0,1,0];
else
    ax = [1,0,0];
end

eucentric = boxsize/2-[0,0,0]*0; %~25 at 12pix registers to vol but desyncs from atlas

% get the transmission wave dose 
if param.phase~=0, phi=0.8; else, phi=1; end % cut 20% of the DQE due to weird PP scattering
DQE = 0.84*0.3*phi;
d = DQE*param.dose/numel(param.tilt)*param.pix^2;
boxsize = param.pix*round(boxsize/param.pix); % adjust for weird pixel sizes

tilt = zeros(round(boxsize(1)/param.pix),round(boxsize(2)/param.pix),numel(param.tilt));
if numel(param.tilterr)~=numel(param.tilt) && param.tilterr==0
    param.tilterr = zeros(size(param.tilt));
end
prog = 0; progdel = ''; % initialize starting vals for progress bar
for t=1:numel(param.tilt)
    prog = prog + 100/numel(param.tilt); %progress update block
    progstr = sprintf('progress %3.1f  current angle: %i', prog,param.tilt(t));
    fprintf([progdel, progstr]);
    progdel = repmat(sprintf('\b'), 1, length(progstr));
    
    angle = param.tilt(t)+param.tilterr(t);
    atomtmp = atoms; %move to keeping only one set and rotating continuously? fewer reassigns
    %tmp2 = atomtmp(:,1:3); tmp2 = tmp2-eucentric; tmp2 = tmp2*rotmat(ax,deg2rad(angle)); 
    %tmp2 = tmp2+eucentric; atomtmp(:,1:3) = tmp2;
    %tmp3 = atomtmp(:,1:3); tmp3 = (tmp3-eucentric)*rotmat(ax,deg2rad(angle))+eucentric; atomtmp(:,1:3) = tmp3;
    atomtmp(:,1:3) = (atomtmp(:,1:3)-eucentric)*rotmat(ax,deg2rad(angle))+eucentric;
    
    % project a set of slices - higher resolution in Z? start with isotropy
    atomtmp(:,3) = (atomtmp(:,3)-min(atomtmp(:,3)*1,[],'all'))/slabthick*-1;
    % fixing boxsize seems to crop out excess slices
    % reformulate to add min to z instead of using an offset value?
    sz = boxsize; sz(3) = max(atomtmp(:,3),[],'all')-min(atomtmp(:,3),[],'all');
    %of = min(atomtmp(:,1:3),[],1);
    of = [0,0,min(atomtmp(:,3),[],'all')];
    [vol,solv] = helper_atoms2vol(param.pix,atomtmp,sz,of);
    vol = vol+solv;
    % get the transmission wave
    %d = param.dose*param.pix^2;
    %dvol = poissrnd((vol*1)*d,size(vol)); %extremely slow with many sections - do at the end?
    
    % propogate transmission
    mid = round(size(vol,3)/2);
    convolved = zeros(size(vol));
    cv = convolved;
    tparam = param;
    for i=1:size(vol,3)
        adj = (tparam.pix*slabthick*(i-mid))/1e4*1e0; %convert from ang to um
        tparam.defocus = param.defocus+adj;
        %tparam.defocus
        tparam.tilt = 0;
        %[convolved(:,:,i), ctf, tparam] = helper_ctf(vol(:,:,i),tparam);
        [convolved(:,:,i), ctf] = flatctf(vol(:,:,i),slabthick,tparam);
    end
    
    %[cv,ctf2] = flatctf(vol,slabthick,param);
    cv2(:,:,t) = sum(cv,3);
    
    tilt(:,:,t) = sum(convolved,3);
end
fprintf('\n%i tilts simulated\n',t)
%ctf = 0;
dtilt = poissrnd((d*rescale(tilt*1,0,1))*01,size(tilt));
cv2 = poissrnd((d*rescale(cv2*1,0,1))*01,size(cv2));
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