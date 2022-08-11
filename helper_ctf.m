function [convolved, ctf, param] = helper_ctf(input,param,pad)
%no help yet, good luck!
%check cts_param for what the param argument should be, this program defaults to it
arguments
    input
    param = {}
    pad = 10; %padding added to volume before any computations
end
if iscell(param), param = cts_param(param{:}); end
if param.ctfoverlap==0, convolved=input; return; end %if overlap==0, skip doing CTF
tilts = param.tilt(1):param.tilt(2):param.tilt(3);

%{
if strcmp(param,'default'), param = struct('volt',3e2,'sphereabb',2.7,'defocus',-6,'sigma',1);
elseif strcmp(param,'gui')
    current = pwd; cd(userpath); cd('Toolbox/tomo_sim/defaults');
    [file, path] = uigetfile('*.mat','Select the desired parameter file');
    load(fullfile(path,file),'param')
    cd(current)
elseif strcmp(param,'manual')
    param = [];
    p = inputdlg( {'Microscope voltage (kV): ','Microscope spherical aberration (mm): ',...
        'Session focus (negative for defocus, nm): ','Sigma value for CTF: ',...
        'Save these parameters as a file for later use? (1 for yes): '},'Microscope and imaging parameters');
    p = str2double(p); %convert to numbers
    param.volt = p(1); param.sphereabb = p(2); param.defocus = p(3); param.sigma = p(4);
    if p(5)==1
        current = pwd;
        q = append('volt',string(param.volt),'_sa',string(param.sphereabb),...
            '_defoc',string(param.defocus),'_sigma',string(param.sigma),'.mat');
        cd(userpath); cd('Toolbox/tomo_sim/defaults');
        save(q,'param'); cd(current); disp('Saved input parameters')
    end
end
%}

V = param.voltage*1000; %convert from KeV to eV
cs = param.aberration/1000; %convert from mm to m
pix = param.pix/1e10; %convert from angstroms to m
Dz = param.defocus/1e6; %convert from microns to m

L = relativistic_electrons(V); %compute wavelength from voltage, correcting for relativistic effects

Ny = 1/(2*pix); B = param.sigma*Ny; q = 0.07; %nyquist, envelope, and amplitude contrast values
%envelope still needs validation and corroboration to our real data

fprintf('Parameters: pixels %g angstroms, voltage %i KeV, aberration %g nm, sigma %g, defocus %d nm\n',...
    param.pix, param.voltage, param.aberration, param.sigma, param.defocus)

k = 1:size(input,1); divs = k(rem(size(input,1),k)==0); %find divisible factors from volume size
binspacing = divs(round(end/2)); %use the middle divisor as the spacing
bins = size(input,1)/binspacing+1; %determine bins from the spacing and vol size

binlength = binspacing*param.ctfoverlap; %this is one-sided, not full length
edge = binlength+binspacing*(param.ctfoverlap-1);
padded = padarray(input,[edge pad],mean(input,'all'));
bincenter = (linspace(binlength,size(padded,1)-binlength,bins+(param.ctfoverlap-1)*2)); %compute strip centres

xl = size(padded,2); %uses dim2 to avoid needing to permute the CTF to the image space
yl = binlength*2;
%zl = size(padded,3);%+pad*2; %not using due to 2d implementation
[x,y] = meshgrid(-Ny:2*Ny/xl:Ny-Ny/xl,-Ny:2*Ny/yl:Ny-Ny/yl);%,-Ny:2*Ny/zl:Ny-Ny/zl);

k = sqrt(x.^2+y.^2);%+z.^2); %evaluate inverse distance, identical for all strips
cv = zeros(size(padded)); %pre-initialize output array

%weight for overlapping portions of bins
inc = 0.5/binlength; %make edges not fall to 0/centre not 1
weight = 1-abs(linspace(-1+inc,1-inc,yl))'; %more simple weight function for strips
weight = repmat(weight/param.ctfoverlap,[1 xl]); %replicate across the strip length
%[weight] = abs(abs(meshgrid(-1+0.5/binlength:1/(binlength):1-0.5/binlength,1:size(padded,2)))-1);
%weight = permute(weight,[2 1])/overlap; %former method, a bit more convoluted to get the right values

for i=1:numel(tilts) %loop through tilts
    shift = tand(tilts(i)); %proportion of length by tilt to compute vertical displacement
    
    for j=1:numel(bincenter)
        six = round([1+bincenter(j)-binlength, bincenter(j)+binlength]);
        sdist = pix*(bincenter(j)-size(padded,1)/2)*1; %calculate horizontal distance
        Dzs = Dz + shift*sdist; %average defocus in the strip given tilt angle and horizontal distance
        in = padded(six(1):six(2),1:end,i);
        [lg, ctf] = internal_ctf(in,cs,L,k,Dzs,B,q); %get ctf-convolved subvolume
        lg = lg.*weight; %scale by weight for overlap of nearby strips
        cv(six(1):six(2),1:end,i) = cv(six(1):six(2),1:end,i)+lg;
        %wgt(six(1):six(2),1:end,i) = wgt(six(1):six(2),1:end,i)+weight;
        
        %verbose real defoc listing
        %fprintf('tiltangle %g ix %g to %g strip distance %g at defoc %g\n',tilts(i),six(1),six(2),sdist,Dzs)
    end
    
end
convolved = cv(1+edge:end-edge,1+pad:end-pad,:); %extract image area from padded dimensions
%wgt = wgt(1+edge:end-edge,1+pad:end-pad,:);

end

function [out,ctf] = internal_ctf(in,cs,L,k,Dz,B,q)
eq = pi/2*(cs*L^3*k.^4 - 2*Dz*L*k.^2); %main equation for each part of CTF
env = exp(-(k./(B)).^2); %envelope function of the overall CTF
ctf = ( sqrt(1-q^2)*sin(eq) - q*cos(eq) ) .*env; %complete CTF evaluation
out = real(ifft2(ifftshift(fftshift(fft2(in)).*ctf))); %fft stack to translate from ctf fourier to realspace
end

function L = relativistic_electrons(V) %for calculating relativistic wavelengths of electrons
m = 9.1093837e-31; %mass of electron Kg
c = 299792458;  %speed of light m/s
e = 1.60217663e-19; %charge of an electron coulombs
h = 6.62607015e-34; %planck constant m^2 Kg/s
L = h*c/sqrt(e*V*(2*m*c^2+e*V)); %calculation of wavelength L from voltage and derived constants
end