function [out,ctf] = math_ctf(in,cs,L,k,Dz,B,A)
%A = -0.07; % scatter amplitude, approx .07 (replace q val?)
phi = pi/2; % phase shift, ideal pi/2 phase imaging, assuming 0 otherwise
eq = pi/2*(cs*L^3*k.^4 - 2*Dz*L*k.^2); % main equation for each part of CTF signal wave
env = exp(-(k./(B)).^2); % envelope function of the overall CTF, radial signal falloff
ctf = sin(phi+eq-A).*env; % combined CTF equation terms and applied envelope
out = real(ifft2(ifftshift(fftshift(fft2(in)).*ctf))); %fft stack to translate from ctf fourier to realspace
end