function [out,ctf] = math_ctf(in,cs,L,k,Dz,B,A,phi)
eq = pi/2*(cs*L^3*k.^4 - 2*Dz*L*k.^2); % main equation for each part of CTF signal wave
env = exp(-(k./(B)).^2); % envelope function of the overall CTF, radial signal falloff
ctf = sin(phi+eq-A).*env; % evaluate CTF terms (phase, defoc/abb, amplitude) and apply envelope falloff
out = real(ifft2(ifftshift(fftshift(fft2(in)).*ctf))); %fft stack to translate from ctf fourier to realspace
end