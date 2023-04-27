function [noised, noise] = helper_noisegen(inputvol,pixang,opt)
%generate various kinds of noise based on an input volume
%is rescaling to a common standard the best way to generate similar noise spectra? or can they be derived
%independent of the initial image intensity?
%currently deprecated legacy code
arguments
    inputvol %must be rescaled before input, maybe rejigger to do rescaling and unscaling internally
    pixang
    opt.verbose = 0
    opt.diag = 0
end

%needs rework to deal with beads a bit better and have more flexibility in generating noise
%gaussian field is simple enough, just might need recalibrating
%chunky noise extrema works relatively well, might need changes
%don't have a way to reproduce very punctate noise from real tomograms. can it even be done on the tilt? why
%do tomograms often look like that?
%add sparse noise to try to simulate it?

inputvol = rescale(inputvol);

stdev = std(inputvol,0,'all'); avg = mean(inputvol,'all'); variance = var(inputvol,0,'all');
if opt.verbose==1, fprintf('Mean %g, stdev %g, var %g \n',avg,stdev,variance); end
invsize = round((10+pixang*2)/pixang); %nonlinearly scale distances to increase smoothness for low-res images

%output or have options for gaussian, uniform, chunky, and large particle noise layers?

%gaussian field noise
gauss = randn(size(inputvol))*stdev;

%chunky noise extrema bubbles
extrema = randn(size(inputvol))*stdev;
noisetop = extrema.*(extrema>stdev); noisebot = -extrema.*(extrema<-stdev); 
diltop = imdilate(noisetop,strel('sphere',invsize)); dilbot = imdilate(noisebot,strel('sphere',invsize));
extrema = diltop+dilbot*-1;

nsglobal = 0.3; nsgauss = 0.8; nsextrema = 0.6;
noise = nsglobal * (gauss*nsgauss + extrema*nsextrema);
noised = inputvol+noise;


%gaussian noise around 0, flattened according to input var
%samplenoise = randn(size(inputvol))*stdev;% + 0*variance*(randi(3,size(inputvol))-2); variance does little
%mask to only values beyond 1 stdev, both top and bottom
%noisetop = samplenoise.*(samplenoise>stdev); noisebot = -samplenoise.*(samplenoise<-stdev); 
%dilate the noise to produce intersecting spheres of noise extrema
%diltop = imdilate(noisetop,strel('sphere',invsize)); dilbot = imdilate(noisebot,strel('sphere',invsize));
%unused secondary dilation step for extra chunkyness
%diltop = imdilate(diltop,strel('sphere',invsize)); dilbot = imdilate(dilbot,strel('sphere',invsize));

%noiseout = diltop+dilbot*-1; noiseout = imgaussfilt3(noiseout,0.2);
%noise = 0.5*(samplenoise*0.8+noiseout*0.5+gaussianfield*0.8);
%noise = imgaussfilt3(samplenoise+noiseout,1);
%samplenoised = inputvol*-1 + noise;

%current settings are (generously) replicating in vitro filament scans, too much noise if anything
%part of why tomos are so gritty might be they are integer values, these are continuous

%{
if opt.diag==1 %diagnostic output figures
subplot(1,5,1); histogram(inputvol), title volume; 
subplot(1,5,2); histogram(samplenoise), title initnoise; 
subplot(1,5,3); histogram(noiseout); title 'noise extrema'; 
subplot(1,5,4); histogram(noise);  title 'applied noise';
subplot(1,5,5); histogram(samplenoised), title 'noised volume'; 
figure(); sliceViewer(noiseout); title 'noise extrema';
figure(); sliceViewer(noise); title 'applied noise';
figure(); sliceViewer(samplenoised); title 'noised sample';
end
%}

end