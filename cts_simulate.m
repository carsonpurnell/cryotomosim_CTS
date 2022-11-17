function [noised, conv, tiltseries, atlas, ctf] = cts_simulate(sampleMRC,param,opt)
%[noised, convolved, tiltseries] = tomosim_simulate(sampleMRC, param, opt)
%simulates tomographic data as if collected from a sample and reconstructs a tomogram
%
%  Inputs
%
%sample MRC     required, preferred input = 'gui'
%an .mrc volume or .mat of a cts struct that provides the simulated sample (the .mat is much better)
%'gui' option allows choosing the file with matlab's basic file browser utility
%an input .mrc must be positive-scale density
%
%param          default = cts_param == {}
%see cts_param for arguments and their usage. to pass arguments to cts_param, {enclose in brackets}
%controls most of the behavior of the simulation through name-value pairs. use {'gui'} for manual input.
%
%  name-val options:
%
%suffix         default '' (empty)
%suffix appended to the output filenames
%
%atlasindividual    default 0
%if 1, in addition to generating a single-volume label atlas a binary volume is generated for each particle
%
%dynamotable         default 0
%if 1, generates dynamo .tbl files for each target particle of the simulation model
%
%
%  Outputs
%
%outputs some workspace variables of intermediates if you want them for something
%makes a folder with several intermediate volumes and a tiltangles.txt
%0_model   the initial model from the input, intensify flipped to EM negative density
%1_tilt    tiltprojection of initial model
%3_ctf     CTF convolved tiltseries
%4_dose    additional noise added after CTF
%5_recon   reconstruction of a tomogram
%if the input was a ts .mat from tomosim_model, also outputs individual models for each constituent class
%each unique target ID counts as a class, as well as beads (if any)

arguments
    sampleMRC (1,1) string %full path to input mrc/ts, or 'gui' for browser
    param = {} %input a cts_param call, or within {} to send to it
    opt.suffix string = ''
    opt.atlasindividual = 0
    opt.dynamotable = 0
end
if iscell(param), param = cts_param(param{:}); end %parse params if given as argument input

if strcmp(sampleMRC,'gui') %load model via GUI or specific filename
    [sampleMRC, path] = uigetfile({'*.mrc;*.mat'},'Select input MRC or generated ts.mat',getenv('HOME')); 
    if sampleMRC==0, error('At least one file must be selected or input'), end
else
    [path,sampleMRC,ext] = fileparts(sampleMRC); sampleMRC = append(sampleMRC,ext);
end

[path, filename, ext] = fileparts(fullfile(path,sampleMRC));
switch ext
    case '.mat'
        q = load(fullfile(path,sampleMRC));
        if ~isfield(q,'cts'), error('Selected mat file is not a tomosim structure'); end
        cts = q.cts; vol = cts.vol; pixelsize = cts.pix(1);
    case '.mrc'
        [vol, head] = ReadMRC(fullfile(path,sampleMRC)); 
        pixelsize = head.pixA; cts = 0;
    otherwise
        error('selected file is not a .mat or a .mrc, aborting')
end

cd(path); %cd to the input file location to prepare session folder
filename = append(filename,'_',opt.suffix); %generate initial filename
runfolder = append('sim_dose_',string(sum(param.dose)),'_',opt.suffix);
mkdir(runfolder); cd(runfolder); delete *.mrc; fprintf('Session folder: %s\n',runfolder);

if param.pix==0, param.pix=pixelsize; end %override pixel size unless 0
param.size = size(vol); %fix for x axis tilt size and for detect thickness
if strcmp(param.tiltax,'X')
    vol = imrotate3(vol,90,[0 1 0]); %rotate X tiltax volume to be not insane (still mirrored and tilted)
else
    param.tiltax = 'Y';
end

%run the simulation itself within the subfunction. might extend 'real' to also 'ideal' later
[noised, conv, tiltseries, ctf] = internal_sim(vol,filename,param,'real');

if isstruct(cts) %if a tomosim formatted .mat struct is selected, generate a particle atlas
    atlas = helper_particleatlas(cts,opt.atlasindividual,opt.dynamotable);
end

cd(userpath) %return to the user directory
end

function [noised, convolved, tilt, ctf] = internal_sim(in,filename,param,type)
pix = param.pix;

base = append(filename,'.mrc'); 
%do inversion differently to get better dose numbers?
%might use matlab imcomplement for simplicity
in = rescale(in*-1,min(in,[],'all'),max(in,[],'all')); %rescale to same range to avoid 0 and clamping
prev = append('0_model_',base);
WriteMRC(in,pix,prev) %write as initial model and for tiltprojection

donoise = 0; convolved = 0; noised = 0;
if strcmp(type,'real')
%future tilt randomization here?
file = fopen('tiltangles.txt','w'); fprintf(file,'%i\n',param.tilt); fclose(file);

if donoise==1
samplenoised = helper_noisegen(in,pixelsize); %add multifactor noise
prev = append('1_noised_',base); WriteMRC(samplenoised,pix,prev);
end
end

%project the tiltseries
tbase = append('1_tilt_',base);
w = string(round(param.size(1)*0.5)); %need to make width useful for avoiding empty ends of tilt
%better to use width or x/y min max extents?
%how to compute the width to use to avoid edge loss? also requires cropping the atlas to match
cmd = append('xyzproj -axis ', param.tiltax, ' -width ',w,' -tiltfile tiltangles.txt ',prev,' ',tbase);
%-ray borks up smaller width completely, -constant makes beads slightly weird
prev = tbase; disp(cmd); [~] = evalc('system(cmd)'); %run command, capture console spam

if strcmp(type,'real') %electron detection and CTF
tilt = ReadMRC(prev); %load the projected tiltseries as a volume

order = 2; %electron detection changable order thing because i still don't know which is better!
if order==1 %dose first, hackjob dose increase to get the scaling to work
    detected = helper_electrondetect(tilt,param);
    WriteMRC(detected,pix,append('2_dosetilt_',base));
    [convolved,ctf] = helper_ctf(detected,param); %per-tilt ctf convolution
    prev = append('3_ctf_',base);
    WriteMRC(convolved,pix,prev); %save the convolved image for review
end
if order==2 %dose second, reduces pixel size effect so more generalized but values are arbitrary?
    [convolved,ctf] = helper_ctf(tilt,param); %per-tilt ctf convolution
    prev = append('3_ctf_',base);
    WriteMRC(convolved,pix,prev); %save the convolved image for review
    convolved = rescale(convolved,min(tilt,[],'all'),max(tilt,[],'all')); %fix negative CTF
    detected = helper_electrondetect(convolved,param);
    prev = append('4_dose_',base);
    WriteMRC(detected,pix,prev);
end

if donoise==1 %hard coded multifactorial noise toggle
noised = helper_noisegen(convolved,pix); %add multifactor noise
WriteMRC(noised,pix,append('4_noised_',base)); %save the noised volume for reconstruction
prev = append('4_noised_',base);
end
end

thick = string(round(param.size(3)*1)); w = string(param.size(1)-200);
%reconstruct and rotate back into the proper space
%radial command for fourier filtering the output, no idea what normal runs use so random numbers
%
%lower second number sharper cutoff? or fill value past cutoff?
%-hamminglikefilter should work similarly but only needs one input
cmd = append('tilt -tiltfile tiltangles.txt -RADIAL 0.35,0.05 -width ',w,' -thickness ',thick,' ',prev,' temp.mrc'); 
disp(cmd); [~] = evalc('system(cmd)'); %run the recon after displaying the command
cmd = append('trimvol -rx temp.mrc ',append('5_recon_',base)); %#ok<NASGU>
[~] = evalc('system(cmd)'); %run the command and capture outputs from spamming the console

delete temp.mrc %remove temporary files after they are used for rotation
end