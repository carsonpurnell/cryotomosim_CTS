function [noised, conv, tiltseries] = cts_simulate(sampleMRC,param,opt)
%{ help block
%[noised, convolved, tiltseries] = tomosim_simulate(sampleMRC, param, opt)
%simulates tomographic data as if collected from a sample and reconstructs a tomogram
%
%  Inputs
%
%sample MRC     required, preferred input = 'gui'
%an .mrc volume or .mat of a ts struct that provides the simulated sample (the .mat is much better)
%'gui' option allows choosing the file with matlab's basic file browser utility
%an input .mrc must be positive-scale density
%
%param          default = cts_param == {}
%see cts_param for arguments and their usage. to pass arguments to cts_param, {enclose in brackets}
%controls most of the behavior of the simulation through name-value pairs. use {'gui'} for manual input.
%
%opt.suffix         default ''
%suffix appended to the output filenames
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
%}
arguments
    sampleMRC (1,1) string %full path to input mrc/ts, more realistically 'gui' for browser
    param = {} %input a cts_param call, or within {} to send to it
    opt.suffix string = ''
    opt.bin = 1 %by default binarize individual model outputs
end
if iscell(param), param = cts_param(param{:}); end

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
    if ~isfield(q,'ts'), error('Selected mat file is not a tomosim structure'); end
    ts = q.ts; vol = ts.vol; pixelsize = ts.pix(1);
case '.mrc'
    [vol, head] = ReadMRC(fullfile(path,sampleMRC)); 
    pixelsize = head.pixA; ts = 0;
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
[noised, conv, tiltseries] = internal_sim(vol,filename,param,'real');

if isstruct(ts) %if a tomosim formatted .mat struct is selected, generate individual particle standards
%if isfield(ts,'splitmodel') %obsolete second check, an invalid ts would already break
    helper_particleatlas(ts);%,'dynamotable',1,'individual',1);
    %{
    splitnames = fieldnames(ts.splitmodel);
    for i=1:numel(splitnames)
        filename = append('ind',string(i),'_',splitnames{i},'.mrc');
        ind = ts.splitmodel.(splitnames{i});
        if opt.bin==1; ind=imbinarize(rescale(ind)); end
        WriteMRC(ind,param.pix,filename)
    end
    %}
%end
%if isfield(ts.model,'beads'), WriteMRC(ts.model.beads,param.pix,'ind_beads.mrc'), end
%if isfield(ts.model,'mem'), WriteMRC(ts.model.mem,param.pix,'ind_membrane.mrc'), end
end

cd(userpath) %return to the user directory
end

function [noised, convolved, tilt] = internal_sim(in,filename,param,type)
pix = param.pix;

base = append(filename,'.mrc'); 
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
w = string(round(param.size(1)*1)); %default no shrinkage
cmd = append('xyzproj -axis ', param.tiltax, ' -width ',w,' -tiltfile tiltangles.txt ',prev,' ',tbase);
%-ray borks up smaller width completely, -constant makes beads slightly weird
prev = tbase; disp(cmd); [~] = evalc('system(cmd)'); %run command, capture console spam

if strcmp(type,'real') %electron detection and CTF
tilt = ReadMRC(prev); %load the projected tiltseries as a volume

order = 2; %electron detection changable order thing because i still don't know which is better!
if order==1 %dose first, hackjob dose increase to get the scaling to work
    detected = helper_electrondetect(tilt,param);
    WriteMRC(detected,pix,append('2_dosetilt_',base));
    convolved = helper_ctf(detected,param); %per-tilt ctf convolution
    prev = append('3_ctf_',base);
    WriteMRC(convolved,pix,prev); %save the convolved image for review
end
if order==2 %dose second, reduces pixel size effect so more generalized but values are arbitrary?
    convolved = helper_ctf(tilt,param); %per-tilt ctf convolution
    prev = append('3_ctf_',base);
    WriteMRC(convolved,pix,prev); %save the convolved image for review
    convolved = rescale(convolved,min(tilt,[],'all'),max(tilt,[],'all')-0); %fix negative CTF
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

thick = string(round(param.size(3)*1)); %w = string(param.size(1));
%reconstruct and rotate back into the proper space
cmd = append('tilt -tiltfile tiltangles.txt -width ',w,' -thickness ',thick,' ',prev,' temp.mrc'); 
disp(cmd); [~] = evalc('system(cmd)'); %run the recon after displaying the command
cmd = append('trimvol -rx temp.mrc ',append('5_recon_',base)); %#ok<NASGU>
[~] = evalc('system(cmd)'); %run the command and capture outputs from spamming the console

delete temp.mrc %remove temporary files after they are used for rotation
end