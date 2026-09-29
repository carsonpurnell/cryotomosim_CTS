function cts_getcif(id,param)
% downloads structures from the rcsb to the local machine for use with CTS.
% downloaded structures are usable immediately, if the destination folder is on/added to path
% inputs:
% id (required) list of structures to download. list of strings as ["1tub","9sln"]
% params: all name-value pairs
%   savepath (default null/[] for matlab default userpath) full path of saved structures
%   format (.cif or .pdb) format to download and save, default is '.cif'

arguments
    id
    param.savepath = []
    param.format {mustBeMember(param.format,{'.cif','.pdb'})} = '.cif'
end
if isempty(param.savepath) 
    param.savepath = fullfile(userpath,'rcsb_cif'); % for cross-OS compatability
end
if isempty(dir(param.savepath))
    mkdir(param.savepath)
end
prev = pwd;
cd(param.savepath)
url = 'https://files.rcsb.org/download/';

for i=1:numel(id)
    dl = append(url,id(i),param.format);
    dat = urlread(dl);
    fn = append(id(i),param.format);
    fid = fopen(fn,'w');
    fprintf(fid,dat);
    fclose(fid);
    
end

cd(prev); % go to previous folder to avoid potential weirdness after running funct
end