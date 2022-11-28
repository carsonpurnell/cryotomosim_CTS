function binned = helper_binvol(input, binning, verbose)
%binned = tomosim_bin(input, binning, verbose)
%input mrc volume required, binning is optional with default 2 - must be integer if supplied
%currently unused demo code
arguments
    input (1,1) string %mustBeFile requires 2020b for this validator
    binning (1,1) {mustBeInteger, mustBePositive} = 2 %ensure positive integer, default to 2
    verbose (1,1) {mustBeInteger} = 0
end

call = 'binvol -binning '; 
root = strrep(input,'.mrc','');
out = append(root,'_bin',num2str(binning),'.mrc'); %append _bin and binning to output filename

cmd = append(call,num2str(binning),' ',input,' ',out);
fprintf('Command passed: %s\n',cmd) %display the passed command for validation
dump = evalc('system(cmd)'); %capture imod output from displaying
if verbose==1, disp(dump), end %display imod output with verbose option

binned = ReadMRC(out);
end