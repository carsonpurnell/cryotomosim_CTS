function [memdat] = gen_mem_atom(sz,pix,param)
%
%
arguments
    sz
    pix
    param.num = 1:6
    param.frac = -1
    param.memsz = 1
end
box = sz*pix;

if numel(param.num)>2
    param.num = param.num(randi(numel(param.num))); % random member of set
elseif numel(param.num)==2
    param.num = sort(param.num); % sort to make sure no negatives
    param.num=randi(param.num(end)-param.num(1)+1)+param.num(1)-1; % random number between two integers
end

% individual mem definitions (mostly garbage still)
mdict(1) = struct('class','vesicle','thick',28,'thickvar',6,'size',1,'sphericity',0.9);
%mdict(2) = struct('class','er','thick',14,'thickvar',4,'size',0.6,'sphericity',0.2);
%mdict(3) = struct('class','membrane','thick',32,'thickvar',4,'size',2,'sphericity',0.1);
%mdict(4) = struct('class','mito','thick',35,'thickvar',3,'size',3,'sphericity',0.8);
% make a fixed dict and use a selector function to grab the target ones?

% derived vars
% also have an auto calc fallback for not using frac? 10% per vesicle?
%if numel(num)>1, num=randi(num(end)-num(1)+1)+num(1)-1; end % target range calc
% make an easier vector expansion randi selector? x:y randi(numel(num)) sort of deal?
if param.frac<0, param.frac = min(sqrt(param.num/10)/2,1); end % fallback computed fraction of vol
param.seeds = round(param.num/param.frac)+0; 
% number of seeds needed for given membrane number and coverage ratio

%param = struct('num',num,'frac',frac,'memsz',memsz,'seeds',seeds);

[minit,blobtable] = voronoiblobcells(box,param,mdict);

[atoms,memcell,normcell] = blob2mem(minit,blobtable,mdict);

memdat = atoms; % needs to not be stupid
end

