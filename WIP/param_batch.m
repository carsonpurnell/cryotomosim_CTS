function param = param_batch(n,varargin)
% param = param_batch(n,varargin)
%inputs:
% n = number of parameter files to generate==number of batch runs
% varargin = same as inputs to param_model and param_simulate
%output:
% param = cell array of parameter structures from param_model or param_simulate
%
%note: must be run TWICE, one to generate each set of parameter structs
%param_simulate is used by inclusion of the 'dose' parameter, otherwise it feeds into param_model

if rem(numel(varargin),2)==1, error('CTS batch params: bad number of args'); end
param = cell(n,1);
ix = find(strcmp(varargin,'pix'),1); if isempty(ix); ix=0; end
doseix = find(strcmp(varargin,'dose'),1); if isempty(doseix); doseix=0; end
for i=1:n
    tmp = batchrand(varargin);
    if ix>0 && doseix==0
        pix = tmp{ix+1};
        param{i} = param_model(pix,tmp{:});
    else
        param{i} = param_simulate(tmp{:});
    end
end
end

%% internal function
function paramcell = batchrand(vars)
paramcell = struct(vars{:}); % place the initial parameters
f = fieldnames(paramcell);
for j=1:numel(f)
    if isequal(size(paramcell.(f{j})),[1,2]) && ~iscell(paramcell.(f{j}))
        % hideous randomization, subfunct for it?
        paramcell.(f{j}) = paramcell.(f{j})(1)+rand*(paramcell.(f{j})(2)-paramcell.(f{j})(1));
        paramcell.(f{j}) = round(paramcell.(f{j}),2,'significant');
    else
        %param{i}.(f{j}) = param{i}.(f{j}); % do nothing same val
    end
end
paramcell = namedargs2cell(paramcell);
end