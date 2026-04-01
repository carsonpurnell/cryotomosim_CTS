function cts_batch(sz,batchmod,batchsim,opt)
arguments
    sz
    batchmod
    batchsim
    opt.ideal = 0
    opt.batchname = ''
end
n = numel(batchmod);
for i=1:n
    %pmod.pix = batchmod{i}.pix; pmod.iters(2) = batchmod{i}.iters;
    suf = append(opt.batchname,'_',string(i));
    % model
    [cts,~,~,~,~,~,~,~,~,outfile] = cts_model_atomic(sz,batchmod{i},'suffix',suf,'dynamotable',1);
    [path,name,ext] = fileparts(outfile);
    outfile = fullfile(path,append(name,'.atom.mat')); %bake into sim function?
    
    % simulation
    batchsim{i}.pix = cts.param.pix;
    %tmp = namedargs2cell(batchsim{i});
    %tsim = param_simulate(tmp{:});
    cts_simulate_atomic(outfile,batchsim{i},'suffix',append('sim_',string(i)));
    % ideal sim run
    if isstruct(opt.ideal) % run ideal sim if argument given
        %should already be a consolidated param? or allow a cell array of values?
        isim = opt.ideal; isim.pix = cts.param.pix;
        cts_simulate_atomic(outfile,isim,'suffix',append('ideal_',string(i)));
    end
    fprintf('CTS batch: finished %i of %i runs\n',i,n)
end
%fprintf('done batch of %i runs\n',n)
end
