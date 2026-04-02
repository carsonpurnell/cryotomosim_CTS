function cts_batch(sz,batchmod,batchsim,opt)
arguments
    sz
    batchmod
    batchsim
    opt.method {mustBeMember(opt.method,['atom','atomic','vol'])} = 'atomic'
    opt.ideal = 0
    opt.batchname = ''
end
n = numel(batchmod);
for i=1:n
    %pmod.pix = batchmod{i}.pix; pmod.iters(2) = batchmod{i}.iters;
    suf = append(opt.batchname,'_',string(i));
    % model
    if strcmp(opt.method,'vol')
        [cts,outfile] = cts_model(zeros(sz),batchmod{i},'suffix',suf);
    else
        [cts,~,~,~,~,~,~,~,~,outfile] = cts_model_atomic(sz,batchmod{i},'suffix',suf,'dynamotable',1);
    end
    
    % simulation
    batchsim{i}.pix = cts.param.pix;
    %tmp = namedargs2cell(batchsim{i}); %tsim = param_simulate(tmp{:});
    if strcmp(opt.method,'vol')
        cts_simulate(outfile,batchsim{i},'suffix',append('sim_',string(i)));
    else
        [path,name,ext] = fileparts(outfile);
        outfile = fullfile(path,append(name,'.atom.mat')); %bake into sim function?
        cts_simulate_atomic(outfile,batchsim{i},'suffix',append('sim_',string(i)));
    end
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
