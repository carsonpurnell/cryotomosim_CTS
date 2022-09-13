function [label,roinames,indvol] = helper_particleatlas(cts,opt)

arguments
    cts struct %might make this work other than needing the cts struct
    opt.individual = 0 %by default don't save individual splits, only global atlas image
    opt.dynamotable = 0 %switch to generate dynamo table
    %opt.path %hoping that this inherits the path from higher-level functions
end

%move other important particles to splitmodel if they exist
if isfield(cts.model,'beads'), cts.splitmodel.beads = cts.model.beads; end
if isfield(cts.model,'mem'), cts.splitmodel.mem = cts.model.mem; end
if isfield(cts.model,'grid'), cts.splitmodel.grid = cts.model.grid; end

roinames = fieldnames(cts.splitmodel); %retrieve component names
indvol = cell(1,numel(roinames));

for i=1:numel(roinames)
    %filename = append('ind',string(i),'_',roinames{i},'.mrc');
    tmp = cts.splitmodel.(roinames{i});
    %if opt.bin==1; indvol=imbinarize(rescale(indvol{i})); end
    
    
    
    
    if opt.individual==1
        filename = append('ind',string(i),'_',roinames{i},'.mrc');
        WriteMRC(tmp,ts.pix,filename)
    end
    
end




end