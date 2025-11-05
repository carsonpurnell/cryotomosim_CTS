%% placing memprots with new atomic membrane structure
pix = 10;
targ = {'ATPS.membrane.complex.mat'};
pmod = param_model(pix,'layers',targ);

sz = [200,200,50];
memdat = gen_mem_atom(sz,pix);
% needs a bit more work, a few vectors (probably due to corners) are not well-oriented

%%
for i=1:numel(pmod.layers)
    for j=1:numel(pmod.layers{i}.id)
        atoms.(pmod.layers{i}.id{j}) = zeros(0,4);
    end
end

%pick a membrane
memsel = randi(numel(memdat));


% start off preselecting coords from it or just start running through them? they are not spatially ordered
% inner loop: random axial rotation, rotation to transmembrane vector, collision test
% if no collision, switch to place subunits as needed