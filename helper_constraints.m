function vol = helper_constraints(vol,edges)
%vol = helper_constraints(vol,edges)
%vol is output as a zero array of the same size, with 1s written into certain borders defined by edges
%edges is a 1x3 char vector of the chars [& +-], defining the borders of those dimensions
%+ makes a border on the high index, - on the low index, & on both, and ' '/space adds no border
vec = cell(1,numel(edges));
for i=1:numel(edges)
    for j=1:3
        vec{j} = 1:size(vol,j); %initialize vector triad
    end
    switch edges(i) %write the subject vector to list the appropriate border indices
        case "+", vec{i} = size(vol,i);
        case "-"; vec{i} = 1;
        case "&", vec{i} = [1,size(vol,i)];
        case " ", vec{i} = []; %empty vector writes nothing when used so if check not needed
    end
    vol(vec{1},vec{2},vec{3}) = 1;
end
end