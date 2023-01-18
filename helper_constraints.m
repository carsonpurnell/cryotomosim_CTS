%edge write for constraint generation function
function vol = helper_constraints(vol,edges)

%edges = '&&&';
%vol = zeros(20,20,10);
% + side, - side, neither, or & for both?
%could instead do single char '&=- ' sort

%make the [1,end] vector into the needed version for each entry first then use them?
%and empty [] vector writes nothing, usefully
vec = cell(1,numel(edges));
for i=1:numel(edges)
    for j=1:3
        vec{j} = 1:size(vol,j);
    end
    switch edges(i)
        case "+"
            vec{i} = size(vol,i);
        case "-"
            vec{i} = 1;
        case "&"
            vec{i} = [1,size(vol,i)];
        case " "
            vec{i} = [];
    end
    vol(vec{1},vec{2},vec{3}) = 1;
    %do the write thing
    
end

end