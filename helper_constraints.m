function vol = helper_constraints(vol,edges)
%vol = helper_constraints(vol,edges)
%outputs an empty copy of vol, with 1s written on certain borders define by edgges input
%1x3 char array as '&&&' in yxz(matlab) and xyz(IMOD) oreientation, + is top, - bottom, space none, & both
vec = cell(1,numel(edges));
for i=1:numel(edges)
    for j=1:3
        vec{j} = 1:size(vol,j); %initialize vector triad
    end
    switch edges(i) %write the subject vector to list the appropriate border indices
        case "+"
            vec{i} = size(vol,i);
        case "-"
            vec{i} = 1;
        case "&"
            vec{i} = [1,size(vol,i)];
        case " "
            vec{i} = []; %empty vector writes nothing when used so if check not needed
    end
    vol(vec{1},vec{2},vec{3}) = 1;
end
end