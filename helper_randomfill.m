function [outarray, split] = helper_randomfill(inarray,set,iters,density,opt)
%[outarray, split] = helper_randomfill(inarray,set,iters,density)
%shared function for adding particles randomly, used for generating models and adding distractors
%
arguments %is validation needed for dedicated helper function?
    inarray (:,:,:) double
    set
    iters
    density = 0.4
    opt.type = 'object'
end
insize = size(inarray); counts = struct('s',0,'f',0); %initialize counts and get input size
counts.first =0; counts.sec = 0; counts.th =0;

if isstruct(set) %prep for set being a struct - wait, what if it isn't? it should just break!
namelist = [set(:).id]; %simplified collection of all ids instead of the former double loop
for i=1:numel(namelist)
    split.(namelist{i}) = zeros(size(inarray)); %initialize split models, duplicates made redundant
end
end

fprintf('Attempting %i %s placements:  ',iters,opt.type)

for i=1:iters
    which = randi(numel(set)); 
    particle = set(which).vol; name = set(which).id{1};
    
    switch set(which).type
        case {'single','complex'} %hopefully this is an OR as expected
            sumvol = sum( cat(4,set(which).vol{:}) ,4); %vectorized sum of all vols within the group
            
            for retry=1:3
            tform = randomAffine3d('Rotation',[0 360]);
            rot = imwarp(sumvol,tform);
            %{
            %might transition to imwarp for slight speed increase if it doesn't cause problems
            randvec = rand(1,3); randang = randi(360); %randomize rotation and alignment
            rot = imrotate3(sumvol,randang,randvec);
%             x = round(randi(size(inarray,1)+(size(rot,1)/3)*0)-(size(rot,1)/2));
%             y = round(randi(size(inarray,2)+(size(rot,2)/3)*0)-(size(rot,2)/2));
%             z = round(randi(size(inarray,3)+(size(rot,3)/3)*0)-(size(rot,3)/2));
%             loc = [x y z];
            %}
            loc = round( rand(1,3).*size(inarray)-size(rot)/2 );

            [inarray,err] = helper_arrayinsert(inarray,rot,loc,'overlaptest');
            %{
            if err==0 && retry==1 %if first placement works, bail out of loop check
                counts.first = counts.first+1;
            elseif err==0 && retry==2 %just checking
                counts.sec = counts.sec+1;
            elseif err==0 && retry==3 %just checking
                counts.th = counts.th+1;
            end
            %}
            if err==0, break; end
            end
            %if err==1, disp(counts), end
            
            counts.f = counts.f + err;
            if err==0 %on success, place in splits and working array
                counts.s=counts.s+numel(set(which).vol);
                [inarray] = helper_arrayinsert(inarray,rot,loc);
                for t=1:numel(set(which).vol) %rotate and place each component of complex
                    %rot = imrotate3(set(which).vol{t},randang,randvec); 
                    rot = imwarp(set(which).vol{t},tform);
                    split.(set(which).id{t}) = helper_arrayinsert(split.(set(which).id{t}),rot,loc);
                end
            end
            
        case 'bundle'
            if i<iters/5 || randi(numel(set(which).vol))==1 %have some scalable value to determine weight?
            [inarray, split, counts] = radialfill(inarray,set(which),18,split,counts);
            %increase iters by fraction of N to reduce runtime?
            %i = i+8;
            end
            
        case 'cluster'
            disp('oh boy this is going to bug out for sure')
    end
    %i = i+1;
    
    if 1<0 %janky internal function call
    [split, err, inarray, counts, loc] = fn_placement(inarray, split, particle, name, insize, [0 0 0], counts);
    end
    
    if strcmp(set(which).type,'cluster') && err==0 %maybe while loop to terminate early and reduce nesting?
        %q = 0;
        %still can only deal with a single input option, not a mixed cluster
        for j=1:12
            sz = max(size(particle));
            [x, y, z] = sph2cart(rand*360,rand*360,rand*sz*1.6+sz*0.8);
            vec = round([x y z]+loc);
            
            [split, ~, inarray, counts] = fn_placement(inarray, split, particle, name, insize, vec, counts);
            %q = err+q;
            %if q>10, disp('early end'), break,  end
        end
        
    end
    
    if rem(i,25)==0, fprintf('%i,',counts.s), end
    if rem(i,500)==0, fprintf('\n'), end
    
    if counts.s>iters/2 && rem(counts.s,6)==0 %filter to prevent the slower IF from running so often
    if nnz(inarray)/numel(inarray)>density, fprintf('Density limit reached.'), break, end
    end
    
end

fprintf('\nPlaced %i particles, failed %i attempted placements\n',counts.s,counts.f)
%verbose test for indicating when placements were successful
%fprintf('%i first successes, %i second successes, %i third successes\n',counts.first,counts.sec,counts.th)

outarray = zeros(insize);
splitnames = fieldnames(split);
for i=1:numel(splitnames)
    outarray = outarray+split.(splitnames{i});
end


end

%radial internal func
function [inarray,split,counts] = radialfill(inarray,bundle,n,split,counts)

for retry=1:4 %attempt placing the initial segment, currently 4 tries
which = randi(numel(bundle.vol));
primary = bundle.vol{which};
tform = randomAffine3d('Rotation',[0 360]);
primary = imwarp(primary,tform);
%{
% x = round(randi(size(inarray,1)+(size(primary,1)/3)*0)-(size(primary,1)/2));
% y = round(randi(size(inarray,2)+(size(primary,2)/3)*0)-(size(primary,2)/2));
% z = round(randi(size(inarray,3)+(size(primary,3)/3)*0)-(size(primary,3)/2));
% init = [x y z];
%}
init = round( rand(1,3).*size(inarray)-size(primary)/2 );

[inarray,err] = helper_arrayinsert(inarray,primary,init,'overlaptest');

if err==0
    counts.s=counts.s+1;
    [inarray] = helper_arrayinsert(inarray,primary,init);
    split.(bundle.id{which}) = helper_arrayinsert(split.(bundle.id{which}),primary,init);
    break
elseif err==1 && retry==3
    counts.f=counts.f+1;
end
end

if err==0 %don't try radial if primary not placed
[vec] = shapeaxis(primary);
vecfix = [vec(2) vec(1) vec(3)]; %infuriating vector refactoring for use in arrays
p = null(vecfix)'; %get horizontal components of normal plane, alternatively could use svd

%theta = rand(1,n)*2*pi; %randomized points %theta = (0:1/n:1)*2*pi; %fixed points
ri = regionprops3(imbinarize(rescale(primary)),'PrincipalAxisLength');
len = ri.PrincipalAxisLength(1);
ri = sum(ri.PrincipalAxisLength(2:3))/4; %initial object radius, 1/4 sum of diameters
rs = zeros(1,numel(bundle.vol));

%loop to tform the parts instead of doing them every time again
for i=1:numel(bundle.vol)
    bundle.vol{i} = imwarp(bundle.vol{i},tform);
    props = regionprops3(imbinarize(rescale(bundle.vol{i})),'PrincipalAxisLength');
    rs(i) = sum(props.PrincipalAxisLength(2:3))/4; %current obj radius, other half of radial distance
end

initcenter = size(primary)/2+init;
cum = 0; r = 0;
for i=1:n
    which = randi(numel(bundle.vol)); %randomize index of the bundle member to place
    sec = bundle.vol{which}; %extract the selected volume
    rot = imrotate3(sec,randi(360),vec); %randomly rotate around its axis for variability
    
    theta = rand*2*pi; %random points for variety
    rad = (r+ri+rs(which));
    radial = rad*(p(1,:)*cos(theta)+p(2,:)*sin(theta)); %generate radial vector with radius and null vecs
    slide = (rand-0.5)*len*vecfix; %random displacement along axis direction
    loc = round(initcenter-size(rot)/2-(radial+slide));
    
    [~,err] = helper_arrayinsert(inarray,rot,loc,'overlaptest');
    if err==0
        [inarray] = helper_arrayinsert(inarray,rot,loc);
        split.(bundle.id{which}) = helper_arrayinsert(split.(bundle.id{which}),rot,loc);
        cum=0; counts.s=counts.s+1;
    elseif cum>(n/2-1)
        break
    else
        cum = cum+1; r=r+(rs(which)+ri)*cum/20; counts.f=counts.f+1;
    end
    
end
end

end

%need efficient and straightforward generic function for single, complex, group/multiple and bundle initial

%superjanky internal function call that i will probably discard
function [split, err, inarray, counts, loc] = fn_placement(inarray, split, particle, name, insize, vec, counts)
    % nargin<6, vec = [0 0 0]; end
    
    randvec = rand(1,3); 
    randang = randi(360);
    rot = imrotate3(particle,randang,randvec);
    
    %{
    x = randi(insize(1)+round(size(rot,1)/3))-round(size(rot,1)/3);
    y = randi(insize(2)+round(size(rot,2)/3))-round(size(rot,2)/3);
    z = randi(insize(3)+round(size(rot,3)/3))-round(size(rot,3)/3);
    loc = [x y z];
    %}
    loc = round( rand(1,3).*size(inarray)-size(rot)/2 );
    if sum(vec)~=0, loc= vec; end
    %disp(loc)
    
    [inarray,err] = helper_arrayinsert(inarray,rot,loc,'nonoverlap');
    %if err==1, disp(counts), end
    %err = 0;
    
    counts.f = counts.f + err;
    if err==0
        counts.s=counts.s+1;
        split.(name) = helper_arrayinsert(split.(name),rot,loc);
        if rem(counts.s,25)==0, fprintf('%i,',counts.s), end
        if rem(counts.s,700)==0, fprintf('\n'), end
    end
    
end