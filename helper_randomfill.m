function [outarray, split] = helper_randomfill(inarray,set,iters,vescen,vesvol,density,opt)
%[outarray, split] = helper_randomfill(inarray,set,iters,density,opt)
%shared function for adding particles randomly, used for generating models and adding distractors
arguments
    inarray (:,:,:) double
    set %either a single array struct or a cell array of array structs
    iters
    vescen = 0 %is definitely janky
    vesvol = 0 %the other part of the jank
    density = 0.4 %either single val or vector 
    opt.type = 'particle'
    opt.graph = 0
    opt.memvol = 0
end
%insize = size(inarray); 
counts = struct('s',0,'f',0); %initialize counts and get input size

if opt.graph==1 %graphical output of particles block
    try %first try to find a cts gui plot to ouput to
        gui = findall(0,'Name','ctsgui').RunningAppInstance.UIAxes;
    catch %if no cts gui found, output to a newly created figure axis
        gui = figure('Name','Placement Progress'); gui = gca;
        xlabel('Failed placements'); ylabel('Successful placements');
    end
    hold(gui,'off'); plot(gui,0,0,'.'); hold(gui,'on'); %clear any prior contents of the graph
end

%probably need to make this a double loop across cells of particle sets

%for ii=1:numel(layers)
namelist = [set(:).id]; %vector collection of all ids instead of the former double loop
for i=1:numel(namelist)
    split.(namelist{i}) = zeros(size(inarray)); %initialize split models of target ids
end
%end



% membrane setup stuff
if iscell(vesvol) %prep skeleton point map if provided for TMprotein
    ismem = 1; 
    memvol = sum( cat(4,vesvol{:}) ,4);
    bw = bwdist(~memvol); %calculate distances inside the shape
    mask = rescale(imgradient3(bw))>0.5; %generate an inverse mask that approximates the border, minus the mid
    skel = (bw.*~mask)>max(bw,[],'all')/2-1; %apply the mask to the distance map and threshold edge noise
    skel = ctsutil('edgeblank',skel,2);
    skel = bwareaopen(skel,20); %clean any remaining outlier points
    %can skel be used to find normal vectors without needing centroid stored array?
    memlocmap = skel; %store map of valid locations inside the membrane
    init = [0,0,1]; %initial required orientation for memprots
    %sliceViewer(skel); %it does work
    
    %inside/outside membrane localization maps
    nonmem = bwdist(inarray)>2; %locmap for all available area
    
    %find the largest component, assumed to be the background space
    CC = bwconncomp(nonmem);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [~,idx] = max(numPixels);
    mout = zeros(size(nonmem));
    mout(CC.PixelIdxList{idx}) = 1; %locmap for outside of vesicles
    
    min = nonmem-mout; %locmap for inside vesicles
    %sliceViewer(nonmem); figure(); sliceViewer(mout); figure(); sliceViewer(min);
    
    %numel(find(skel))/numel(skel) %check occupancy
else
    ismem = 0;
end
% membrane setup stuff end


diagout = zeros(size(inarray,1),size(inarray,2),0);

%make a double loop, possibly making the internal loop an internal function?
%for ww=1:numel(layers)
%set = layer{ww}
%layeriters = iters(min(ww,end));
%do minor cleanup of locmaps - removing islands, subtract the working array?
fprintf('Layer 1, attempting %i %s placements:  \n',iters,opt.type)
%etc
for i=1:iters
    which = randi(numel(set)); 
    particle = set(which).vol; 
    %precall more things so structs aren't called into so many times
    vols = set(which).vol; 
    flags = set(which).flags; 
    
    
    %put split/group placement box after the type switch for efficiency and to make complex/memplex/assembly
    %more general schemes
    
    %locmap switch for getting randomized lists of potentially valid points
    %might need to use if/else instead to fulfill multiple conditions - lacking membranes etc
    
    %locmaps into a struct so they can be directly selected by flags name indexing?
    %or use a funct for writing/reading the appropriate locmap? falls back to all when no membrane?
    %{
    switch set(which).type
        case {'memplex','membrane'} %placement into membrane
            locmap = memlocmap>0;
            
        case 'inmem' %inside vesicle volume
            locmap = min>0;
                
        case 'outmem' %only outside vesicles
            locmap = mout>0;
            
        otherwise %everything else goes into the global locmap
            locmap = bwdist(inarray)>2; %surprisingly very slow, just skip it? or faster method?
            %ind2sub not surprisingly is slow - fast vector replacement method?
            [x,y,z] = ind2sub(size(locmap),find(locmap>0)); %don't need >0, minor speed loss
            pts = [x,y,z];
    end
    %}
    %
    if ismem==1 && ismember(set(which).type,{'memplex','membrane','inmem','outmem'})
        switch set(which).type
            case {'memplex','membrane'} %placement into membrane
                locmap = memlocmap==1;
            case 'inmem' %inside vesicle volume
                locmap = min==1;
            case 'outmem' %only outside vesicles
                locmap = mout==1;
        end
        
        %do membrane stuff
    else
        %locmap = bwdist(inarray)>2; %surprisingly very slow, just skip it? or faster method?
        %ind2sub not surprisingly is slow - fast vector replacement method?
        %locmap = ~inarray; %avoiding bwdist for speed here
        locmap = inarray==0; %faster than logical somehow
        
        %do placement testing and verification here?, place into splits and working array separately
    end
    %}
    
    %do final placement based on if there is a tform~=0 or theta/ax ~=0?
    %fallthrough for complex/assembly too
    %inherit complex/assembly/sum from some earlier check, assembly should check variable sumvol anyway
    %place the stuff at the end, either the final rot sumvol for non-plex and the tf/rots/ax/theta for plexes
    
    %org for new flag-based system:
    %switch for placement location to get locmap (mem,ves,cytosol, or any)
    rflags = (flags(randperm(length(flags)))); %randomize flag order for multiloc usage
    locpick = matches(rflags,{'membrane','vesicle','cytosol'}); %check if any special loc found
    if ismem==1 && any(locpick)
        %this switch needs a better expression, multiple locations should work (for membrane+else)
        %membrane-centering not placed inside membrane does a neat near-membrane localization
        locpick = rflags(locpick); locpick = locpick{1};
        switch locpick
            case 'membrane' %placement into membrane
                locmap = memlocmap==1;
            case 'vesicle' %inside vesicle volume
                locmap = min==1;
            case 'cytosol' %only outside vesicles
                locmap = mout==1;
        end
    else
        locmap = inarray==0; %faster than logical somehow
    end
            
    %switch for group class (bundle, cluster, or single) for placing - also need one for mem?
    %thing to catch whichever flag is relevant?
    if any(matches(rflags,{'bundle','cluster'}))
        specialflag = rflags{matches(rflags,{'bundle','cluster'})};
    end
    %specialflag = rflags(matches(rflags,{'bundle','cluster'}));
    %can be empty, if statement or fill with something?
    
    %switch or if to fallthrough for placing either sumvol or individuals and with rot, tform, or ax/theta
    
    %placement switch for each particle class
    switch set(which).type
        %bundle first because it's going to break all the flags and needs an overhaul
        case 'bundle' %bundle placement got complicated, need to refactor the internal function
            if i<iters/5 || randi(numel(set(which).vol))==1 %have some scalable value to determine weight?
            [inarray, split, counts] = radialfill(inarray,set(which),18,split,counts);
            %increase iters by fraction of N to reduce runtime? can't modify i inside for loop
            end
            
        %cluster second because it also breaks flags and needs reworking
        case 'cluster' %need to move into call to cluster function like bundle has
            sub = randi(numel(particle)); %get random selection from the group
            [rot,~,loc,err] = testplace2(inarray,locmap,set(which).vol{sub},3);
            counts.f = counts.f + err;
            if err==0 %on success, place in splits and working array
                counts.s=counts.s+1;
                [inarray] = helper_arrayinsert(inarray,rot,loc);
                split.(set(which).id{sub}) = helper_arrayinsert(split.(set(which).id{sub}),rot,loc);
                %generate list of random points near loc
                num = 12;
                r = randn(1,num).*mean(size(set(which).vol{sub}))/2+mean(size(set(which).vol{sub}));
                az = rand(1,num)*360; el = rand(1,num)*360;
                [x,y,z] = sph2cart(az,el,abs(r));
                clusterpts = round([x;y;z])+loc';
                
                %loop through points and try to place them
                for j=1:size(clusterpts,2)
                    sub = randi(numel(particle)); %new random member
                    
                    tform = randomAffine3d('Rotation',[0 360]); %generate random rotation matrix
                    rot = imwarp(set(which).vol{sub},tform); %generated rotated particle
                    [inarray,errc] = helper_arrayinsert(inarray,rot,clusterpts(:,j),'nonoverlap'); %test place
                    if errc==0 %if nonoverlap record and add to split
                        counts.s=counts.s+1;
                        split.(set(which).id{sub}) = helper_arrayinsert(split.(set(which).id{sub}),rot,clusterpts(:,j));
                    end
                end
            end
            
        case {'inmem','outmem','single','group'} %universal for non-special non-complexes
            sub = randi(numel(particle));
            
            %{
            for retry=1:4 %implement as general tester again?
                tform = randomAffine3d('Rotation',[0 360]); %generate random rotation matrix
                rot = imwarp(set(which).vol{sub},tform); %generated rotated particle
                
                %r = randi(size(pts,1)); loc = pts(r,:); %get a test point
                loc = ctsutil('findloc',locmap);
                
                loc = round(loc-size(rot)/2); %shift to place by the COM
                %loc = round( rand(1,3).*size(inarray)-size(rot)/2 ); %randomly generate test position
                [~,err] = helper_arrayinsert(inarray,rot,loc,'overlaptest');
                if err==0, break; end
            end
            %}
            
            [rot,~,loc,err] = testcyto(inarray,locmap,set(which).vol{sub},4);
            
            counts.f = counts.f + err;
            if err==0 %on success, place in splits and working array
                counts.s=counts.s+1;
                [inarray] = helper_arrayinsert(inarray,rot,loc);
                %tmp = split.(set(which).id{sub}); %negligible
                %tmp = helper_arrayinsert(tmp,rot,com); %was ~25 with 506 iters? slower to split assignments
                %split.(set(which).id{sub}) = tmp; %~6 s extra
                split.(set(which).id{sub}) = helper_arrayinsert(split.(set(which).id{sub}),rot,loc); %faster
                if ismem==1 && strcmp(set(which).type,'inmem')
                    [min] = helper_arrayinsert(min,-rot,loc);
                elseif ismem==1 && strcmp(set(which).type,'outmem')
                    [mout] = helper_arrayinsert(mout,-rot,loc);
                end
            end
            
        case {'complex','assembly'} %all or multiple structured components of a protein complex
            %sumvol = sum( cat(4,set(which).vol{:}) ,4); %vectorized sum of all vols within the group
            sumvol = set(which).sumvol;
            
            [rot,tform,loc,err] = testplace2(inarray,locmap,sumvol,4);
            %[rot,tform,loc,err] = testplace(inarray,sumvol,3);
            
            counts.f = counts.f + err;
            if err==0 %on success, place in splits and working array
                counts.s=counts.s+numel(set(which).vol);
                [inarray] = helper_arrayinsert(inarray,rot,loc);
                members = 2:numel(set(which).vol);
                if strcmp(set(which).type,'assembly') 
                    members = members(randperm(length(members))); 
                    if numel(members)>1, members = members(randi(numel(members)+1):end); end
                end
                members = [1,members]; %#ok<AGROW>
                for t=members %rotate and place each component of complex
                    rot = imwarp(set(which).vol{t},tform);
                    split.(set(which).id{t}) = helper_arrayinsert(split.(set(which).id{t}),rot,loc);
                end
            end
            
        %{    
        case {'single','group'} %randomly select one particle from the group (including single)
            sub = randi(numel(particle)); %get random selection from the group
            [rot,~,loc,err] = testplace(inarray,set(which).vol{sub},3);
            
            counts.f = counts.f + err;
            if err==0 %on success, place in splits and working array
                counts.s=counts.s+1;
                [inarray] = helper_arrayinsert(inarray,rot,loc);
                split.(set(which).id{sub}) = helper_arrayinsert(split.(set(which).id{sub}),rot,loc);
            end
           %}
        %{
        case {'inmem','outmem'}
            sub = randi(numel(particle));
            if strcmp(set(which).type,'inmem') %get locmap depending on target location
                [x,y,z] = ind2sub(size(min),find(min>0));
            else
                [x,y,z] = ind2sub(size(mout),find(mout>0));
            end
            pts = [x,y,z]; %r = randi(size(pts,1));
            
            for retry=1:3
                tform = randomAffine3d('Rotation',[0 360]); %generate random rotation matrix
                rot = imwarp(set(which).vol{sub},tform); %generated rotated particle
                r = randi(size(pts,1)); loc = pts(r,:); %get a test point
                com = round(loc-size(rot)/2); %shift to place by the COM
                %loc = round( rand(1,3).*size(inarray)-size(rot)/2 ); %randomly generate test position
                [~,err] = helper_arrayinsert(inarray,rot,com,'overlaptest');
                if err==0, break; end
            end
            counts.f = counts.f + err;
            if err==0 %on success, place in splits and working array
                counts.s=counts.s+1;
                [inarray] = helper_arrayinsert(inarray,rot,com);
                split.(set(which).id{sub}) = helper_arrayinsert(split.(set(which).id{sub}),rot,com);
            end
        %}
        
        case {'memplex','membrane'}
            [particle] = ctsutil('trim',particle); %trims all vols according to their sum
            %centered = centervol(molc); 
            %centered = ctsutil('centervol',particle); %centers all vols in cells to the COM of the first
            %sumvol = sum( cat(4,centered{:}) ,4); %get sum volume (should pregenerate in _input)
            sumvol = set(which).sumvol;
            
            %{
            [x,y,z] = ind2sub(size(memlocmap),find(memlocmap>0)); %don't need >0, minor speed loss
            pts = [x,y,z]; %probably need to replace this block with something much faster
            r = randi(size(pts,1));
            loc = pts(r,:);
            %}
            %need testplacing for membrane localization
            loc = ctsutil('findloc',locmap);
            
            [k] = dsearchn(vescen,loc); %nearest vesicle center and distance to it
            
            %[ax,theta] = sphrot(init,fin,ori)
            targ = loc-vescen(k,:); %get target location as if from origin
            targ=targ/norm(targ); init = init(:)/norm(init); %unitize for safety
            
            rotax=cross(init,targ); %compute the normal axis from the rotation angle
            theta = acosd( dot(init,targ) );
            %end
            %[rotax,theta] = sphrot(init,loc,vescen(k,:));
            
            %sel = particles(randi(numel(particles))).tmvol;
            sel = sumvol;
            spinang = randi(180);
            spin = imrotate3(sel,spinang,init'); %rotate axially before transform to target location
            rot = imrotate3(spin,theta,[rotax(2),rotax(1),rotax(3)]);
            
            %tdest = inarray+memvol*0-vesvol{k}*1; %slow
            tdest = inarray-vesvol{k}; %faster
            
            %if i==500, sliceViewer(rescale(tdest)+skel*0+memlocmap); end
            
            com = round(loc-size(rot)/2);
            [~,err] = helper_arrayinsert(tdest,rot,com,'overlaptest');
            counts.f = counts.f+err; %counts.s = counts.s-err; %increment fails, should refactor s
            
            if err==0
                [inarray] = helper_arrayinsert(inarray,rot,com); %write sum to working array
                [memlocmap] = helper_arrayinsert(memlocmap,-imbinarize(rot),com); %reduce mem loc map
                if ismem==1 %&& strcmp(set(which).type,'inmem') %reduce inmem/outmem maps if present
                    [min] = helper_arrayinsert(min,-rot,com);
                    [mout] = helper_arrayinsert(mout,-rot,com);
                %elseif ismem==1 %&& strcmp(set(which).type,'outmem')
                end
                counts.s = counts.s+1; %increment success, bad old way need to deprecate
                %actually write to the split arrays
                if strcmp(set(which).type,'memplex')
                    members = 1:numel(particle);
                    %assembly rejigger member nums here
                    for t=members %rotate and place each component of complex
                        spin = imrotate3(set(which).vol{t},spinang,init'); 
                        rot = imrotate3(spin,theta,[rotax(2),rotax(1),rotax(3)]);
                        %rot = imwarp(set(which).vol{t},tform);
                        %need to do the rotation for each individual component
                        split.(set(which).id{t}) = helper_arrayinsert(split.(set(which).id{t}),rot,com);
                    end
                else %membrane only designation
                    split.(set(which).id{1}) = helper_arrayinsert(split.(set(which).id{1}),rot,com);
                    %just write the sumvol to the split array
                end
            end
            
            %update the locmaps at the end?
            
    end
    
    %if rem(i,25)==0, fprintf('%i,',counts.s), end
    if rem(i,round(iters/20))==0, fprintf('%i,',counts.s), end
    %if rem(i,600)==0, fprintf('\n'), end
    
    if rem(i,5)==0 && rem(counts.s,3)==0 %filter to prevent the slower IF from running so often
    if nnz(inarray)/numel(inarray)>density, fprintf('Density limit reached.'), break, end, end
    
    if opt.graph==1 %draw progress graph continuously when used
        plot(gui,counts.f,counts.s,'.'); drawnow;
    end
    
    %diagnostic filler image
    %{
    if err==0
        diagout(:,:,end+1) = inarray(:,:,end/2);
    end
    %}
end

%WriteMRC(diagout,10,'diagaccumarray.mrc');

fprintf('\nPlaced %i particles, failed %i attempted placements, final density %g\n',...
    counts.s,counts.f,nnz(inarray)/numel(inarray))

outarray = zeros(size(inarray)); splitnames = fieldnames(split);
for i=1:numel(splitnames)
    outarray = outarray+split.(splitnames{i});
end

end

%placement testing for cytosol proteins
function [rot,tform,loc,err] = testcyto(inarray,locmap,particle,retry)
for retry=1:retry
    tform = randomAffine3d('Rotation',[0 360]); %generate random rotation matrix
    rot = imwarp(particle,tform); %generated rotated particle
    %loc = round( rand(1,3).*size(inarray)-size(rot)/2 ); %randomly generate test position
    loc = ctsutil('findloc',locmap);
    loc = round(loc-size(rot/2));
    [inarray,err] = helper_arrayinsert(inarray,rot,loc,'overlaptest');
    if err==0, break; end
end
end

%placement testing for membrane proteins - TBD


%preliminary internal function for initial placement testing
function [rot,tform,loc,err] = testplace(inarray,particle,retry)
for retry=1:retry
    tform = randomAffine3d('Rotation',[0 360]); %generate random rotation matrix
    rot = imwarp(particle,tform); %generated rotated particle
    loc = round( rand(1,3).*size(inarray)-size(rot)/2 ); %randomly generate test position
    [inarray,err] = helper_arrayinsert(inarray,rot,loc,'overlaptest');
    if err==0, break; end
end
end

%updated testplace subfunct, need renaming. old one still needed for radial/cluster old stuff
%replace again to get flag system working efficiently and neatly
function [rot,tform,loc,err] = testplace2(inarray,locmap,particle,retry)
for retry=1:retry
    tform = randomAffine3d('Rotation',[0 360]); %generate random rotation matrix
    rot = imwarp(particle,tform); %generated rotated particle
    %loc = round( rand(1,3).*size(inarray)-size(rot)/2 ); %randomly generate test position
    loc = ctsutil('findloc',locmap);
    loc = round(loc-size(rot/2));
    [inarray,err] = helper_arrayinsert(inarray,rot,loc,'overlaptest');
    if err==0, break; end
end
end

%need to update radial to start normally, then only call the subfunct for the radials
%radial internal func
function [inarray,split,counts] = radialfill(inarray,bundle,n,split,counts)
which = randi(numel(bundle.vol));

%replacement for redundant initial test code, using 4 tries (might need more to balance out bundles)
[primary,tform,init,err] = testplace(inarray,bundle.vol{which},4);

if err==0
    counts.s=counts.s+1;
    [inarray] = helper_arrayinsert(inarray,primary,init);
    split.(bundle.id{which}) = helper_arrayinsert(split.(bundle.id{which}),primary,init);
elseif err==1
    counts.f=counts.f+1;
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
    rad = (r+ri+rs(which)); %get radial distance for the attempt from component radii and accumulated disp
    radial = rad*(p(1,:)*cos(theta)+p(2,:)*sin(theta)); %generate radial vector with radius and null vecs
    slide = (rand-0.5)*len*vecfix; %random displacement along axis direction
    loc = round(initcenter-size(rot)/2-(radial+slide));
    
    %[inarray,err] = helper_arrayinsert(inarray,rot,loc,'nonoverlap');
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