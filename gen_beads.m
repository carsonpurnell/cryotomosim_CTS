function [beadstrc] = gen_beads(pix,radius)
%no help yet you're on your own
%but use angstroms
%radius default 50, outputs the same number of different beads as input radii
arguments
    pix (1,1) {mustBePositive}
    radius (1,:) = 50 %can use multiple as [40 60 120]
end
if isempty(radius), radius=50; end

%calculate gold atomic density in atom/pix from 19g/cm^3
density = (19/1e8^3)*6.022e23/197; %calculate atoms/A^3
density = density*pix^2.7; %average atoms per pixel
Au = 79; %atomic number of gold
%probably don't need to do density/point tradeoff with mass/points being static across pixel sizes

beadstrc(numel(radius)).id = {''};
for j=1:numel(radius)
    d = round(radius(j)*2/pix)+2; %rad in A, diam in pixels
    cen = radius(j)/pix+1;
    
    gold = zeros(d,d,d);
    atoms = 2*round(density*numel(gold)*1); % double points at half mass for smoothness
    %probably don't need this now? can make it worse on purpose?
    
    pts = rand(3,atoms).*size(gold)'; %pregenerate points
    atomdist = sqrt(sum( (pts-cen).^2, 1)); %marginally faster
    pts = round(pts(:,atomdist<radius(j)/pix)+0.5); %remove atom points that would be outside sphere
    
    for i=1:size(pts,2)
        x = pts(1,i); y = pts(2,i);  z = pts(3,i);
        gold(x,y,z) = gold(x,y,z) + Au/4; %reducing gold signal to prevent protein overlap, needs work
    end
    
    gold = gold(any(gold ~= 0,[2 3]),any(gold ~= 0,[1 3]),any(gold ~= 0,[1 2])); 
    %{
    gold = gold(:,any(gold ~= 0,[1 3]),:); 
    gold = gold(any(gold ~= 0,[2 3]),:,:); 
    gold = gold(:,:,any(gold ~= 0,[1 2]));
    %}
    
    beadstrc(j).vol{1} = gold;
    beadstrc(j).id = append('bead__',string(radius(j)));
    beadstrc(j).type = 'single';
end

end