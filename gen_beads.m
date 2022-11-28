function [beadstrc] = gen_beads(pix,radius)
%no help yet you're on your own
%but use angstroms
%radius default 50, outputs the same number of different beads as input radii
arguments
    pix (1,1) {mustBePositive}
    radius (1,:) = 50 %can use multiple as [40 60 120]
end
%calculate gold atomic density in atom/pix from 19g/cm^3
%atomdensity = d/(cm-A conv)*avagadro/mass
density = (19/1e8^3)*6e23/197; %calculate atoms/A^3, ~.06
density = density*pix^3; %average atoms per pixel
Au = 79; %atomic number of gold
%probably don't need to do density/point tradeoff with mass/points being static across pixel sizes

%are real beads uniform density? hollow core might be from edge effect of uniform opacity
%lower border density might help close the hole

beadstrc(numel(radius)).type = 'single';
for j=1:numel(radius)
    d = round(radius(j)*2/pix)+2; %rad in A, diam in pixels
    cen = radius(j)/pix+1;
    
    gold = zeros(d,d,d);
    atoms = 1*round(density*numel(gold)*1); % double points at half mass for smoothness
    %probably don't need this now? can make it worse on purpose?
    
    pts = rand(3,atoms).*size(gold)'; %pregenerate points
    atomdist = sqrt(sum( (pts-cen).^2, 1)); %marginally faster
    pts = round(pts(:,atomdist<radius(j)/pix)+0.5); %remove atom points that would be outside sphere
    
    for i=1:size(pts,2)
        x = pts(1,i); y = pts(2,i);  z = pts(3,i);
        gold(x,y,z) = gold(x,y,z) + Au/1; %reducing gold signal to prevent protein overlap, needs work
    end
    gold = ctsutil('trim',gold);
    
    beadstrc(j).vol{1} = gold;
    beadstrc(j).id = append('bead__',string(radius(j)));
end

end