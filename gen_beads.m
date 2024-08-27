function [beadstrc] = gen_beads(pix,radius)
%[beadstrc] = gen_beads(pix,radius)
%generates a CTS-format particleset of gold beads based on computed atomic density
%pix = pixel resolution of the beads in angstroms
%radius = vector of radii of different beads required. default 50
arguments
    pix (1,1) {mustBePositive}
    radius (1,:) = 50 %can use multiple as [40 60 120]
end
if isempty(radius), radius = 45; end %backstop for empty but extant radius

density = (19*6e23/197)*(pix/1e8)^3; %d = (d*mol/mass)*(pixel/m-a conv)^3 average atom/pix for gold 19g/cm3
Au = 79/2.5; %atomic number of gold, until atom opacities figured out
%probably don't need to do density/point tradeoff with mass/points being static across pixel sizes
%are real beads uniform density? hollow core might be from edge effect of uniform opacity
%lower border density might help close the hole

beadstrc.type = 'single'; beadstrc.flags = {'cytosol'};
for j=1:numel(radius)
    d = round(radius(j)*2/pix)+2; %generate empty volume for bead points
    cen = radius(j)/pix+2.5; %center for distance to be calculated from
    gold = zeros(d,d,d); %generate temporary array to fill with bead
    pts = rand(3,round(density*numel(gold))).*size(gold)'; %pregenerate points from computed atomic density
    atomdist = sqrt(sum( (pts-cen).^2, 1)); %calculate distance of each point for pruning
    pts = round(pts(:,atomdist<radius(j)/pix)+0.5); %remove atom points that would be outside sphere
    for i=1:size(pts,2)
        x = pts(1,i); y = pts(2,i);  z = pts(3,i);
        gold(x,y,z) = gold(x,y,z) + Au/2; %add gold intensity for each generated atom in the bead
        %half density to avoid giant CTF hole (and sometimes also peak!) in high mag beads
    end
    beadstrc.vol{j} = ctsutil('trim',gold); %write vol and id to the struct to pass to cts_model
    beadstrc.id{j} = append('bead__',string(radius(j)));
end

end