function [iced, ice] = gen_ice(vol,pix)

w = 8+2; %atomic number total in water molecule, until i can get electron opacity numbers
%amorphous ice density .94g/cm^3, but unit vol of liquid h20 molecule is 29.7a^3, estimate 20-30
%unit vol of water is 29.7a^3, dramatically larger than atoms, diameter ~2.75A.
density = 0.94/(1e8)^3/18*6.022e23; %convert from cm3 to a3, then g to molecules via mol/dalton

denspix = density*(pix^2.8);%/24; %^3 theoretical calculation correct, might be as low as 2.7
%different pixel size densities don't scale exactly cubic due to protein folding/surfacing
mol = round(denspix*numel(vol)*0.2); % 20% of ice mass randomly distributed as molecules (half molecules?)
ice = round(vol*0+denspix*0.80*w); %80% of ice mass as flat background for speed

pts = rand(3,mol).*size(vol)'; ix = round(pts+0.5); %fast pregeneration of all points
for i=1:mol
    co = ix(:,i); x = co(1); y = co(2); z = co(3); %slightly faster column indexing
    ice(x,y,z) = ice(x,y,z)+w*1;
   % if rem(i,round(mol/10))==0, fprintf(',%g',i/mol); end
end %fprintf('\n')

iced = max(ice,vol);
end