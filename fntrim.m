function vol = fntrim(vol)
vol = vol(any(vol ~= 0,[2 3]),any(vol ~= 0,[1 3]),any(vol ~= 0,[1 2]));
end