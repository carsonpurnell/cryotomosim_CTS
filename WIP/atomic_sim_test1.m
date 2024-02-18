
% collect atomic model into single set of points
fn = fieldnames(split);
atoms = zeros(0,4);
for i=1:numel(fn)
    atoms = [atoms;split.(fn{i})];
end

% rotate the model to an angle (eucentric adjustment? rotate about 0?)

% project a set of slices - higher resolution in Z? start with isotropy

% get the transmission wave

% propogate transmission

