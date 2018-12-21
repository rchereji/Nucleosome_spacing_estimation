function Occ = extend_dyads_to_footprint(Dyads, a)
% Converts a cell array containing the nucleosome dyad distribution into a
% cell array containing the coverage profile obtained by stacking extended
% footprints, symmetricaly extended from the dyad position (nucleosome
% occupancy profiles, if a = 147 bp).
%
% Inputs:
%   Dyads - cell array, each cell corresponding the nucleosome dyad counts
%           for a separate chromosome
%   a     - footprint of the particle (typical size for nucleosomes = 147)
%
% Output:
%   Occ   - cell array, each cell corresponding the nucleosome occupancy
%           for a separate chromosome
%
% Examples:
% Occ = extend_dyads_to_footprint(Dyads, 147) % compute nucleosome occupancy
% Occ = extend_dyads_to_footprint(Dyads, 101) % extend dyads to 101 bp

halfParticle = floor(a/2);
noChr = numel(Dyads);
Occ = cell(size(Dyads));
for c = 1:noChr
    Occ{c} = filter(ones(1,a), 1, [Dyads{c}(halfParticle + 1 : end), zeros(1, halfParticle)]);
end
