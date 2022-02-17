function [projected_coords, distances] = to_unit_sphere(coords)
% TO_UNIT_SPHERE - Project points onto the unit sphere.
%
% Usage:
%   projected_coords = to_unit_sphere(coords)
%
% Input:
%   coords              Three-column matrix of arbitrary points
%
% Output:
%   projected_coords    Coordinates of those points projected onto the unit sphere,
%                       with [0 0 0] mapped to [0 0 0].
%   distances           Euclidean distances of the points to the origin.

%------------------------------------------------------------------------------
% RAZR engine for Mathwork's MATLAB
%
% Version 0.96.2
%
% Author(s): Michael Schutte
%
% Copyright (c) 2014-2021, Torben Wendt, Steven van de Par, Stephan Ewert,
% University of Oldenburg, Germany.
%
% This work is licensed under the
% Creative Commons Attribution-NonCommercial-NoDerivs 4.0 International
% License (CC BY-NC-ND 4.0).
% To view a copy of this license, visit
% http://creativecommons.org/licenses/by-nc-nd/4.0/ or send a letter to
% Creative Commons, 444 Castro Street, Suite 900, Mountain View, California,
% 94041, USA.
%------------------------------------------------------------------------------


projected_coords = zeros(size(coords));
distances = sqrt(sum(coords .^ 2, 2));

scale = repmat(distances, 1, 3);
nonzero_mask = (distances ~= 0);

if any(nonzero_mask)
    projected_coords(nonzero_mask, :) = coords(nonzero_mask, :) ./ scale(nonzero_mask, :);
end
