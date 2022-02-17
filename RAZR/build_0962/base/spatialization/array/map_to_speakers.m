function [panning_matrix, closest_speakers] = map_to_speakers(source_coordinates, speaker_coordinates)
% MAP_TO_SPEAKERS - Map virtual sound sources to loudspeakers in an array
% based on proximity.
%
% Usage:
%	panning_matrix = MAP_TO_SPEAKERS(source_coordinates, speaker_coordinates, op)
%
% Input:
%   source_coordinates  Coordinates of the sound sources.
%   speaker_coordinates Coordinates of the loudspeakers.
%   op                  Options structure (see GET_DEFAULT_OPTIONS).
%
% Output:
%   panning_matrix      Mapping from sound sources to target loudspeakers,
%                       where a value of g in row i, column j means
%                       that loudspeaker j will play the signal of
%                       sound source i with a gain factor of g.

%------------------------------------------------------------------------------
% RAZR engine for Mathwork's MATLAB
%
% Version 0.96.2
%
% Author(s): Michael Schutte, Josef Poppitz, Torben Wendt
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


num_sources = size(source_coordinates, 1);
num_speakers = size(speaker_coordinates, 1);

if num_speakers == 0
    error('op.array_pos must not be empty if op.spat_mode is array');
end

panning_matrix = zeros(num_sources, num_speakers);

squared_distances = (1:num_speakers)';
closest_speakers = zeros(num_sources,1);

for i = 1:num_sources
    squared_distances(:, 2) = sum(bsxfun(...
        @minus, source_coordinates(i, :), speaker_coordinates) .^ 2, 2);
    sorted_squared_distances = sortrows(squared_distances,2);
    
    closest_speakers(i) = sorted_squared_distances(1,1);
    
    panning_matrix(i,closest_speakers(i)) = 1;
end

end
