function signals = render_array(signal_matrix, panning_matrix, speaker_distances, op, do_distce)
% RENDER_ARRAY - Construct a per-loudspeaker signal matrix from a per-source matrix
% and a mapping between sources and loudspeakers.
%
% Usage:
%	signals = RENDER_ARRAY(signal_matrix, panning_matrix, speaker_distances, op)
%
% Input:
%   signal_matrix       Signal matrix whose rows represent samples and whose columns
%                       represent individual (virtual) sound sources.
%   panning_matrix      Mapping from sound sources to target loudspeakers,
%                       where a value of g in row i, column j means
%                       that loudspeaker j will play the signal of
%                       sound source i with a gain factor of g.
%   speaker_distances	Column vector of distances between the listening position
%                       and each target loudspeaker.
%   op                  Options structure (see GET_DEFAULT_OPTIONS).
%
% Output:
%   signals             The output signals per loudspeaker (rows: samples;
%                       columns: speakers); in essence, it is the matrix product
%                       signal_matrix * panning_matrix, but with gains and delays
%                       applied to compensate for different distances between
%                       the listening position and each loudspeaker.

%------------------------------------------------------------------------------
% RAZR engine for Mathwork's MATLAB
%
% Version 0.96.2
%
% Author(s): Michael Schutte, Torben Wendt
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


num_samples  = size(signal_matrix,  1);
num_speakers = size(panning_matrix, 2);

if do_distce
    [gain_factors, delay_samples] = get_speaker_distce_compen(speaker_distances, op);
else
    gain_factors = ones(num_speakers, 1);
    delay_samples = zeros(num_speakers, 1);
end

if op.jp_reproduction_numerics && do_distce
    % Schutte:
    % compensate for the attenuation of sound level with distance as per the inverse-square law,
    % with a gain factor of 1 for speakers 1 m away from the listening position
    gain_factors = speaker_distances .^ 2;
end

% assemble the output matrix
signals = zeros(num_samples, num_speakers);
for i = 1:num_speakers
    start_sample = delay_samples(i) + 1;
    cut_sample = num_samples - start_sample + 1;
    signals(start_sample:end, i) = gain_factors(i) * ...
        (signal_matrix(1:cut_sample, :) * panning_matrix(:, i));
end
