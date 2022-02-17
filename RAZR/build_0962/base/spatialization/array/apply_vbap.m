function outmat = apply_vbap(inmat, vbap_setup, speaker_distances, op, do_apply_distce)
% APPLY_VBAP - Apply vector base amplitude panning to a signal matrix
%
% Usage:
%   outmat = APPLY_VBAP(inmat, vbap_setup, speaker_distances, op, do_apply_distce)
%
% Input:
%   inmat               Input signal matrix
%   vbap_setup          Output of GET_VBAP_SETUP
%   speaker_distances   Distances (in m) of loudspeakers in array to receiver
%   op                  Options structure (see RAZR)
%   do_apply_distce     If true, apply speaker-distance delay and attenuation
%
% Output
%   outmat              Processed signal matrix
%
% See also: GET_VBAP_SETUP

%------------------------------------------------------------------------------
% RAZR engine for Mathwork's MATLAB
%
% Version 0.96.2
%
% Author(s): Josef Poppitz, Torben Wendt
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


len = size(inmat, 1);
num_speakers = length(speaker_distances);
outmat = zeros(len, num_speakers);

num_active_speakers = length(vbap_setup.speaker_idx);

if op.vbap_color_compen && ~isempty(vbap_setup.colcomp_ir)
    filtfun = @(h, n, x) (fftfilt(h(:, n), x));
else
    filtfun = @(h, n, x) (x);  % identity operation on last input parameter
end

if op.jp_reproduction_numerics
    outmat = apply_vbap_legacy(...
        inmat, vbap_setup, speaker_distances, op, do_apply_distce);
else
    if do_apply_distce
        [gain_factors, delay_samples] = get_speaker_distce_compen(...
            speaker_distances(vbap_setup.speaker_idx), op);
    else
        gain_factors = ones(num_active_speakers, 1);
        delay_samples = zeros(num_active_speakers, 1);
    end
    
    start_samples = delay_samples + 1;
    cut_samples = len - start_samples + 1;
    
    for n = 1:num_active_speakers
        outmat(start_samples(n):end, vbap_setup.speaker_idx(n)) = ...
            outmat(start_samples(n):end, vbap_setup.speaker_idx(n)) + ...
            filtfun(...
            vbap_setup.colcomp_ir, n, ...
            inmat(1:cut_samples(n), vbap_setup.speaker_input_idx(n)) * ...
            vbap_setup.gains(n) * gain_factors(n) ...
            );
    end
end
end


function outmat = apply_vbap_legacy(inmat, vbap_setup, speaker_distances, op, do_apply_distce)
% JP's original code, up to minor changes

len = size(inmat, 1);
num_speakers = length(speaker_distances);

outmat = zeros(len, num_speakers);

[numDirctns, numGainComponents] = size(vbap_setup.triangles);

for gainComponent = 1:numGainComponents
    for dirctn = 1:numDirctns
        out_tmp = inmat(:, dirctn) * vbap_setup.gains_2d(dirctn, gainComponent);
        
        if op.vbap_color_compen
            out_tmp = fftfilt(vbap_setup.h_filt(:, gainComponent, dirctn), out_tmp);
        end
        
        outmat(:, vbap_setup.triangles(dirctn, gainComponent)) = ...
            outmat(:, vbap_setup.triangles(dirctn, gainComponent)) + out_tmp;
    end
end

if do_apply_distce
    % compensate different loudspeaker distances to listener
    outmat = apply_distances(outmat, speaker_distances, op);
end
end


function signals = apply_distances(signal_matrix, speaker_distances, op)
% RENDER_ARRAY - Construct a per-loudspeaker signal matrix from a per-source matrix
% and a mapping between sources and loudspeakers.
%
% Usage:
%	signals = RENDER_ARRAY(signal_matrix, speaker_distances, op)
%
% Input:
%   signal_matrix       Signal matrix whose rows represent samples and whose columns
%                       represent individual (virtual) sound sources.
%   speaker_distances   Column vector of distances between the listening position
%                       and each target loudspeaker.
%   op                  Options structure (see GET_DEFAULT_OPTIONS).
%
% Output:
%   signals             The output signals per loudspeaker (rows: samples;
%                       columns: speakers); in essence, it is the matrix product
%                       signal_matrix * panning_matrix, but with gains and delays
%                       applied to compensate for different distances between
%                       the listening position and each loudspeaker.
%#lic:schutte_poppitz_wendt

[num_samples, num_speakers] = size(signal_matrix);

if op.jp_reproduction_numerics
    op.legacy_distce_atten = true;
end

[gain_factors, delay_samples] = get_speaker_distce_compen(speaker_distances, op);

% assemble the output matrix
signals = zeros(num_samples, num_speakers);
for i = 1:num_speakers
    start_sample = delay_samples(i) + 1;
    cut_sample = num_samples - start_sample + 1;
    signals(start_sample:end, i) = gain_factors(i) .* signal_matrix(1:cut_sample, i);
end
end
