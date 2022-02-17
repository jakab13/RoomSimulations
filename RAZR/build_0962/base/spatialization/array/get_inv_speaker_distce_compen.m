function [gains, delay_samples] = get_inv_speaker_distce_compen(speaker_distances, op)
% GET_INV_SPEAKER_DISTCE_COMPEN - Get gain and delay to compensate for
% speaker-listener distance.
%
% Usage:
%   [gains, delay_samples] = GET_SPEAKER_DISTCE_COMPEN(speaker_distances, op)
%
% Input:
%   speaker_distances   Column vector of distances between the listening position
%                       and each target loudspeaker.
%   op                  Options structure (see GET_DEFAULT_OPTIONS).
%
% Output:
%   gains               Gain factors
%   delay_samples       Delays in samples

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


gains = distce_atten(speaker_distances, op.legacy_distce_atten);


% compensate for the delays due to different distances between the individual speakers
% and the listening positions
delays = speaker_distances / speedOfSound(op.array_TCelsius);
delays = delays - min(delays);
delay_samples = round(op.fs * delays);
