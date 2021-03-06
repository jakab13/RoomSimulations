% COMPOSEDPEQ - Get filter coefficients for composed parametric EQs, to approximate specified
% frequency response. The resulting filter has the order length(freq)*2+1.
%
% Usage:
%   [B, A, b, a] = COMPOSEDPEQ(freq, dGain, fs, linGain, [do_optimization])
%
% Input:
%   freq            Frequency base in Hz, must have the same length as dGain
%   dGain           Desired frequency response, values in dB or linear (specified by input parameter 'linGain')
%   fs              Sampling rate
%   linGain         Set true, if dGain values are linear attenuation factors, false if they are specified in dB
%   do_optimization If true, apply performance optimization by summarize similar sequent dGain values
%                   (optional, default: true)
%
% Output:
%   B, A            Filter coefficients of the shelving filters as Matrix (each filter in one row)
%   b, a            Filter coefficient vectors of the whole EQ

%------------------------------------------------------------------------------
% RAZR engine for Mathwork's MATLAB
%
% Version 0.96.2
%
% Author(s): Torben Wendt
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
