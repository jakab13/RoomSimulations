function pulse = pink_pulse(fs, len)
% PINK_PULSE - Generates an impulse response corresponding to a pink shaped
% frequency response.
%
% Usage:
%   pulse = PINK_PULSE([fs], [len])
%
% Input:
%   fs      Sampling rate in Hz (default: 44100)
%   len     Length in samples (default: fs)
%
% Output:
%   pulse   Impulse response signal

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


if nargin < 2
    len = [];
    if nargin < 1
        fs = 44100;
    end
end

if isempty(len)
    len = fs;
end

lcut = [1 50];
hcut = fs;
ntype = 1;
cweight = 0;
minphase = 1;
level = 0;
circ = 1;

pulse = mpulse(len, lcut, hcut, ntype, cweight, minphase, level, circ, fs);
