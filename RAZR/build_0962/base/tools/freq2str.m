function str = freq2str(freq, do_Hz)
% FREQ2STR - Convert vector of frequency values into cell array of strings,
% where numbers above 1000 are written as 1k etc.
%
% Usage:
%	str = FREQ2STR(freq, [do_Hz])
%
% Input:
%   freq    Frequencies as vector of doubles
%   do_Hz   If true, add unit Hz (optional, default: false)
%
% Output:
%   str     Frequencies as cell array of strings
%
% Example:
%   >> freq = [62.5, 125, 1e3, 16e3, 44100];
%   >> str = freq2str(freq)
%   str = 
%       '62.5'    '125'    '1k'    '16k'    '44.1k'
%   >> str = freq2str(freq, true)
%   str = 
%       '62.5 Hz'    '125 Hz'    '1 kHz'    '16 kHz'    '44.1 kHz'

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


if nargin == 1
    do_Hz = false;
end

tsd = {'', 'k'};

numFrq = length(freq);
str = cell(size(freq));

if do_Hz
    for n = 1:numFrq
        str{n} = sprintf('%g %sHz', ...
            freq(n)*(1e-3)^(freq(n)>=1e3), tsd{(freq(n) >= 1e3) + 1});
    end
else
    for n = 1:numFrq
        str{n} = sprintf('%g%s', ...
            freq(n)*(1e-3)^(freq(n)>=1e3), tsd{(freq(n) >= 1e3) + 1});
    end
end
