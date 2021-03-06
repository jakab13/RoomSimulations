function t = timevec(arg1, fs)
% TIMEVEC - Time vector for RIR structure or time signal
%
% Usage:
%   t = TIMEVEC(ir)
%   t = TIMEVEC(len, fs)
%
% Input:
%   ir      RIR structure (see RAZR)
%   len     Signal length in samples
%   fs      Sampling rate
%
% Output:
%   t       Time vector (in sec. if fs is specified in Hz)

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


if isstruct(arg1)
    len = size(arg1.sig, 1);
    fs = arg1.fs;
else
    len = arg1;
end

t = (0:len-1)'/fs;
