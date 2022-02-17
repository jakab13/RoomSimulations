function ir = shape_irspec(ir, type)
% SHAPE_IRSPEC - Filter an ir to shape its frequency response.
%
% Usage:
%   ir = SHAPE_IRSPEC(ir, type)
%
% Input:
%   ir      Impulse response structure (see RAZR)
%   type    Filter type. Available:
%           'none', 'n': No filtering (default)
%           'pink', 'p': Pink shaped frequency response
%
% Output:
%   ir      Filtered impulse response
%
% See also: SOUNDIR, APPLY_RIR

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


switch type
    case {'none', 'n'}
        return;
    case {'pink', 'p'}
        for n = 1:length(ir)
            ir(n).sig = fftfilt(pink_pulse(ir(n).fs), ir(n).sig);
        end
    otherwise
        error('Unknown filter specification: %s', type);
end
