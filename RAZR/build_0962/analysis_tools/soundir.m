function soundir(ir, filtr)
% SOUNDIR - Plays RIR using soundsc.
%
% Usage:
%   SOUNDIR(ir, [filtr])
%
% Input:
%   ir      ir structure (see RAZR)
%   filtr   Optional specification of filter type to filter the ir. Available:
%           'none', 'n': No filtering (default)
%           'pink', 'p': Pink shaped frequency response
%
% See also: APPLY_RIR, SHAPE_IRSPEC

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
    filtr = 'none';
end

ir = shape_irspec(ir, filtr);
soundsc(ir.sig, ir.fs);
