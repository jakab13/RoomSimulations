function op = get_default_options_v092
% GET_DEFAULT_OPTIONS_V092 - Returns default options used in RAZR 0.92
%
% Usage:
%   op = GET_DEFAULT_OPTIONS_V092
%
% See also: GET_DEFAULT_OPTIONS, GET_DEFAULT_OPTIONS_*

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


op = get_default_options;
op.legacy_distce_atten = 1;
