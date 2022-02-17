function op = get_default_options_v091
% GET_DEFAULT_OPTIONS_V091 - Returns default options used in RAZR 0.91
%
% Usage:
%   op = GET_DEFAULT_OPTIONS_V091
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
op.len_rt_factor = 2;
op.ism_jitter_type = 'cart_legacy';
op.ism_jitter_factor = 0.2;
op.ism_rand_start_order = 2;
op.spat_limit_filtrange = 1;
op.legacy_distce_atten = 1;
