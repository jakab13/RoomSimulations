function op = get_default_options_jp
% GET_DEFAULT_OPTIONS_JP - Returns default options to restore behavior from
% J. Poppitz' Master's thesis.
%
% Usage:
%   op = GET_DEFAULT_OPTIONS_JP
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


op = get_default_options_v091;  % the version JP forked from

op.len_rt_factor = 1.5;
op.len_rt_measure = 'mean';
op.spat_limit_filtrange = 0;

op.hrtf_options.mk2_legacy_angle_conversion = 1;

op.jp_reproduction_numerics = 1;
