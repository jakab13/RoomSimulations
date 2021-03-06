function op = get_default_options_2014
% GET_DEFAULT_OPTIONS_2014 - Returns options structure for BRIR synthesis according to the method
% described in
% Wendt, van de Par, and Ewert (2014): A Computationally-Efficient and Perceptually-Plausible
% Algorithm for Binaural Room Impulse Response Simulation, JAES vol. 62.
%
% Usage:
%   op = GET_DEFAULT_OPTIONS_2014
%
% Output:
%   op      Options structure (see RAZR)
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

%% general options

op.len_rt_factor = 2;
op.spat_mode = 'hrtf';
op.hrtf_database = 'mk2';
op.hrtf_filtfunc = 'filter';
op.spat_limit_filtrange = 1;
op.filtCreatMeth = 'cs';
op.rt_estim = 'sabine_legacy';
op.legacy_rand = 1;
op.legacy_distce_atten = 1;

%% ism

op.ism_enableAirAbsFilt = 0;
op.ism_enableToneCorr = 0;
op.ism_enableBP = 1;
op.legacy_bp = 1;

op.ism_jitter_factor = 0.2;
op.ism_jitter_type = 'sph';
op.ism_rand_start_order = 0;

%% fdn

op.fdn_legacy_input = 1;
op.fdn_enableBP = 1;
op.fdn_delays_choice = 'roomdim2014';
