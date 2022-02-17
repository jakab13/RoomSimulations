% CREATE_ISM_OUTPUT - Calculation of image source signals
%
% Uasge:
%   [out, isd] = CREATE_ISM_OUTPUT(isd, ism_setup, op)
%
% Input:
%   isd             Image source data; output of SCALE_IS_PATTERN
%   ism_setup       Output of GET_ISM_SETUP
%   op              Options structure (see GET_DEFAULT_OPTIONS)
%
% Output:
%   out             Structure with following fields:
%       sig             Spatialized time signal of image source pulses (contains
%                       only those specified by isd.idx_auralize)
%       sig_post_spat   If op.enable_post_spat == true: Binaural synthesis of
%                       array-spatialized signal. Otherwise: empty matrix.
%       sigmat          Not spatialized time signals of specular reflections
%                       (one image source per column). If
%                       op.ism_diffr_mc_output == false, diffraction was neither
%                       applied on these signals
%       sigmat_diffuse  Time signals of diffuse reflections, same format as
%                       signalmat
%   isd             Same as input, but with some additional fields.
%       filter_ranges   Restriction of filtering to specified signal intervals
%                       in order to save computations. Matrix of the form
%                       [start_ch1, end_ch1, ...; start_chN, end_chN].

%------------------------------------------------------------------------------
% RAZR engine for Mathwork's MATLAB
%
% Version 0.96.2
%
% Author(s): Torben Wendt, Nico Goessling, Oliver Buttler
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
