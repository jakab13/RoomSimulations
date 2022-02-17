% GET_REFL_BLEND_SETUP - Generate virtual reverb source positions for blended
% reflection coefficients.
%
% Usage:
%   sph = GET_REFL_BLEND_SETUP(room, fdn_setup, op, [do_plot])
%
% Input:
%   room        Room structure (see RAZR)
%   fdn_setup   Output of GET_FDN_SETUP
%   op          Options structure (see RAZR)
%   do_plot     If true, plot reverb source positions (morphed, at intermediate
%               state of calculation)
%
% Output:
%   sph         Structure containing source positions and weighting setup (using
%               VBAP)

%------------------------------------------------------------------------------
% RAZR engine for Mathwork's MATLAB
%
% Version 0.96.2
%
% Author(s): Josef Poppitz, Torben Wendt
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
