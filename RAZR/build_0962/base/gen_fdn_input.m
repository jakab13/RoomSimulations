% GEN_FDN_INPUT - Perform mapping of ISM-output matrices to FDN input channels. For matrices
% containing diffuse reflections, corresponding processing is performed.
%
% Usage:
%   [FDNinputmat, fdn_setup] = gen_fdn_input(...
%       ir, ism_setup, fdn_setup, ism_data, isd, room, op)
%
% Input:
%   ir              ir structure (see RAZR)
%   ism_setup       Output of GET_ISM_SETUP
%   fdn_setup       Output of GET_FDN_SETUP
%   ism_data        Struct containing metadata for image sources (output of
%                   IMAGE_SOURCE_MODEL)
%   room            Room structure (see RAZR)
%   op              Options structure (see GET_DEFAULT_OPTIONS)
%
% Output:
%   FDNinputmat     Signal matrix to be fed into FDN
%   fdn_setup       See input; some fields added

%------------------------------------------------------------------------------
% RAZR engine for Mathwork's MATLAB
%
% Version 0.96.2
%
% Author(s): Torben Wendt, Nico Goessling
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
