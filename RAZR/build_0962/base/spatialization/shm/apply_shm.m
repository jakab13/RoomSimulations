function [out_L, out_R] = apply_shm(inmat, setup, op)
% APPLY_SHM - Apply spherical head model to a signal matrix along columns
% This is just a wrapper for hsfilterMod.
%
% Usage:
%   out = APPLY_SHM(in, setup, op)
%
% Input:
%   in      Input signal matrix, containing time signals in columns
%   setup   Structure containing the following fields:
%           Required:
%               azim - Azimuth angles   (0: front, 90: left, -90: right)
%               elev - Elevation angles (0: at ears, 90: up, -90: down)
%           Optional fields:
%               filtrange_start_spls - Vector of length size(in,2) containing
%               time samples (for each column of 'in') at which the channel-wise
%               filtering starts. From this sample on, filtering is applied on
%               the doubled length of one HRIR. Default behavior: filtering on
%               the whole signal.
%   op      RAZR-options structure
%
% Output:
%   out_L, out_R    Filtered signal matrices for left and right ear, respectively
%
% See also: HSFILTERMOD, APPLY_HRTF

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


if ~isfield(setup, 'filtrange_start_spls')
    setup.filtrange_start_spls = [];
end

out_L = hsfilterMod(setup.azim, setup.elev, 90, ...
    op.shm_warpMethod, op.fs, inmat, setup.filtrange_start_spls);

out_R = hsfilterMod(setup.azim, setup.elev, -90, ...
    op.shm_warpMethod, op.fs, inmat, setup.filtrange_start_spls);
