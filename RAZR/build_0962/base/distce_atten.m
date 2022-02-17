function atten = distce_atten(distce, use_legacy_method)
% DISTCE_ATTEN - "1/distance" attenuation factor.
%
% Usage:
%   atten = DISTCE_ATTEN(distce, [use_legacy_method])
%
% Input:
%   distce              Source-receiver distances in meters
%   use_legacy_method   If true, use legacy smooth transition to constant
%                       neutral gain for distances <= 1 (Default: false)
%
% Output:
%   atten               Attenuation factors

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
    use_legacy_method = false;
end

if use_legacy_method
    distc_shift = 0;  % in meters
    fadeConst = 0.9./(1 + exp(4*(distce + distc_shift - 1)));
    fadeHyperb = 1./(1 + exp(4*(1 - distce + distc_shift)));
    atten = min(fadeHyperb./(distce + distc_shift) + fadeConst, 1);
else
    atten = min(1 ./ distce, 1);
end
