function [rt, surfarea] = estimate_rt(room, measure)
% ESTIMATE_RT - Reverberation time estimation after Eyring or Sabine.
%
% Usage:
%   [rt, surfarea] = ESTIMATE_RT(room, [measure])
%
% Input:
%   room        room strcuture (see RAZR)
%   measure     'eyring'/'e' or 'sabine'/'s' for specifying the estimation measure
%               (optional, default: 'eyring')
%
% Output:
%   rt          Reverberation time in seconds for those frequency bands for which room.abscoeff is
%               specified
%   surfarea    Total wall surface area

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
    measure = 'eyring';
end

[Aeq, surfarea, room] = eq_abs_surfarea(room);

const = 24*log(10)/speedOfSound(room);

switch measure
    case {'e', 'eyring'}
        rt = const*prod(room.boxsize)./(-surfarea.*log(1 - mean(room.abscoeff, 1)));
    case {'s', 'sabine'}
        rt = const*prod(room.boxsize)./Aeq;
    case 'sabine_legacy'
        rt = 55.3/speedOfSound(room)*prod(room.boxsize)./aeq_legacy(room);
end

end % eof


function Aeq = aeq_legacy(room)
% A legacy code snippet to calculate the equivalent absorption area.
% The newer code in eq_abs_surfarea.m performs a matrix multiplication instead
% of a loop. The numeric results of both methods differ slightly from each other.

k = [1 2; 3 1; 2 3; 3 2; 1 3; 2 1];
Aeq = zeros(1, size(room.abscoeff, 2));
for n = 1:size(k, 1)
    Aeq = Aeq + room.abscoeff(n, :).*room.boxsize(k(n, 1))*room.boxsize(k(n, 2));
end

end % eof
