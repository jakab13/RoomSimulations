function init_rng(seed)
% INIT_RNG - Initialize random number generator independently from your current
% MATLAB version.
%
% Usage:
%   INIT_RNG(seed)
%   INIT_RNG(state)
%   INIT_RNG
%
% Input:
%   seed        Seed for the random number generator (Mersenne Twister algorithm)
%   state       Output of GET_RNG_STATE
%   (no input)  seed will be set to 5489 (the MATLAB default for any version)
%
% Use INIT_RNG(sum(100*clock)) to get non-pseudo-random numbers.
%
% For Matlab version 7.12 or higher, calling this function is equivalent to
% RNG(seed). For older versions, it is equivalent to rand('twister', seed).
%
% See also: GET_RNG_STATE, RNG, RAND

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


if nargin == 0
    seed = 5489;	% Matlab default
end

% seed random generator (both methods should produce the same random numbers):
if exist('rng', 'file') || exist('rng', 'builtin')
    rng('default');             % prevents error when old method has been used before
    rng(seed);
else
    if size(seed, 1) > 1        % i.e. input is probably state and no seed
        rand('state', seed);
    else
        rand('twister', seed);
    end
end
