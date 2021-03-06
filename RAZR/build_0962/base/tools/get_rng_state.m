function scurr = get_rng_state
% GET_RNG_STATE - Get current state of the shared random number generator.
%
% Usage:
%   scurr = GET_RNG_STATE
%
% Output:
%   scurr       Current random number generator state
%               (its class depends on the current MATLAB version)
%
% See also: INIT_RNG, RNG, RAND

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


% Account for different syntax, depending on current MATLAB version:
if exist('rng', 'file') || exist('rng', 'builtin')
    scurr = rng;
else
    scurr = rand('state');
end
