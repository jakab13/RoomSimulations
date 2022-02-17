function one_warning(msg)
% ONE_WARNING - Displays a warning only if the last warning had another message.
% Useful if the same warning thrown at multiple calls of the same code snippet.
%
% Usage:
%   ONE_WARNING(msg)
%
% Input:
%   msg     Warning message
%
% See also: WARNING, LASTWARN

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


if ~strcmp(lastwarn, msg)
    warning(msg);
end