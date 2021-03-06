function varargout = adjust_len(varargin)
% ADJUST_LEN - Zeropad smaller matrices to match size of largest matrix, with respect to first
% dimension.
%
% Usage:
%   [A_out, B_out, ...] = adjust_len(A_in, B_in, ...)
%
% Input:
%   A_in, B_in, ...     Matrices of different sizes
%
% Output:
%   A_out, B_out, ...   Zero-padded input matrices
%

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


lens = zeros(nargin, 1);
numcols = zeros(nargin, 1);

for n = 1:nargin
    [lens(n), numcols(n)] = size(varargin{n});
end

maxlen = max(lens);
varargout = varargin;

for n = 1:nargin
    varargout{n} = [varargin{n}; zeros(maxlen - lens(n), numcols(n))];
end
