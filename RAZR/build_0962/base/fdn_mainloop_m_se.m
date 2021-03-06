% FDN_MAINLOOP_M - Matlab implementation of feeback delay network (FDN) mainloop, same algorithm
% as implented in the MEX-file FDN_MAINLOOP.
%
% Usage:
%   [outmat, outsig] = FDN_MAINLOOP_M(in, m, b, a, A);
%
% Input:
%   in      Input signal (1 or length(m) channels, i.e. columns)
%   m       Vector containing delays in samples
%   b, a    Matrix containing coefficients for absorption filters (one filter per row)
%   A       Feedback matrix of size length(m) x length(m)
%
% Output:
%   outmat  Output signal as matrix; columns represent FDN channels
%   outsig  sum(outmat, 2)
%
% See also: FDN_MAINLOOP

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
