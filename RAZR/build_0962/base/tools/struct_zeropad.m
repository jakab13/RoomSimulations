function [s, maxlen] = struct_zeropad(s, fldnames, ignore_empty)
% STRUCT_ZEROPAD - For a struct containing vectors and/or matrices zeropad all of
% them to match size along first dimension.
%
% Usage:
%   [s_out, maxlen] = STRUCT_ZEROPAD(s_in, [fldnames], [ignore_empty])
%
% Input:
%   s_in            Input structure
%   fldnames        Names of fields to be zeropadded. If an empty (default), all
%                   fields are taken into account
%   ignore_empty    If true, ignore all fields f with isempty(s.f) == true.
%                   Default: false.
%
% Output:
%   s_out       Output structure, all fields with adjusted lengths
%   maxlen      Actual lengths of fiels of s_out
%
% See also: STRUCT_REPMAT

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


if nargin < 3
    ignore_empty = false;
    if nargin < 2
        fldnames = {};
    end
end

if isempty(fldnames)
    fldnames = fieldnames(s);
end

% keep only existing fields:
fldnames = fldnames(isfield(s, fldnames));

if ignore_empty
    fldnames = rm_empty(fldnames, s);
end

for n = length(fldnames):-1:1
    [lens(n), numCols(n)] = size(s.(fldnames{n}));
end

maxlen = max(lens);

for n = find(lens < maxlen)
    s.(fldnames{n}) = [s.(fldnames{n}); zeros(maxlen - lens(n), numCols(n))]; 
end
end


function fldnames = rm_empty(fldnames, s)
% Remove names of empty fields

numFld = length(fldnames);
keep_idx = true(numFld, 1);

for f = 1:numFld
    if isempty(s.(fldnames{f}))
        keep_idx(f) = false;
    end
end

fldnames = fldnames(keep_idx);

end
