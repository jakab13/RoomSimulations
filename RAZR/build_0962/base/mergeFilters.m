% merge filter coeffs to filters of maxOrder (default 11)
%
% SE 4/3/2020 5:22:00 PM

%------------------------------------------------------------------------------
% RAZR engine for Mathwork's MATLAB
%
% Version 0.96.2
%
% Author(s): Stephan D. Ewert
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

function [bb, aa] = mergeFilters(b, a, maxOrder, direction)

% ck 20-04-08 quick and dirty format matching
a = a';
b = b';

if nargin < 4
   % start from end
   direction = -1;
end

if nargin < 3
  maxOrder = 11;
end

if ~iscell(b)
    b = mat2cell(b, ones(1, size(b, 1)), size(b, 2));
else
    b = b(:);
end
if ~iscell(a)
    a = mat2cell(a, ones(1, size(a, 1)), size(a, 2));
else
    a = a(:);
end

if direction == -1
    b = flipud(b);
    a = flipud(a);
end

done = 0;
mergeIdx = 1;
startIdx = 1;
maxIdx = length(b);

while done == 0
    
    bb{mergeIdx,1} = b{startIdx };
    aa{mergeIdx,1} = a{startIdx };
    
    if startIdx == maxIdx
        break;
    end
    
    for k = startIdx+1:maxIdx
        
        % next round gets too long
        if length(bb{mergeIdx,1}) + length(b{k}) - 1 > maxOrder
            break;
        end
        
        bb{mergeIdx,1} = conv(bb{mergeIdx,1}, b{k});
        aa{mergeIdx,1} = conv(aa{mergeIdx,1}, a{k});
        
        % got until the end, done
        if k == maxIdx
            done = 1;
            break;
        end

    end
    
%     if done == 1
%         break;
%     end

    mergeIdx = mergeIdx + 1;
    startIdx = k;
end

if direction == -1    
    bb = flipud(bb);
    aa = flipud(aa);
end

% ck 20-04-08 quick and dirty format matching
bb = bb';
aa = aa';
