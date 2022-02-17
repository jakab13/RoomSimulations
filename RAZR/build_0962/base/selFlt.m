function [y, zf] = selFlt(b,a,x,zi,op)
% Function enables selecting from different filtering strategies. 
% Not all of these strategies
% are optimised for performance, e.g. filter coefficients are convolved 
% repeatedly. In the current implementation, flexibility is key.
% If no commitment can be made in future relases, a different data passover
% structure could be implemented in order to regain some performance.
% Input:
% b  -  nominator coefficients, cell array with n x 1, n denoting different
%       stages of the cascade
% a  -  denominator coefficients, cell array n x 1
% x  -  Input signal, vector or matrix, will be filtered along first dimension
% zi -  filter states, currently only works for 'conv'
% op -  RAZR options struct, important parameter op.fltMode
%       'conv' - convolve coefficients
%       'casc' - run cascaded filters in a loop
%       'sos'  - tbd
%       'smartCasc' - tbd
% Output
% y  - filtered signal
% zf - filter states (only works for 'conv' atm)
% TO DO:
% - make states work for cascade (save state for every part of the cascade)
% 20-03-30 CK

%------------------------------------------------------------------------------
% RAZR engine for Mathwork's MATLAB
%
% Version 0.96.2
%
% Author(s): Christoph Kirsch
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



%     if strcmp(op.fltMode,'conv')
%         bb = 1;
%         aa = 1;
%         for k = 1:length(a)
%             bb = conv(bb, b{k});
%             aa = conv(aa, a{k});
%         end
%         [y,zf] = filter(bb,aa,x,zi);
    if strcmp(op.fltMode,'casc') || strcmp(op.fltMode,'conv') || strcmp(op.fltMode,'smart')
        for k = 1:length(a)
            x = filter(b{k},a{k},x);
        end
        zf = [];%%%%%% !!!
        y = x;
    elseif strcmp(op.fltMode,'sos')
        % CK 03-20 check options for plausibility (yw and cs might not work
        % as SOS
%         for k = 1:length(b)
%             [sos(k,:),g(k)] = tf2sos(b{k},a{k}); % Generating SOS matrix
%         end
%     sos(:,1:3) = sos(:,1:3).*(g'); % Gain Adjustment
    y = sosfilt(makeSOS(b,a),x); 
    else
        error(['unknown filter mode: ' op.fltMode])
    end
    
end
