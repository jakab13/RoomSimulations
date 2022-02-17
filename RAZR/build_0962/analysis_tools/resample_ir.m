function ir = resample_ir(ir, new_fs)
% RESAMPLE_IR - Resample impulse response
%
% Usage:
%   ir = RESAMPLE_IR(ir, new_fs)
%
% Input:
%   ir      RIR structure (see RAZR)
%   new_fs  Desired sampling rate in Hz
%
% Output:
%   ir      RIR structure with resampled time signal and updated field "fs"

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


if all([ir.fs] == new_fs)
    return;
end

for n = 1:length(ir)
    
    old_fs = ir(n).fs;
    
    ir(n).sig  = resample(ir(n).sig, new_fs, old_fs);
    ir(n).fs   = new_fs;
    
    if isfield(ir(n), 'len')
        ir(n).len  = size(ir(n).sig, 1);
    end
    
    if isfield(ir(n), 'start_spl')
        ir(n).start_spl = round(ir(n).start_spl*new_fs/old_fs);
    end
end
