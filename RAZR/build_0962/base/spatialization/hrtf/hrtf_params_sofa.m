function dbase = hrtf_params_sofa(dbase, cfg)
% HRTF_PARAMS_SOFA

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


%% Paths and filenames

if isfield(cfg, 'sofa_api_path')
    oldpath = addpath(cfg.sofa_api_path);
else
    error('Path to SOFA API not specified in razr.cfg.m.');
end

SOFAstart('silent');
dbase.sofaobj = SOFAload(dbase.sofa_filename);

path(oldpath);
dbase.did_addpath = false;

dbase.len = size(dbase.sofaobj.Data.IR, 3);
dbase.fs = dbase.sofaobj.Data.SamplingRate;
