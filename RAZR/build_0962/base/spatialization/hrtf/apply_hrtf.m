function out_mat = apply_hrtf(in, setup, op)
% APPLY_HRTF - Convolution of an input matrix with HRIRs along columns
%
% Usage:
%   out = APPLY_HRTF(in, setup, op)
%
% Input:
%   in      Input signal matrix, containing time signals in columns
%   setup   Structure containing the following fields:
%           Required:
%               azim - Azimuth angles   (0: front, 90: left, -90: right)
%               elev - Elevation angles (0: at ears, 90: up, -90: down)
%           Optional fields:
%               filtrange_start_spls - Vector of length size(in,2) containing
%               time samples (for each column of 'in') at which the channel-wise
%               filtering starts. From this sample on, filtering is applied on
%               the doubled length of one HRIR. Default behavior: filtering on
%               the whole signal.
%   op      RAZR-options structure
%
% Output:
%   out_mat    Filtered signal matrices for all microphones
%
% See also: APPLY_SHM


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



%% input parameters

[len, numCh] = size(in);

if ~isfield(op, 'hrtf_options')
    op.hrtf_options = struct;
end

if length(setup.azim) ~= numCh || length(setup.elev) ~= numCh
    error('Wrong number of angles.');
end

%% get database informations

% dbase is passed as string for standalone usage of this function.
% In razr, get_hrtf_dbase is called earlier, and dbase is set as a field of op.
if ~isfield(op, 'hrtf_dbase')
    force_loading = true;
    op.hrtf_dbase = get_hrtf_dbase(op.hrtf_database, op.hrtf_options, force_loading);
elseif ~isfield(op.hrtf_dbase,'fs')
    force_loading = true;
    op.hrtf_dbase = get_hrtf_dbase(op.hrtf_database, op.hrtf_options, force_loading);
end

if op.fs ~= op.hrtf_dbase.fs
    error('Sampling rate of %g Hz not supported by HRTF database "%s". Must be %g Hz.', ...
        op.fs, op.hrtf_dbase.input_string, op.hrtf_dbase.fs);
end

%% do the filtering

try
    hrir = op.hrtf_dbase.pick_hrir_func(...
        setup.azim, setup.elev, op.hrtf_dbase, op.hrtf_options);
catch exc
    if strcmp(exc.identifier, 'MATLAB:UndefinedFunction')
        error(['Function "%s" not found but needed to support the HRTF', ...
            ' database "%s".\nSee the RAZR README for details.'], ...
            op.hrtf_dbase.fname_pick, op.hrtf_dbase.input_string);
    else
        error(exc.message);
    end
end

fr = get_filtranges(len, numCh, setup, op);

nb_mics = size(hrir, 3);
out_mat = zeros(nb_mics,size(in,1),numCh); %[nb_mics, time, nb_chan]

for mic = 1:nb_mics
     out_mat(mic,:,:) = in;

    switch op.hrtf_filtfunc
        case 'fftfilt'
            if op.spat_limit_filtrange
                for n = 1:numCh
                    out_mat(mic, fr(n, 1):fr(n, 2), n) = fftfilt(squeeze(hrir(:, n, mic)), out_mat(mic, fr(n, 1):fr(n, 2), n));
                end
            else
                out_mat(mic,:,:) = fftfilt(squeeze(hrir(:, :, mic)), squeeze(out_mat(mic, :, :)));
            end
        case 'filter'        
            for n = 1:numCh
                out_mat(mic, fr(n, 1):fr(n, 2), n) = filter(squeeze(hrir(:, n, mic)), 1, out_mat(mic, fr(n, 1):fr(n, 2), n));
            end
        otherwise
            error('Unknown op.hrtf_filtfunc: %s', op.hrtf_filtfunc);
    end
end
end % eof


function filtranges = get_filtranges(len, numCh, setup, op)

if ~isfield(setup, 'filtrange_start_spls') || ...
        isempty(setup.filtrange_start_spls) || ~op.spat_limit_filtrange
    filtranges = repmat([1, len], numCh, 1);
elseif length(setup.filtrange_start_spls) == numCh
    filtranges = [...
        setup.filtrange_start_spls, ...
        min(setup.filtrange_start_spls + 2*op.hrtf_dbase.len, len)];
else
    error('Length of filtrange_start_spls must match the number of input channels.');
end

end % eof
