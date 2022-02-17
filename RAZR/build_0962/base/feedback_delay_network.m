function [out, outmat, out_post] = feedback_delay_network(insig, fdn_setup, op)
% FEEDBACK_DELAY_NETWORK
%
% Usage:
%   [out, outmat] = FEEDBACK_DELAY_NETWORK(in, fdn_setup, op)
%
% Input:
%   in          Input signal (1 or length(m) channels, i.e. columns)
%   fdn_setup   Structure containg FDN setup specifications returned by GET_FDN_SETUP.
%   op          Options struct (complete, i.e. custom options already merged with defaults!)
%
% Output:
%   out         Output signal
%   outmat      Output signal (monaural, multichannel)

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


% This file will be released as a p file.

%% allpass cascade

if op.fdn_enable_apc
    % fdn_setup.max_filtrange_input was written in gen_fdn_input:
    if isfield(fdn_setup, 'max_filtrange_input')
        rnge = fdn_setup.max_filtrange_input(1):fdn_setup.max_filtrange_input(2);
    else
        rnge = 1:size(insig, 1);
    end
    insig(rnge, :) = filter(fdn_setup.b_apc, fdn_setup.a_apc, insig(rnge, :));
end

%% adjust matrix size

[len, numCh] = size(insig);

if ~(numCh == 1 || numCh == fdn_setup.numDel)
    error('Input signal must have either 1 column or length(m) columns.');
elseif numCh == 1
    insig = repmat(insig, 1, fdn_setup.numDel);
end

if len < fdn_setup.len
    insig = [insig; zeros(fdn_setup.len - len, fdn_setup.numDel)];
else
    insig = insig(1:fdn_setup.len, :);
end

if isfield(fdn_setup, 'num_pre_zeros')  % set in gen_fdn_input
    insig = [insig; zeros(fdn_setup.num_pre_zeros, fdn_setup.numDel)];
end

%%

if any(size(fdn_setup.fmatrix) ~= fdn_setup.numDel)
    error('Feedback matrix must be quadratic.')
end

%% main loop

% try
%     [outmat] = fdn_mainloop(...
%         insig, fdn_setup.delays, fdn_setup.b_abs, fdn_setup.a_abs, fdn_setup.fmatrix);
% catch
%     warning(['Mex file ''fdn_mainloop'' either not found or not executable. ', ...
%         'fdn_mainloop_m is called instead.']);
%     [outmat] = fdn_mainloop_m(...
%         insig, fdn_setup.delays, fdn_setup.b_abs, fdn_setup.a_abs, fdn_setup.fmatrix);
% end

% ACHTUNG CK 20-03-30
% run SE's "fast" M-file FDN instead of .mex in order to be able to run
% cascaded filters

if (strcmp(op.filtCreatMeth,'yw')||strcmp(op.filtCreatMeth,'cs')) && strcmp(op.fltMode,'sos')
    warning('Warning, this combination of CreatMeth and fltMode is known to be buggy, consider using "casc" instead!')
end

    [outmat] = fdn_mainloop_m_se(...
        insig, fdn_setup.delays, fdn_setup.b_abs', fdn_setup.a_abs', fdn_setup.fmatrix);
    
    % ----------------------------------- end CK 20-03-30
%% pre-trunc

if isfield(fdn_setup, 'num_pre_zeros')  % was set in gen_fdn_input
    outmat = outmat((fdn_setup.num_pre_zeros + 1):end, :);
end

%% reflection filters and downmixing

if op.fdn_reflBeforeDownmix
    outmat = apply_reflfilt(outmat, fdn_setup, op);
end

outmat = downmix_fdn(outmat, op.fdn_numDownmixDirections);

if ~op.fdn_reflBeforeDownmix
    outmat = apply_reflfilt(outmat, fdn_setup, op);
end

%% spatialization
[out, out_post] = spatialize(outmat, fdn_setup, op);

%% security check
if any(any(isnan(out)))
    error('nan detected in FDN output.');
end
end


function sigmat = apply_reflfilt(sigmat, fdn_setup, op)

if ~op.fdn_enableReflFilt
    return
end

numCh = size(sigmat, 2);

if any([length(fdn_setup.b_refl), length(fdn_setup.a_refl)] ~= numCh)
    error('Number of reflection filters must equal number of signal channels.');
end

for ch = 1:numCh
    sigmat(:, ch) = selFlt(fdn_setup.b_refl(:,ch), fdn_setup.a_refl(:,ch), sigmat(:, ch),[],op);
end
end


function sigmat = downmix_fdn(sigmat, numCh_mix)
% Channel downmixing from 96 FDN-Channels to number of directions 
% specified in numCh_mix

numCh_in = size(sigmat, 2);

if numCh_mix == numCh_in
    return
elseif numCh_mix > numCh_in
    error(['Reverb downmix: Number of downmix channels too high. ', ...
        'Must be smaller than number of delays.']);
end

numCh_allowed = [6, 12, 24, 48, 96];  % both, mix and input

if all(numCh_mix ~= numCh_allowed)
    error('Number of downmix channels must be one of [%s]', num2str(numCh_allowed));
end

if all(numCh_in ~= numCh_allowed)
    error('Number of channels must be one of [%s]', num2str(numCh_allowed));
end

if numCh_in == 96 && numCh_mix <= 48
    sigmat = [...
        sigmat(:,  1:24) + sigmat(:, 25:48), ...
        sigmat(:, 49:72) + sigmat(:, 73:96)];
end

if numCh_in >= 48 && numCh_mix <= 24
    sigmat = sigmat(:, 1:24) + sigmat(:, 25:48);
end

if numCh_in >= 24 && numCh_mix <= 12
    sigmat = sigmat(:, 1:12) + sigmat(:, 13:24);
end

if numCh_in >= 12 && numCh_mix == 6
    sigmat = sigmat(:, 1:6) + sigmat(:, 7:12);     
end
end  % eof
