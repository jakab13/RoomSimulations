function fdn_setup = get_fdn_setup(room, op)
% GET_FDN_SETUP - Calculation of internal FDN parameters (i.e. not to be controlled by the user - in
% contrast to options).
%
% Usage:
%   [fdn_setup, room] = GET_FDN_SETUP(room, op)
%
% Output:
%   fdn_setup   Struct containing FDN parameters

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

%% signal length
if op.debugMode
    global bugOut;
    global whereAmI;
    whereAmI = 'FDNsetup';
end
if isempty(op.len)
    switch op.len_rt_measure
        case 'max'
            rt_operation = @max;
        case 'mean'
            rt_operation = @mean;
        otherwise
            error('Unknown op.len_rt_measure: %s', op.len_rt_measure);
    end
    
    fdn_setup.len = round(...
        op.len_rt_factor * rt_operation(estimate_rt(room, op.rt_estim)) * op.fs);

elseif op.len > op.len_max
    warning([sprintf('Automatically estimated RIR length higher than op.len_max!\n'), ...
        'Length will be reduced to that value. Specify op.len_max if needed.']);
    fdn_setup.len = op.len_max;

else
    fdn_setup.len = op.len;
end

%%

rec_outside = any(room.recpos < 0 | room.recpos > room.boxsize);
%src_outside = any(room.srcpos < 0 | room.srcpos > room.boxsize);

%% Delays

fdn_setup.delays = get_fdn_delays(room, op);
fdn_setup.numDel = length(fdn_setup.delays);

%% Global bandpass

% global bandpass within absorption loops:
if op.fdn_enableBP
    [fdn_setup.b_bp, fdn_setup.a_bp] = get_ism_bp(room, op);
else
    fdn_setup.b_bp = 1;
    fdn_setup.a_bp = 1;
end

%% Absorption filters

if op.fdn_enableAbsFilt
    [fdn_setup.b_abs, fdn_setup.a_abs] = getAbsfiltCoeff(room, fdn_setup, op);
else
    fdn_setup.b_abs = num2cell(ones(fdn_setup.numDel, 1)'); % HS 17.04.2020
    fdn_setup.a_abs = num2cell(ones(fdn_setup.numDel, 1)'); % HS 17.04.2020
end

%% Smearing

if op.enableSR(end)
    fdn_setup = add_sr_filtcoeffs(fdn_setup, room, op);
end

%% Feedback matrix

seed = razr_get_seed(0, op);

if ischar(op.fdn_fmatrix)
    switch op.fdn_fmatrix
        case 'eye'
            fdn_setup.fmatrix = eye(fdn_setup.numDel);
        case 'house'
            % Householder matrix (see Jot, 1997)
            fdn_setup.fmatrix = (2/fdn_setup.numDel*ones(fdn_setup.numDel) - eye(fdn_setup.numDel));
        case 'hadam'
            fdn_setup.fmatrix = 1/sqrt(fdn_setup.numDel)*hadamard(fdn_setup.numDel);
        case 'randOrth'
            fdn_setup.fmatrix = RandOrthMat3(fdn_setup.numDel, seed);
        case 'ones'
            fdn_setup.fmatrix = ones(fdn_setup.numDel);
        otherwise
            error('Feedback matrix not valid: %s', op.fdn_fmatrix);
    end
else
    fdn_setup.fmatrix = op.fdn_fmatrix;
end

%% Allpass cascade

if isoutside(room, room.recpos) || isoutside(room, room.srcpos)
    fdn_setup.do_apc = true;
else
    fdn_setup.do_apc = false;
end

if op.fdn_enable_apc && fdn_setup.do_apc
    apc_max_decaytime = 400e-3;%400e-3;
    
    [dpos, door_cntr] = doorpos(room);
    distce_src = norm(door_cntr - room.srcpos);
    distce_rec = norm(door_cntr - room.recpos);
    
    apc_decaytime = apc_max_decaytime * (1 - exp(-0.5*distce_src)) * (1 - exp(-0.5*distce_rec));
    
    [fdn_setup.b_apc, fdn_setup.a_apc] = getAllpassCascade(apc_decaytime*op.fs);
else
    fdn_setup.b_apc = 1;
    fdn_setup.a_apc = 1;
end

%% Spatialization angles

fdn_setup.spat_mode = op.spat_mode{end};

if ~strcmp(fdn_setup.spat_mode, 'diotic')
    twk.legacy_calc1 = op.jp_reproduction_numerics;
    
    [angles, points] = get_reverb_directions(room, op.fdn_numDownmixDirections, ...
        op.fdn_legacy_angles, op.fdn_twist_directions, twk, false);
    
    fdn_setup.lrscale = panscale(room, points);
    
    if op.fdn_door_idx_to_shift_angles
        [angles] = aperture_distce_atten(angles, room, ...
            op.fdn_door_idx_to_shift_angles, rec_outside, op.plot_shifted_angles);
    end
    
    fdn_setup.azim = angles(:, 1);
    fdn_setup.elev = angles(:, 2);
    fdn_setup.relpos = points;
end

%% Reflection filters

if op.fdn_blend_reflcoeffs
    fdn_setup.reflcoeff_blend = blend_reflcoeffs(room, fdn_setup, op);

    [fdn_setup.b_refl, fdn_setup.a_refl] = getReflFiltCoeff(room.freq, fdn_setup.reflcoeff_blend, ...
        fdn_setup.b_bp, fdn_setup.a_bp, op);
else
    [fdn_setup.b_refl, fdn_setup.a_refl] = getReflFiltCoeff(room.freq, room.reflcoeff, ...
        fdn_setup.b_bp, fdn_setup.a_bp, op);
    
    % coefficients for all FDN channels:
    
        fdn_ch_idx = wall2sub(fdn2wall(fdn_setup.numDel));
    if op.fdn_enableReflFilt % HS 17.04.2020
        fdn_setup.b_refl = fdn_setup.b_refl(:,fdn_ch_idx);
        fdn_setup.a_refl = fdn_setup.a_refl(:,fdn_ch_idx);
    else   
        fdn_setup.b_refl = num2cell(ones(1, fdn_setup.numDel));
        fdn_setup.a_refl = fdn_setup.b_refl;
    end
    

end
end



function reflcoeff_blend = blend_reflcoeffs(room, fdn_setup, op)

sph = get_refl_blend_setup(room, fdn_setup, op);

numPos = size(sph.pos, 1);
reflcoeff_blend = zeros(numPos, size(room.reflcoeff, 2));

for k = 1:numPos
    refl_tmp = bsxfun(@times, room.reflcoeff(sph.vbap_triangles(k, :), :), ...
        (sph.vbap_gain(k, :)').^2);
    reflcoeff_blend(k, :) = sum(refl_tmp, 1);
    
    if ~isreal(reflcoeff_blend(k, :))
        error('Blended reflection coefficient is complex!');
    end
end
end % eof
