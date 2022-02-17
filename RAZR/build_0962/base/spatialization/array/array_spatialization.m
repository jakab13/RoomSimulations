function [out, out_post] = array_spatialization(inmat, setup, op)
% ARRAY_SPATIALIZATION - Perform spatialization for loudspeaker array
% Usage, input and output are the same as in SPATIALIZE. Actually, this function
% is called by SPATIALIZE.

%------------------------------------------------------------------------------
% RAZR engine for Mathwork's MATLAB
%
% Version 0.96.2
%
% Author(s): Torben Wendt, Josef Poppitz
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


array.src_coords = to_unit_sphere(setup.relpos);
[array.speaker_coords, array.speaker_distces] = to_unit_sphere(op.array_pos);

% for loudspeaker mapping: nearest speaker or VBAP
switch op.array_render
    case 'nearest'
        [panning_matrix, array.active_spk_idx] = ...
            map_to_speakers(array.src_coords, array.speaker_coords);
        out = render_array(...
            inmat, panning_matrix, array.speaker_distces, op, true);
        
        if setup.do_post_spat
            out_no_distce = render_array(...
                inmat, panning_matrix, array.speaker_distces, op, false);
        end
        
    case 'vbap'
        array.vbap_setup = get_vbap_setup(...
            array.speaker_coords, array.src_coords, op);
        out = apply_vbap(...
            inmat, array.vbap_setup, array.speaker_distces, op, true);
        
        if setup.do_post_spat
            out_no_distce = apply_vbap(inmat, array.vbap_setup, ...
                array.speaker_distces, op, false);
            array.active_spk_idx = unique(array.vbap_setup.speaker_idx);
        end
        
    otherwise
        error('spatialize: Unknown op.array_render: %s', op.array_render);
end

if setup.do_post_spat
    % SE 13.12.2018/13.02.2020 FIXME: the array positions in op.array_pos are in world
    % coordinates. The original positions are in receiver coordinates (x to front, y to left, z to up), and
    % for a rotated receiver (in the examples it is ussually oriented to
    % pos y direction (90 deg counterclockwise), the array_pos are rotated.) For headphone
    % spatialization we need, however, the not rotated original positions in receiver coords
    
    post_setup = get_postspat_setup(...
        setup.post_spat_mode, array.active_spk_idx, op);
    
    if op.jp_reproduction_numerics
        % SE WARNING these might be broken because the angles are now in
        % receiver space as required for SHM
        switch op.array_render
            case 'nearest'
                out_post = legacy_nearest_postspat(...
                    inmat, post_setup, panning_matrix, op);
            case 'vbap'
                out_post = legacy_vbap_postspat(...
                    inmat, post_setup, array.vbap_setup, op);
        end
    else
        out_post = spatialize(...
            out_no_distce(:, array.active_spk_idx), post_setup, op);
    end
else
    out_post = [];
end
end


function setup = get_postspat_setup(post_spat, spk_idx, op)
% Setup for recursive call of spatialize() in order to perform binaural
% post-spatialization of a multichannel signal

not_implemented = {'ild', 'maplr1', 'maplr2', 'array'};

if ismember(post_spat, not_implemented)
    error('post spatialization ''%s'' not implemented.', post_spat);
end

% [az, el] = cart2sph(...
%     op.array_pos(spk_idx, 1), ...
%     op.array_pos(spk_idx, 2), ...
%     op.array_pos(spk_idx, 3));

% SE 13.12.2018/13.02.2020 here we need the rec space coords
[az, el] = cart2sph(...
    op.array_pos_rec(spk_idx, 1), ...
    op.array_pos_rec(spk_idx, 2), ...
    op.array_pos_rec(spk_idx, 3));

setup.azim = rd2dg(az);
setup.elev = rd2dg(el);

setup.spat_mode = post_spat;
setup.do_post_spat = false;
setup.filtrange_start_spls = [];

% ILDs are not implemented:
setup.lrscale = [];

end


function out_post = legacy_nearest_postspat(inmat, post_setup, panning_matrix, op)
% JP's original code snippet

% However, I think the sorting introduced by find() is wrong in the
% original code:
[~, col] = find(panning_matrix);
[az, el] = cart2sph(op.array_pos(col, 1), op.array_pos(col, 2), op.array_pos(col, 3));
post_setup.azim = rd2dg(az);
post_setup.elev = rd2dg(el);
out_post = spatialize(inmat, post_setup, op);

end


function out = legacy_vbap_postspat(inmat, setup, vbap_setup, op)
% JP's original code (up to minor changes)

if ~strcmp(setup.spat_mode, 'hrtf')
    error(['If op.jp_reproduction_numerics == true, then ', ...
        'op.post_spat_mode must be ''hrtf''.']);
end

[numDirctns, numGainComponents] = size(vbap_setup.gains_2d);

for gainComponent = 1:numGainComponents
    [az, el] = cart2sph(...
        op.array_pos(vbap_setup.triangles(:, gainComponent), 1), ...
        op.array_pos(vbap_setup.triangles(:, gainComponent), 2), ...
        op.array_pos(vbap_setup.triangles(:, gainComponent), 3));
    
    setup.azim = rd2dg(az);
    setup.elev = rd2dg(el);
    
    [out_tmp_l, out_tmp_r] = apply_hrtf(inmat, setup, op);
    
    out_tmp_l2 = bsxfun(@times, out_tmp_l, vbap_setup.gains_2d(:, gainComponent)');
    out_tmp_r2 = bsxfun(@times, out_tmp_r, vbap_setup.gains_2d(:, gainComponent)');
    
    if op.vbap_color_compen
        for dirctn = 1:numDirctns
            out_tmp_l3(:, gainComponent, dirctn) = fftfilt(...
                vbap_setup.h_filt(:, gainComponent, dirctn), out_tmp_l2(:, dirctn));
            out_tmp_r3(:, gainComponent, dirctn) = fftfilt(...
                vbap_setup.h_filt(:, gainComponent, dirctn), out_tmp_r2(:, dirctn));
        end
    else
        for dirctn = 1:numDirctns
            out_tmp_l3(:, gainComponent, dirctn) = fftfilt(1, out_tmp_l2(:, dirctn));
            out_tmp_r3(:, gainComponent, dirctn) = fftfilt(1, out_tmp_r2(:, dirctn));
        end
    end
end

out_l = sum(out_tmp_l3,3);
out_r = sum(out_tmp_r3,3);

out = [sum(out_l, 2), sum(out_r, 2)];

end
