function [out, out_post] = spatialize(inmat, setup, op)
% SPATIALIZE - Spatial rendering of a signal matrix
%
% Usage:
%   out = SPATIALIZE(outmat, setup, op)
%
% Input:
%   inmat   Input signal matrix
%   setup   Structure whose required fields depend on spatialization mode:
%           - spat_mode:  'shm' (spherical head model), 'hrtf', 'diotic', 'ild',
%                         'array' (loudspeaker array)
%           - azim, elev: Azimuth and elevation angles for binaural rendering
%           - lrscale:    Scaling factors for ILDs
%           - relpos:     Loudspeaker positions, relative to listener, in
%                         cartesian coordinates; for loudspeaker rendering
%   op      Options structure (see RAZR)
%
% Output:
%   out     Spatialized signal


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



if ~isfield(setup, 'do_post_spat')
    setup.do_post_spat = op.enable_post_spat;
end

if ~isfield(setup, 'post_spat_mode')
    setup.post_spat_mode = op.post_spat_mode;
end

out_post = [];

switch setup.spat_mode
    case 'shm'
        [outmat_l, outmat_r] = apply_shm(inmat, setup, op);
        out = [sum(outmat_l, 2), sum(outmat_r, 2)];
        
    case 'hrtf'
        [outmat] = apply_hrtf(inmat, setup, op);
%         out = [sum(outmat_l, 2), sum(outmat_r, 2)];
        out = squeeze(sum(outmat,3))';
    case 'diotic'
        out = repmat(sum(inmat, 2), 1, 2);
        
    case 'ild'
        outmat_l = bsxfun(@times, inmat, setup.lrscale(:, 1)');
        outmat_r = bsxfun(@times, inmat, setup.lrscale(:, 2)');
        out = [sum(outmat_l, 2), sum(outmat_r, 2)];
        
    case 'maplr1'
        % full mapping to channels L/R, round scaling factors to [0, 1, 2]
        idxl = logical(round(setup.lrscale(:, 1)));
        idxr = logical(round(setup.lrscale(:, 2)));
        out = 2*[sum(inmat(:, idxl), 2), sum(inmat(:, idxr), 2)];
        
    case 'maplr2'
        % full mapping to channels L/R, round scaling factors to [0, 2]
        idxl = logical(round(setup.lrscale(:, 1)/2));
        idxr = logical(round(setup.lrscale(:, 2)/2));
        out = 2*[sum(inmat(:, idxl), 2), sum(inmat(:, idxr), 2)];
        
    case '1stOrdAmb'
        % first order ambisonics panning
        out = zeros(size(inmat, 1), 4);
        
        for ch = 1:size(inmat, 2)
            out(:, 1) = out(:, 1) + sqrt(0.5) * inmat(:, ch);
            out(:, 2) = out(:, 2) + cosd(setup.elev(ch)) * cosd(setup.azim(ch)) * inmat(:, ch);
            out(:, 3) = out(:, 3) + cosd(setup.elev(ch)) * sind(setup.azim(ch)) * inmat(:, ch);
            out(:, 4) = out(:, 4) + sind(setup.elev(ch)) * inmat(:, ch);
        end
        
    case 'array'
        % create multichannel signal for speaker array
        [out, out_post] = array_spatialization(inmat, setup, op);
        
    otherwise
        error('spatialize: Unknown spat_mode: %s', setup.spat_mode);
end
end
