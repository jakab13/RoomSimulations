function setup = get_vbap_setup(speaker_coordinates, source_coordinates, op)
% GET_VBAP_SETUP
%
% Usage:
%   setup = GET_VBAP_SETUP(speaker_coords, source_coords, op)
%
% Input:
%   speaker_coords  Cartesian coordinates of loudspeakers in array, relative to
%                   receiver
%   source_coords   Cartesian coordinates of virtual sources relative to receiver
%   op              Options structure (see RAZR)
%
% Output:
%   setup           Setup structure
%
% See also: APPLY_VBAP

%------------------------------------------------------------------------------
% RAZR engine for Mathwork's MATLAB
%
% Version 0.96.2
%
% Author(s): Josef Poppitz, Torben Wendt
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



[azim, elev] = cart2sph(...
    speaker_coordinates(:, 1), ...
    speaker_coordinates(:, 2), ...
    speaker_coordinates(:, 3));

azim = rd2dg(azim);
elev = rd2dg(elev);

ls_dirs = [azim'; elev']';

[array.ls_groups, array.layout] = findLsTriplets(ls_dirs);

layoutInvMtx = invertLsMtx(ls_dirs, array.ls_groups);

[source_azim, source_elev] = cart2sph(...
    source_coordinates(:, 1), ...
    source_coordinates(:, 2), ...
    source_coordinates(:, 3));

source_azim = rd2dg(source_azim);
source_elev = rd2dg(source_elev);
source_angles = [source_azim, source_elev];
setup.gain_all_speakers = vbap(source_angles, array.ls_groups, layoutInvMtx);

[row, col] = find(setup.gain_all_speakers);

q = [row' col'];
q2 = sortrows(q);

nK = size(array.ls_groups, 1);
nSC = size(speaker_coordinates,1);

K_tetrahedron = [array.ls_groups, ones(nK,1)*(nSC+1)];
array_pos_tetrahedron = [speaker_coordinates; 0 0 0];

ti = tsearchn(array_pos_tetrahedron,K_tetrahedron,source_coordinates*0.5);

setup.triangles = array.ls_groups(ti, :);

setup.gains_2d = zeros(size(source_coordinates, 1), size(array.ls_groups, 2));

for kk = 1:size(source_coordinates,1)
    setup.gains_2d(kk, :) = setup.gain_all_speakers(kk, setup.triangles(kk, :));
end

if op.vbap_color_compen
    setup.h_filt = get_coloration_compensation(setup, op.fs);
end

[numDir, numGainComp] = size(setup.triangles);

% new fields:
% TODO: remove those above that are not needed any longer,
% and provide the new ones directly in the new format
setup.gains = setup.gains_2d(:)';
setup.speaker_idx = setup.triangles(:)';
input_idx = repmat(1:numDir, numGainComp, 1)';
setup.speaker_input_idx = input_idx(:)';

if op.vbap_color_compen
    hlen = size(setup.h_filt, 1);
    setup.colcomp_ir = reshape(permute(setup.h_filt, [1 3 2]), hlen, numGainComp*numDir);
else
    setup.colcomp_ir = [];
end

if op.plot_vbap_arrangement
    plot_vbap_arrangement(speaker_coordinates, source_coordinates, setup, array);
end

end


function h_filt = get_coloration_compensation(stp, fs)
% VBAP filters with frequency-dependent P-value

Lfilt = 512;
f = (0:Lfilt/2)*fs/Lfilt;
DTT = 1; % 1 for anechoic conditions, ~0.5 for listening rooms, 0 for standard power normalization
pValue = getPValueResponse(f, DTT);
[ch_num, ls_num] = size(stp.triangles);

h_filt = zeros(Lfilt, ls_num, ch_num);
H_filt = zeros(Lfilt/2 + 1, ls_num);

for ch_idx = 1:ch_num
    for nf = 1:Lfilt/2+1
        pv_f = pValue(nf); % p-value for this frequency
        
        H_filt(nf, :) = ...
            stp.gain_all_speakers(ch_idx, stp.triangles(ch_idx, :)) ./ ...
            (...
            (sum(stp.gain_all_speakers(ch_idx, stp.triangles(ch_idx, :)).^pv_f)).^(1/pv_f) * ...
            ones(1, ls_num)...
            );
    end
    
    %H_filt_all(:, :, ch_idx) = H_filt;
    h_filt(:, :, ch_idx) = fftshift(ifft([H_filt; H_filt(end-1:-1:2, :)]), 1);
end
end


function plot_vbap_arrangement(speaker_coordinates, source_coordinates, vbap_setup, array)

figure;
plot3(speaker_coordinates(:,1), speaker_coordinates(:,2), speaker_coordinates(:,3), 'rx','LineWidth',2)
hold on;
trisurf(array.ls_groups, ...
    speaker_coordinates(:,1), speaker_coordinates(:,2), speaker_coordinates(:,3), ...
    'FaceAlpha',0)

plot3(source_coordinates(:,1), source_coordinates(:,2), source_coordinates(:,3), 'gx','LineWidth',2);
plot3(...
    speaker_coordinates(vbap_setup.triangles, 1), ...
    speaker_coordinates(vbap_setup.triangles, 2), ...
    speaker_coordinates(vbap_setup.triangles, 3), 'bx','LineWidth',2,'MarkerSize',15);
plot3(0,0,0,'mx');
%                     plot3(op.array_pos(:,1),op.array_pos(:,2),op.array_pos(:,3),'rx','LineWidth',2)
origin = repmat([0],size(source_coordinates,1),1);
quiver3(origin,origin,origin,source_coordinates(:,1),source_coordinates(:,2),source_coordinates(:,3),1.2);
axis([-1.2 1.2 -1.2 1.2 -1.2 1.2])

end
