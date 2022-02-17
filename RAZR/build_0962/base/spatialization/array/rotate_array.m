function op = rotate_array(room, op)
% Rotate array to account for receiver orientation

%------------------------------------------------------------------------------
% RAZR engine for Mathwork's MATLAB
%
% Version 0.96.2
%
% Author(s): Josef Poppitz
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


if isempty(op.array_pos)
    return;
end

[az, el, r] = cart2sph(op.array_pos(:, 1), op.array_pos(:, 2), op.array_pos(:, 3));
azim = rd2dg(az);
elev = rd2dg(el);

azim2 = azim + room.recdir(1);
azim2(azim2 > 180) = azim2(azim2 > 180) - 360;

elev2 = elev + room.recdir(2);
elev2(elev2 > 180) = elev2(elev2 > 180) - 360;

az2 = dg2rd(azim2);
el2 = dg2rd(elev2);

[op.array_pos(:, 1), op.array_pos(:, 2), op.array_pos(:, 3)] = sph2cart(az2, el2, r);

end
