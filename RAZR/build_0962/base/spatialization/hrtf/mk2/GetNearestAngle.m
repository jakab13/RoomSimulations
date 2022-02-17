function [Azim,Elev] = GetNearestAngle(az,el,matrix)

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

% calculate the distance from [az el] to all angles in the matrix
N = length(matrix);
difference = repmat( [az el] , [N 1] ) - matrix;
distance = sum( difference.^2 , 2 );

% get the matrix index of the angle with the smallest distance
[ans0, idx] = min( distance );          % TW: added ans0 for compatibility, 2014-12-02
angles = matrix( idx , : );
Azim = angles(1);
Elev = angles(2);


end
