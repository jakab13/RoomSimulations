function [spherty, DT, K] = sphericity(points)
% SPHERICITY - For a set of points in the three-dimensional space, their
% sphericity is calculated, that is, the ratio of the surface area of a sphere
% (having the same volume as created by the set of points) to the surface area
% assigned to the set of points.
%
% Usage:
%   [spherty, DT, K] = SPHERICITY(points)
%
% Input:
%   points      Matrix of points [x1, y1, z1; ... xN, yN, zN]
%
% Output:
%   spherty     Sphericity
%   DT          Delaunay triangulation
%   K           Indices of the vertices on the convex hull
%   
% See also: DELAUNAYTRIANGULATION, CONVEXHULL

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


pn = points./repmat(sqrt(sum(points.^2, 2)), 1, 3); % normalize to unit length
DT = delaunayTriangulation(pn);
[K, vol] = convexHull(DT); % v = volume

% calculate volume
verts = DT.Points;
faces = K;

a = verts(faces(:, 2), :) - verts(faces(:, 1), :);
b = verts(faces(:, 3), :) - verts(faces(:, 1), :);
c = cross(a, b, 2);

area = 1/2 * sum(sqrt(sum(c.^2, 2)));

% sphericity: surface area of sphere with same volume in relation to surface area
sphereSurfaceArea = (pi*(6*vol)^2)^(1/3);
spherty = sphereSurfaceArea/area;
