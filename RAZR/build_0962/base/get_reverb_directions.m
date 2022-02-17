function [angles, p] = get_reverb_directions(room, num_ang, legacy_mode, twist, tweaks, do_plot)
% GET_REVERB_DIRECTIONS - Calculates angles (relative to receiver) for virtual reverb sources
%
% Usage:
%   [angles, points] = GET_REVERB_DIRECTIONS(room, num_ang, [tweaks], [do_plot])
%
% Input:
%   room            room structure (see RAZR)
%   num_ang         Number of angles being returned (order according to WALL2FDN).
%   tweaks
%   legacy_mode     If true, use legacy mode to calc reverb source positions
%   twist           If legacy_mode == true AND num_ang == 12, this flag (true or false)
%                   specifies whether the reverb sources are placed on the first or second
%                   choice of two possible diagonals on the cube surfaces.
%   do_plot         If true, create a test plot.
%
% Output:
%   angles          [azim, elev] angles in degrees
%   p               positions of reverb sources centered around [0 0 0]

%------------------------------------------------------------------------------
% RAZR engine for Mathwork's MATLAB
%
% Version 0.96.2
%
% Author(s): Torben Wendt, Josef Poppitz, Stephan D. Ewert
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



if nargin < 6
    do_plot = false;
    if nargin < 5
        tweaks = [];
        if nargin < 4
            twist = false;
            if nargin < 3
                legacy_mode = false;
            end
        end
    end
end

if isempty(tweaks)
    tweaks = struct;
end

%% Tweaks

% Some deviations in numerical accuracy to get same values as in JP's work:
default_twk.legacy_calc1 = false;  % Use when called from get_fdn_setup.m
default_twk.legacy_calc2 = false;  % Use when called from get_refl_blend_setup.m

twk = overwrite_merge(tweaks, default_twk, 1, 1);

%%
num_ang_allowed = [6, 12, 18, 24, 48, 49, 96];

if legacy_mode
    if all(num_ang ~= [12, 24])
        error('op.fdn_legacy_angles only available for 12 or 24 directions.');
    end
    % original
    % (order is in line with the standard order for walls used in all other functions,
    % especially in fdn-channel-to-wall matching):
    
    % TODO: use wall2fdn instead of hard-coded subscripts
    
    if ~twist || num_ang == 24
        p1 = zeros(12, 3);
        p1(1, :)  = [1 1 0] * 1/3;              % -z
        p1(2, :)  = [1 0 1] * 1/3;              % -y
        p1(3, :)  = [0 1 1] * 1/3;              % -x
        p1(4, :)  = [1 0 1] + [0 1 -1]*1/3;     % +x
        p1(5, :)  = [1 1 0] + [-1 0 1]*1/3;     % +y
        p1(6, :)  = [1 0 1] + [-1 1 0]*1/3;     % +z
        p1(7, :)  = [1 1 0] * 2/3;              % -z
        p1(8, :)  = [1 0 1] * 2/3;              % -y
        p1(9, :)  = [0 1 1] * 2/3;              % -x
        p1(10, :) = [1 0 1] + [0 1 -1]*2/3;     % +x
        p1(11, :) = [1 1 0] + [-1 0 1]*2/3;     % +y
        p1(12, :) = [1 0 1] + [-1 1 0]*2/3;     % +z
    else
        p1 = [];
    end
    
    % points on diagonals opposing to the first ones:
    if twist || num_ang == 24
        p2 = zeros(12, 3);
        p2(1, :)  = [1 0 0] + [-1 1 0]*1/3;     % -z
        p2(2, :)  = [1 0 0] + [-1 0 1]*1/3;     % -y
        p2(3, :)  = [0 1 0] + [0 -1 1]*1/3;     % -x
        p2(4, :)  = [1 0 0] + [0  1 1]*1/3;     % +x
        p2(5, :)  = [0 1 0] + [1  0 1]*1/3;     % +y
        p2(6, :)  = [0 0 1] + [1  1 0]*1/3;     % +z
        p2(7, :)  = [1 0 0] + [-1 1 0]*2/3;     % -z
        p2(8, :)  = [1 0 0] + [-1 0 1]*2/3;     % -y
        p2(9, :)  = [0 1 0] + [0 -1 1]*2/3;     % -x
        p2(10, :) = [1 0 0] + [0  1 1]*2/3;     % +x
        p2(11, :) = [0 1 0] + [1  0 1]*2/3;     % +y
        p2(12, :) = [0 0 1] + [1  1 0]*2/3;     % +z
    else
        p2 = [];
    end
    
    p = [p1; p2];
    p = 2*(p - 0.5);  % center points around [0 0 0]
    
else
    if twist
        one_warning('"twist" has only an effect for op.fdn_legacy_angles == true.');
    end
    % SE for 12, 18, 24
    % JP for 6, 48, 96
    switch num_ang
        case 6
            p = zeros(6, 3);
            p(1, :)  = [0 0 -1];    % -z
            p(2, :)  = [0 -1 0];    % -y
            p(3, :)  = [-1 0 0];    % -x
            
            p(4, :)  = [1 0 0];     % +x
            p(5, :)  = [0 1 0];     % +y
            p(6, :)  = [0 0 1];     % +z
            
        case 12
            % new diag mode: more omnidirectional, equally spaced to 8 directions, avoids clustering of
            % original implementation on tetrahedron
            p = zeros(12, 3);
            offset = tan(360/8/2/180*pi); %0.5;
            p(1, :)  = [offset -offset -1];    % -z
            p(2, :)  = [offset -1 offset];     % -y
            p(3, :)  = [-1 -offset offset];    % -x
            
            p(4, :)  = [1 offset offset];      % +x
            p(5, :)  = [-offset 1 offset];     % +y
            p(6, :)  = [-offset -offset 1];    % +z
            
            p(7, :)  = [-offset offset -1];    % -z
            p(8, :)  = [-offset -1 -offset];   % -y
            p(9, :)  = [-1 offset -offset];    % -x
            
            p(10, :)  = [1 -offset -offset];   % +x
            p(11, :)  = [offset 1 -offset];    % +y
            p(12, :)  = [offset offset 1];     % +z
            
        case 18
            % triangle per surface
            p = zeros(18, 3);
            % CK 20-04-03: Wäre es hier nicht viel cooler, die VRS so zu
            % drehen, dass sie dort sind, wo unsere Lautsprecher hängen?
            % Das könnte den Einfluss von VBAP Artfakten verringern!
            offsetLR = tan(360/6/2/180*pi);
            offsetUP = 0.5;
            p(1,:)  = [offsetUP -offsetLR -1];    % -z
            p(2,:)  = [offsetLR -1 offsetUP];     % -y
            p(3,:)  = [-1 -offsetLR -offsetUP];   % -x
            
            p(4,:)  = [1 offsetLR -offsetUP];     % +x
            p(5,:)  = [-offsetLR 1 offsetUP];     % +y
            p(6,:)  = [-offsetUP offsetLR 1];     % +z
            
            p(7,:)  = [offsetUP offsetLR -1];     % -z
            p(8,:)  = [-offsetLR -1 offsetUP];    % -y
            p(9,:)  = [-1 offsetLR -offsetUP];    % -x
            
            p(10,:)  = [1 -offsetLR -offsetUP];   % +x
            p(11,:)  = [offsetLR 1 offsetUP];     % +y
            p(12,:)  = [-offsetUP -offsetLR 1];   % +z
            
            p(13,:)  = [-offsetUP 0 -1];          % -z
            p(14,:)  = [0 -1 -offsetUP];          % -y
            p(15,:)  = [-1 0 offsetUP];           % -x
            
            p(16,:)  = [1 0 offsetUP];            % +x
            p(17,:)  = [0 1 -offsetUP];           % +y
            p(18,:)  = [offsetUP 0 1];            % +z
            
        case {24, 48, 96}
            p = snubcube(num_ang);
            
            if twk.legacy_calc1
                p = p/2 + 0.5 - 0.5;
            end
            
            if twk.legacy_calc2
                p = (p/2 + 0.5);
                p = p + repmat(room.recpos - [1 1 1]/2, num_ang, 1);
                p = p - repmat(room.recpos - [1 1 1]/2, num_ang, 1);
                p = (p - 0.5).*2;
            end
            
        otherwise
            error('Allowed number of directions: %s', num2str(num_ang_allowed));
    end
end

[azim, elev] = cart2sph(p(:, 1), p(:, 2), p(:, 3));

azim = rd2dg(azim);
elev = rd2dg(elev);

% account for reciever orientation:
azim = azim - room.recdir(1);
elev = elev - room.recdir(2);

angles = wrap_angles([azim, elev], 'deg');

if nargout == 0
    clear angles;
end

if ~do_plot
    return;
end

%% Plot

[spherty, DT, K] = sphericity(p);

figure;
plot3(p(:, 1), p(:, 2), p(:, 3), '*');
for k = 1:size(p,1)
    text(p(k,1),p(k,2),p(k,3)+0.1,int2str(k))
end
hold on
plotbox(gca, 2*ones(1, 3), [-1 -1 -1]);
plot3(0, 0, 0, 'k+');  % center position
xlabel('x')
ylabel('y')
zlabel('z')
set(gca, 'DataAspectRatio', [1 1 1]);
grid on;
view(-17, 16);
title(sprintf('%d reverb sources. Sphericity: %g', num_ang, spherty));

% plot
figure
hp = patch('Vertices', DT.Points, 'Faces', K);
set(hp, 'FaceColor', 'red','FaceAlpha', 0.5, 'EdgeColor', 'red', 'AlignVertexCenters', 'on');
daspect([1 1 1])
view(3);
camlight
lighting gouraud
xlabel('x')
ylabel('y')
zlabel('z')
title(sprintf('polyhedron of normalized directions, sphericity %g', spherty))
grid on;
