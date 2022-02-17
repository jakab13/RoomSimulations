function hrir = pick_hrir_mk2(azim, elev, dbase, options)
% PICK_HRIR_MK2
% 
% See also: HRTF_PARAMS_MK2

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


% SE WARNING THIS DOES NOT WORK mk2_legacy_angle_conversion not existing
% (missing in get_default_options, FIXED here
% angle conversion -90 -> +270:

% for debugging insert this in apply_hrtf before line 67
% % SE Test interpolation
%     for kk = 0:-0.1:-2
%         azim = kk;
%         elev = 0;
%         hrir = op.hrtf_dbase.pick_hrir_func(...
%         azim, elev, op.hrtf_dbase, op.hrtf_options);
%     end

% SE interpolate or not (old code)
interpolate = 1;

% % SE triangulate and save to options if not existing
% if ~isfield(options, 'sehack') && interpolate
%    % convert to 
%    [verts(:, 1), verts(:, 2), verts(:, 3)] = sph2cart(dbase.grd.hrtf_grid_deg(:,1)/180*pi, dbase.grd.hrtf_grid_deg(:,2)/180*pi,1);
%    
%    % for convex hull
%         DT = delaunayTriangulation(verts);
%         [faces,v] = convexHull(DT); % v = volume
%     options.sehack.mesh = faces;
%     options.sehack.verts = verts;
%         
% end

% if (0)
%    figure(1)
% hold off
% % hp = patch('Vertices',verts,'Faces',mesh(:,1:3));
% % set(hp,'FaceColor','none','EdgeColor','black');
% % oc = 1;
% % oc2 = 5;
%  hp = patch('Vertices',options.sehack.verts,'Faces',options.sehack.mesh(:,1:3));
%  set(hp,'FaceColor','none','EdgeColor','black');
% % hold on
% % oc2 = 3;
% % hp = patch('Vertices',verts,'Faces',mesh(octList(1:octListLen(oc,oc2),oc,oc2),1:3));
% % set(hp,'FaceColor','none','EdgeColor','red');
% % oc2 = 4;
% % hp = patch('Vertices',verts,'Faces',mesh(octList(1:octListLen(oc,oc2),oc,oc2),1:3));
% % set(hp,'FaceColor','none','EdgeColor','blue');
% 
% daspect([1 1 1])
% view(3);
% camlight
% lighting gouraud
% xlabel('x')
% ylabel('y')
% zlabel('z')
% %title(['image sources marked by asterisks, squares, triangles, red if occluded']);
% grid on;
%     
%     
% end

if isfield(options, 'mk2_legacy_angle_conversion') && options.mk2_legacy_angle_conversion
    azim = azim - (sign(azim) - 1)*180;
else
    azim = pm180_to_360(azim);
end

num = length(azim);
hrir = zeros(dbase.len, num, 2);


if ~interpolate
    for n = 1:num
        [Azim, Elev] = GetNearestAngle(azim(n), elev(n), dbase.grd.hrtf_grid_deg);
        load(fullfile(dbase.path, GetSaveStr(Azim, Elev)));  % loads y_hrir
        hrir(:, n, :) = y_hrir;
    end
    return;
end

for n = 1:num
    %[Azim, Elev] = GetThreeNearestAngles(azim(n), elev(n), dbase.grd.hrtf_grid_deg);
    
    %[Azim, Elev] = GetTriangle(azim(n), elev(n), options.sehack, dbase.grd.hrtf_grid_deg)
    
    [Azim, Elev] = GetThreeNearestAngles2(azim(n), elev(n), dbase.grd.hrtf_grid_deg);
    
% [dirs(:, 1), dirs(:, 2), dirs(:, 3)] = sph2cart(dbase.grd.hrtf_grid_deg(:,1)/180*pi,dbase.grd.hrtf_grid_deg(:,2)/180*pi, dbase.grd.hrtf_grid_deg(:,1)*0+1);
% figure(1)
% hold off
% plot3(dirs(:, 1), dirs(:, 2), dirs(:, 3), 'k.');
% hold on
% [source(:, 1), source(:, 2), source(:, 3)] = sph2cart(azim(n)/180*pi,elev(n)/180*pi, azim(n)*0+1);
% plot3(source(:, 1), source(:, 2), source(:, 3), 'ro');
% verts = [];
% [verts(:, 1), verts(:, 2), verts(:, 3)] = sph2cart(Azim/180*pi,Elev/180*pi, Azim*0+1);
% plot3(verts(:, 1), verts(:, 2), verts(:, 3), 'g>');
% pause()    
    
    for kk = 1:3
        load(fullfile(dbase.path, GetSaveStr(Azim(kk), Elev(kk))));  % loads y_hrir
        hrirRaw(:,:,kk) = y_hrir;
    end
    
    % interpolate
    y_hrir = interpolateHRIR(azim(n), elev(n), Azim', Elev', hrirRaw);
    
%     % SE convert to vectors
%     [verts(:,1), verts(:,2), verts(:,3)] = sph2cart(Azim, Elev, 1);
%     [dir(:,1), dir(:,2), dir(:,3)] = sph2cart(azim(n), elev(n), 1);
%     y_hrir = interpolateHRIR(dir, verts, hrirRaw);
    
    hrir(:, n, :) = y_hrir;
end

end

%function hrir = interpolateHRIR(dir, verts, hrirRaw);
function hrir = interpolateHRIR(azim, elev, Azim, Elev, hrirRaw);

% OMG use awful angle stuff just to reuse the original vbap

% should be columns
% Azim = rd2dg(Azim);
% Elev = rd2dg(Elev);

ls_dirs = [Azim Elev];

%[array.ls_groups, array.layout] = findLsTriplets(ls_dirs);

layoutInvMtx = invertLsMtx(ls_dirs, [1 2 3]);

% source_azim = rd2dg(azim);
% source_elev = rd2dg(elev);
% source_angles = [source_azim, source_elev];

source_angles = [azim, elev];

% use VBAP
gains = vbap(source_angles, [1 2 3], layoutInvMtx);

% convert to lin gains, this should be barycentric coords in the triangle
gains = gains.^2;

if ~any(gains)
    hello=1
%     source_angles
%     ls_dirs
%     Triangles
    
%    pause();
% else
%     hello=0
end


% interpolate magnitude spectrum and unwrapped phase
for chan = 1:2
    
    hrtfRaw = fft(squeeze(hrirRaw(:,chan,:)));
    
    hrtfAbs = abs(hrtfRaw(:,1))*gains(1) + abs(hrtfRaw(:,2))*gains(2) + abs(hrtfRaw(:,3))*gains(3);
    hrtfAngle = unwrap(angle(hrtfRaw(:,1)))*gains(1) + unwrap(angle(hrtfRaw(:,2)))*gains(2) + unwrap(angle(hrtfRaw(:,3)))*gains(3);
    
    %hrtfAngle = (angle(hrtfRaw(:,1)))*gains(1) + (angle(hrtfRaw(:,2)))*gains(2) + (angle(hrtfRaw(:,3)))*gains(3);
    
    
    len = size(hrirRaw,1);
    mask = [1; 2*ones(len/2-1,1); 1; zeros(len/2-1,1)];
    
    hrir(:,chan) = real(ifft((hrtfAbs.*mask).*exp(i*hrtfAngle)));
end
% %hrir = ifft((abs(hrtfRaw(:,1)).*mask).*exp(i*(angle(hrtfRaw(:,3)))));
% figure(1)
% plot(hrir)
% figure(2)
% plot(abs(fft(hrir)))
% hold on
% plot(mask)
% 
% %hrir = ifft((hrtfAbs.*mask).*exp(i*unwrap(angle(hrtfRaw(:,3)))));
% %hrir = ifft((hrtfAbs.*mask).*exp(i*hrtfAngle));
% figure
% plot(20*log10(abs(hrtfRaw(:,1))))
% hold on
% plot(20*log10(abs(hrtfRaw(:,2))))
% plot(20*log10(abs(hrtfRaw(:,3))))
% plot(20*log10(hrtfAbs),'k')
% pause

end

function [Azim, Elev] = GetTriangle(az, el, surface, angles)

N = length(surface.verts);
[dir(1), dir(2), dir(3)] = sph2cart(az/180*pi, el/180*pi, 1);

intersect = 0;
for kk = 1:N
    tri(1,:) = surface.verts(surface.mesh(kk,1),:);
    tri(2,:) = surface.verts(surface.mesh(kk,2),:);
    tri(3,:) = surface.verts(surface.mesh(kk,3),:);
    
    % this intersection does not work safely
    %[intersect, t, u, v, xcoor] = TriangleRayIntersection ([0 0 0], dir, tri(1,:), tri(2,:), tri(3,:));
    
    % use vbap code
    Azim = angles(surface.mesh(kk,:),1);
    Elev = angles(surface.mesh(kk,:),2);
    ls_dirs = [Azim Elev];

%[array.ls_groups, array.layout] = findLsTriplets(ls_dirs);

layoutInvMtx = invertLsMtx(ls_dirs, [1 2 3]);

% source_azim = rd2dg(azim);
% source_elev = rd2dg(elev);
% source_angles = [source_azim, source_elev];

source_angles = [az, el];

% use VBAP, same missing detections as with triangle collision. Degenerate
% mesh?
gains = vbap(source_angles, [1 2 3], layoutInvMtx);
    
if any(gains)
     break;
    end
    
%     if intersect
%         break;
%     end
end

% use original angles for return
Azim = angles(surface.mesh(kk,:),1)';
Elev = angles(surface.mesh(kk,:),2)';
    
end

function [Azim,Elev] = GetThreeNearestAngles2(az,el,matrix)

% SE get five nearest angles forming three triangles
% we now use knowledge of the regular topology and extract a quad
% with special cases above elvation of 88 and below -64

% calculate the distance from [az el] to all angles in the matrix
N = length(matrix);
difference = repmat( [az el] , [N 1] ) - matrix;
% % SE 19.02.2020 this distance calculations in spherical coords is not
% % correct. At least the most obvious problem with the circularity of the azimuth 
% % has to be fixed
% unwrapIdx = find(abs(difference(:,1)>180));
% difference(unwrapIdx,1) = 360 - difference(unwrapIdx,1);
% % 2nd, we need to consider the elevation
% difference(:,1) = difference(:,1).*sin(difference(:,2)/180*pi);
% % -SE
% distance = sum( difference.^2 , 2 );

% SE exact squared distance d^2/r^2 =
distance = (2 - 2*cos(el/180*pi)*cos(matrix(:,2)/180*pi).*cos(difference(:,1)/180*pi)- 2 *sin(el/180*pi)*sin(matrix(:,2)/180*pi));


% get the matrix index of the angle with the smallest distance
[ans0, idx] = sort( distance );          % TW: added ans0 for compatibility, 2014-12-02


for kk = 1:N
    % start with the nearest one
    angles = matrix( idx(kk) , : );
    % avoid picking of degenerate triangles, WARNING not nice
    if el <= -64 && angles(2) ~= -64
        continue;
    elseif el >= 88 && angles(2) ~= 88
        continue;
    end
    
    break;
end
Azim(1) = angles(1);
Elev(1) = angles(2);
   
% second one
for kk = 1:N
        angles = matrix( idx(kk) , : );
        % now select the next with the same elevation, special cases for
        % >88 and < -64 are implicitly taken care of because of the first
        % vert
        if angles(1) == Azim(1)
            continue;
        end
        
        if angles(2) ~= Elev(1)
            continue;
        end
        
         break;
end
Azim(2) = angles(1);
Elev(2) = angles(2);

% third one
for kk = 1:N
    angles = matrix( idx(kk) , : );
    % now we have two with the same elevation
    if el <= -64
        % in this case we want to flip to the other side
        az2 = az + 180;
        if az2 >= 360
            az2 = az2 - 360;
        end
        difference2 = repmat( [az2 el] , [N 1] ) - matrix;
        distance2 = sum( difference2.^2 , 2 );
        [ans0, idx2] = sort( distance2 );
        angles = matrix( idx2(1) , : );
    elseif el >= 88
        % in this case we want to make sure we select elvation 90
        if angles(2) ~= 90
            continue;
        end
    else
        % now we want to make sure we a different elevation
%         if angles(1) ~= Azim(1)
%             continue;
%         end
        if angles(2) == Elev(1)
            continue;
        end
        
        % SE why was this removed?? put back 5/27/2020 3:54:48 PM
		% make sure we step in the right elevation direction
        if angles(2) > Elev(1)
            if el <= Elev(1)
                continue
            end
        else
            if el > Elev(1)
                continue
            end
        end
    end
    
     break;
end
Azim(3) = angles(1);
Elev(3) = angles(2);

% forth one
for kk = 1:N
    angles = matrix( idx(kk) , : );
    if el <= -64
        % now get another at the opposite side, the calculations must
        % have been done already
        %                 az2 = az + 180;
        %                 if az2 >= 360
        %                     az2 = az2 - 360;
        %                 end
        %                 difference2 = repmat( [az2 el] , [N 1] ) - matrix;
        %                 distance2 = sum( difference2.^2 , 2 );
        %                 [ans0, idx2] = sort( distance2 );
        angles = matrix( idx2(2) , : );
        
        % SE fixed 5/27/2020 8:40:36 AM
        % this could produce a degenerate triangle, check
        for tryit = 3:N
            if angles(1) == Azim(3)
                if abs(Azim(1) - angles(1)) == 180
                    angles = matrix( idx2(tryit) , : );
                    continue;
                end
            end
            break;
        end
        % end fix
        
    elseif el >= 88
        % get another one at elevation 88, not matching the first two
        % azimuths
        if angles(2) ~= 88
            continue;
        end
        if angles(1) == Azim(1) || angles(1) == Azim(2)
            continue;
        end
    else
        % find one matching the last elevation with another azimuth
        % SE only possible if the last elevation is not 90, fixed 5/27/2020 3:56:38 PM
        if Elev(3) == 90
            if angles(1) == Azim(3) || angles(1) == Azim(2) || angles(1) == Azim(1)
                continue;
            end
            if angles(2) ~= 88
                continue;
            end
        else
            if angles(1) == Azim(3)
                continue;
            end
            if angles(2) ~= Elev(3)
                continue;
            end
        end
    end
    
     break;
end

Azim(4) = angles(1);
Elev(4) = angles(2);

% fifth one
for kk = 1:N
    angles = matrix( idx(kk) , : );
    % now we want to have one with another elavation than the third
    
        if angles(2) == Elev(3) || angles(2) == Elev(1)
            continue;
        end
    
     break;
end
Azim(5) = angles(1);
Elev(5) = angles(2);


% supposed to be that
Triangles = [1 2 3; 1 3 4;1 2 5];

ls_dirs = [Azim' Elev'];
layoutInvMtx = invertLsMtx(ls_dirs, Triangles);

source_angles = [az, el];

gains = vbap(source_angles, Triangles, layoutInvMtx);

if any(gains)
use = find(gains>0);
if max(use) == 4
    % second
    Azim = ls_dirs(Triangles(2,:),1)';
    Elev = ls_dirs(Triangles(2,:),2)';
elseif max(use) == 5
    % second
    Azim = ls_dirs(Triangles(3,:),1)';
    Elev = ls_dirs(Triangles(3,:),2)';
else
    % first
    Azim = ls_dirs(Triangles(1,:),1)';
    Elev = ls_dirs(Triangles(1,:),2)';
end

else
   
    % something went wront, anyway return first triangle as output
		hello = 2
%     source_angles
%     ls_dirs
%     Triangles
%     %pause();
    
    Azim = ls_dirs(Triangles(1,:),1)';
    Elev = ls_dirs(Triangles(1,:),2)';

end

end

function [Azim,Elev] = GetThreeNearestAngles(az,el,matrix)

% SE get three nearest angles, 
% WARNUNG should be done with triangle intersection on a fixed mesh
% we now use knowledge of the regular topology and extract a quad

% calculate the distance from [az el] to all angles in the matrix
N = length(matrix);
difference = repmat( [az el] , [N 1] ) - matrix;
% % SE 19.02.2020 this distance calculations in spherical coords is not
% % correct. At least the most obvious problem with the circularity of the azimuth 
% % has to be fixed
% unwrapIdx = find(abs(difference(:,1)>180));
% difference(unwrapIdx,1) = 360 - difference(unwrapIdx,1);
% % 2nd, we need to consider the elevation
% difference(:,1) = difference(:,1).*sin(difference(:,2)/180*pi);
% % -SE
% distance = sum( difference.^2 , 2 );

% SE exact squared distance d^2/r^2 =
distance = (2 - 2*cos(el/180*pi)*cos(matrix(:,2)/180*pi).*cos(difference(:,1)/180*pi)- 2 *sin(el/180*pi)*sin(matrix(:,2)/180*pi));


% WE WORK WITH ANGLES, SO WE HAVE TO WRAP AZIMUTH AROUND
% N = length(matrix);
% distance2 = ones(N,1)*1000;
% for kkk = 1:N
%     
%     dist = sum(([az el] - matrix(kkk,:)).^2);
%     
%     if matrix(kkk,1) == 0 % also check for azimuth 360
%         dist = min( dist, sum(([az el] - [360 matrix(kkk,2)]).^2));
%     end
%     
%     distance2(kkk,1) = dist;
% end



% get the matrix index of the angle with the smallest distance
[ans0, idx] = sort( distance );          % TW: added ans0 for compatibility, 2014-12-02

for kk = 1:3
    angles = matrix( idx(kk) , : );
    % avoid picking of degenerate triangles, WARNING not nice
    if el <= -64 && angles(2) ~= -64
        step = 1;
            while angles(2) ~= -64
                angles = matrix( idx(kk+step) , : );
                step = step+1;
            end
    end
    
    if kk == 3
        % now check whether the first two share the same azim
        if angles(1) == Azim(2) && angles(1) == Azim(1)
            step = 1;
            while angles(1) == Azim(2)
                angles = matrix( idx(kk+step) , : );
                step = step+1;
            end
        elseif angles(2) == Elev(2) && angles(2) == Elev(1)
            if angles(2) == -64 && el <= -64
                az2 = az + 180;
                if az2 >= 360
                    az2 = az2 - 360;
                end
                difference2 = repmat( [az2 el] , [N 1] ) - matrix;
                distance2 = sum( difference2.^2 , 2 );
                [ans0, idx2] = sort( distance2 );
                angles = matrix( idx2(1) , : );
                
            else
                step = 1;
                while angles(2) == Elev(2)
                    angles = matrix( idx(kk+step) , : );
                    step = step+1;
                end
            end
        end
    end
    
    Azim(kk) = angles(1);
    Elev(kk) = angles(2);
end

end
