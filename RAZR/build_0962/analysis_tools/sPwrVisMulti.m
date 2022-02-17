function [] = sPwrVisMulti(inSig, lsPos, VRSangles, dRange, dMax)
% sPowerVis - Plot the distribution of power amongst loudspeakers in UOL
% VR-LAB (86 loudspeakers, otherwise other LS.positions have to be loaded). 
% Optionally, plot VRS-Positions.
%
% Usage:
%   sPowerVis(inSig, lsPos)
%   sPowerVis(inSig, lsPos,angles)
%   sPowerVis(inSig, lsPos,angles, dRange)
%   sPowerVis(inSig, lsPos,angles, dRange,dMax)
%   sPowerVis(inSig, lsPos,[],dRange)
%   
%
% Input:
%   inSig           Loudspeaker signal matrix (M x N)
%                   N: number of loudspeakers, M = number of samples
% Optional input:
%   lsPos           Loudspeaker positions as carthesian coordinates
%                   Nx3 matrix, N: number of loudspeakers
%   angles          Angles for VRS-visualisation (Az, El, in Deg.)
%   dRange          Dynamic Range for visualisation [dB], default: 30
%   dMax            Maximum Power of entire IR when visualising only parts
%                   of the IR in multiple plots, comparing them w/ eachother. 
%
% 2020-03 Christoph Kirsch

%%%%%%%%%% TO DO:
% - Implement for distance compensation options (dependent on global RAZR
%   dist. comp?), so that visualisation on a sphere regardless of LS
%   positions can be achieved with power compensated for distance
% - Plot VRS Positions
%%%%%%%%%%

%------------------------------------------------------------------------------
% RAZR engine for Mathwork's MATLAB
%
% Version 0.96.2
%
% Author(s): Christoph Kirsch
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




if size(inSig,1)<size(inSig,2)
    inSig = inSig';
%     warning('Input signal array invalid?')
end

    if isempty(VRSangles)
        VRSangles = [];
    end
    if isempty(dRange)
        dRange = 30;
    end
    if isempty(dMax)
        dMax = [];
    end

    if dRange > 99 || dRange <= 0
        error('Dynamic Range allowed between 1 and 99 dB')
    end    
    % calculating level per channel, offset for dotsize calculation
    inPwr = sum(inSig.^2,1);
    if isempty(dMax)
        dMax = max(inPwr);
    end
    rangeOffset = 100;
    pwrVec = 10*log10(inPwr./(dMax))+rangeOffset;
    pwrVec(isinf(pwrVec))=0;  % handling silent channels
    
    nCmap = 128; % resolution of color map
    cMap = jet(nCmap);
    
    % mapping power numbers into ranges relevant for plotting
    a = rangeOffset-dRange;     % lower bound for visualisation
    b = 3;                      % minimum ball size
    c = 100;                    % upper limit of mapping input
    d = 25;                     % maximum ball size
    pwrVecSize = ((d-b).*pwrVec+(b*c-a*d))./(c-a);  % size mapping
    pwrVecColor = round(((nCmap-1).*pwrVec+(1*c-a*nCmap))./(c-a)); % color mapping
    pwrVecSize(pwrVecSize<=0) = b;   % setting minum ball size
    pwrVecColor(pwrVecColor<=0) = 1; % avoiding errors
%     load('ls_array_vrroom_uol.mat')  % 

    sphX = lsPos(:,1);  % no mapping of directions onto a sphere
    sphY = lsPos(:,2);
    sphZ = lsPos(:,3);
%     [az, el, ~] = cart2sph(lsPos(:,1),lsPos(:,2),lsPos(:,3));
%     distSph = 2.5*ones(length(array_pos),1); % hardcode: spherical visualisation in 2.5 meters distance
%     [sphX,sphY,sphZ] = sph2cart(az,el,distSph);

%     figure()
    hold on
    if VRSangles == 1
        vrsflag = 1;
    elseif ~isempty(VRSangles)
        distSphV = 2.5*ones(1,length(VRSangles)); %%%%% HARDCODE FIX ME
        [vrsX,vrsY,vrsZ] = sph2cart(deg2rad(VRSangles(:,1)),deg2rad(VRSangles(:,2)),distSphV);
        plot3(vrsX,vrsY,vrsZ,'xb')
        vrsflag = 0;
    else
        vrsflag = 0;
    end
    if vrsflag
        marker = 'd';
    else
        marker = 'o';
    end
    for k = 1:length(sphX)
        if pwrVecSize(k)>b
            plot3(sphX(k),sphY(k),sphZ(k),marker,'MarkerSize',pwrVecSize(k),'Color','r','MarkerFaceColor',cMap(pwrVecColor(k),:))
        else
            plot3(sphX(k),sphY(k),sphZ(k),marker,'MarkerSize',pwrVecSize(k),'Color','k','MarkerFaceColor',[1 1 1])
        end
    end

    % visualising receiver within array, inspired by 'scene'-function from
    % razr/analysis_tools
            recpos = [0 0 0];

            plot3(recpos(1), recpos(2), recpos(3), 'o', ...
                'Linewidth', 2, 'markersize', 8, ...
                'color', 'k', 'markerfacecolor', 'k');

                noselen = 0.4;
                [rx, ry, rz] = sph2cart(deg2rad(0), deg2rad(0), noselen);
                nose = [recpos; recpos + [rx, ry, rz]];

                plot3(nose(:, 1), nose(:, 2), nose(:, 3), ...
                    '-', 'Linewidth', 3, 'color', 'k');
    axLimits = [-2.75 2.75];
    xlim(axLimits)
    ylim(axLimits)
    zlim(axLimits)
    view([-95 50])
    axis equal
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    colormap jet
    tTicks = linspace(0,-dRange,6);
    tLabels= {};
    for k = 1:6
        tLabels = {tLabels{:} int2str(tTicks(end-k+1))};
    end
    
    colorbar('Ticks',linspace(0,1,6),...
         'TickLabels',tLabels)
end
