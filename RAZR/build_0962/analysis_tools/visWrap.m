function [] = visWrap(inStruct,dRange,skip)
    % Wrapper function for sPwrVis(Multi)
    % Input:
    % inStruct - either RAZR-IR struct or RAZR-debug struct
    % (optional input):
    % dRange   - dynamicRange for plots in [dB]
    % skip     - in order to skip user input and just plot the whole IR.
    % further wishes are entered through the console
    % 20-03, ck
    %-------------------------------------------------------------------------

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


    % explanations for console dialogue
    texts{1,1} = 'S_whl';
    texts{1,2} = 'Single plot (individual normation), whole impulse response';
    texts{2,1} = 'S_dir';
    texts{2,2} = 'Single plot (individual normation), direct sound';
    texts{3,1} = 'S_ear';
    texts{3,2} = 'Single plot (individual normation), early reflections';
    texts{4,1} = 'S_lat';
    texts{4,2} = 'Single plot (individual normation), late reverberation';
    texts{5,1} = 'S_is1';
    texts{5,2} = 'Single plot (individual normation), ISM 1st Order Output';
    texts{6,1} = 'S_is2';
    texts{6,2} = 'Single plot (individual normation), ISM 2nd Order Output';
    texts{7,1} = 'S_is3';
    texts{7,2} = 'Single plot (individual normation), ISM 3st Order Output';
    texts{8,1} = 'M_rir';
    texts{8,2} = 'Cross-normed multiplot, IR parts (entire, direct, early, late)';
    texts{9,1} = 'M_ism';
    texts{9,2} = 'Cross-normed multiplot, ISM stages (direct, 1, 2, 3)';
    texts{10,1} = 'M_deb';
    texts{10,2} = 'Cross-normed multiplot, (entire, direct, ISM1/2/3, late)';
    texts{11,1} = 'S_lat_vrs';
    texts{11,2} = 'late reverberation w/ VRS';
    availPlot = {};  % listing of available plots
    if nargin<2
        dRange = 36;  % set default dynamic range for plots
    end
    if nargin<3
        skip = [];
    end
    % analyse input
    if isfield(inStruct,'ir')
        % enter debug struct analysis mode
        visMode = 'deb';
        if ~strcmp(inStruct.ir.op.spat_mode,'array')
            error('Power Visualisation Tool only available for speaker array renderings!')
        end
        availPlot = {availPlot{:} 'S_whl'};
        if inStruct.ir.op.return_rir_parts
            availPlot = {availPlot{:} 'S_dir'};
            availPlot = {availPlot{:} 'S_ear'};
            availPlot = {availPlot{:} 'S_lat'};
            availPlot = {availPlot{:} 'S_lat_vrs'};
            availPlot = {availPlot{:} 'M_rir'};
        end
        if isfield(inStruct,'ISM_out_order')
            availPlot = {availPlot{:} 'S_is1'};
            availPlot = {availPlot{:} 'S_is2'};
            availPlot = {availPlot{:} 'S_is3'};
            availPlot = {availPlot{:} 'M_ism'};   
            availPlot = {availPlot{:} 'M_deb'};        
        end
    elseif isfield(inStruct,'sig')
        % enter IR struct analysis mode
        visMode = 'ir';
        if ~strcmp(inStruct.op.spat_mode,'array')
            error('Power Visualisation Tool only available for speaker array renderings!')
        end
        availPlot = {availPlot{:} 'S_whl'};
        if isfield(inStruct.op,'return_rir_parts')
            if inStruct.op.return_rir_parts
                availPlot = {availPlot{:} 'S_dir'};
                availPlot = {availPlot{:} 'S_ear'};
                availPlot = {availPlot{:} 'S_lat'};
                availPlot = {availPlot{:} 'M_rir'};
            end
        end
    else
        error('Invalid input structure!')
    end

    if strcmp(visMode,'ir')
        inStruct.ir = inStruct; % unify data structure
    end
    
    if isempty(skip)
    % List available plots
    fprintf('Available plots:\n')
    fprintf('Handle:   Description:\n')

    for k = 1:length(availPlot)
        ix = find(contains(texts,availPlot(k)));
        prtTxt = join([availPlot(k),'   ', texts{ix,2},'\n']);
        fprintf(prtTxt{1})
    end
    end
    
    % ask user what to plot
    cont = 1;
    while cont
        if isempty(skip)
            resp = input('Enter handle for new plot, x for exit!','s');
        else
            cont = 0;
            resp = 'S_whl';
        end
        if strcmp(resp,'x')
            cont = 0;
        elseif ~any(find(contains(availPlot,resp)))||isempty(resp)
            fprintf('invalid handle!\n')
        else
            switch resp
                case 'S_whl'
                    sPwrVis(inStruct.ir.sig,inStruct.ir.op.array_pos_rec,[],dRange)
                    title('Entire IR')
                case 'S_dir'
                    sPwrVis(inStruct.ir.sig_direct,inStruct.ir.op.array_pos_rec,[],dRange)
                    title('Direct Sound')
                case 'S_ear'
                    sPwrVis(inStruct.ir.sig_early,inStruct.ir.op.array_pos_rec,[],dRange)
                    title('Early Reflections')     
                case 'S_lat'
                    sPwrVis(inStruct.ir.sig_late,inStruct.ir.op.array_pos_rec,[],dRange)
                    title('Late Reverberation')                     
                case 'S_is1'
                    sPwrVis(inStruct.ISM_out_order.ISM_order_1.sig,inStruct.ir.op.array_pos_rec,[],dRange)
                    title('1st Order ISM') 
                case 'S_is1'
                    sPwrVis(inStruct.ISM_out_order.ISM_order_2.sig,inStruct.ir.op.array_pos_rec,[],dRange)
                    title('2nd Order ISM')   
                case 'S_is3'
                    sPwrVis(inStruct.ISM_out_order.ISM_order_3.sig,inStruct.ir.op.array_pos_rec,[],dRange)
                    title('3rd Order ISM')
                case 'M_rir'
                    dMax = max(sum(inStruct.ir.sig.^2,1));
                    figure()
                    spr1 = subplot(2,2,1);
                    sPwrVisMulti(inStruct.ir.sig,inStruct.ir.op.array_pos_rec,[],dRange,dMax);
                    title('Entire IR')
                    spr2 = subplot(2,2,2);
                    sPwrVisMulti(inStruct.ir.sig_direct,inStruct.ir.op.array_pos_rec,[],dRange,dMax);
                    title('Direct Sound')
                    spr3 = subplot(2,2,3);
                    sPwrVisMulti(inStruct.ir.sig_early,inStruct.ir.op.array_pos_rec,[],dRange,dMax);
                    title('Early Reflections')
                    spr4 = subplot(2,2,4);
                    sPwrVisMulti(inStruct.ir.sig_late,inStruct.ir.op.array_pos_rec,[],dRange,dMax);
                    title('Late Reverberation')
                    hlink_rir = linkprop([spr1,spr2,spr3,spr4],{'CameraPosition','CameraUpVector'}); % link rotations of subplots
                case 'M_ism'
                    dMax = max(sum(inStruct.ir.sig.^2,1));
                    figure()
                    spi1 = subplot(2,2,1);
                    sPwrVisMulti(inStruct.ISM_out_order.ISM_order_0.sig,inStruct.ir.op.array_pos_rec,[],dRange,dMax);
                    title('Direct Sound (0th order)')
                    spi2 = subplot(2,2,2);
                    sPwrVisMulti(inStruct.ISM_out_order.ISM_order_1.sig,inStruct.ir.op.array_pos_rec,[],dRange,dMax);
                    title('ISM, 1st Order')
                    spi3 = subplot(2,2,3);
                    sPwrVisMulti(inStruct.ISM_out_order.ISM_order_2.sig,inStruct.ir.op.array_pos_rec,[],dRange,dMax);
                    title('ISM, 2nd Order')
                    spi4 = subplot(2,2,4);
                    sPwrVisMulti(inStruct.ISM_out_order.ISM_order_3.sig,inStruct.ir.op.array_pos_rec,[],dRange,dMax);
                    title('ISM, 3rd Order')
                    hlink_ism = linkprop([spi1,spi2,spi3,spi4],{'CameraPosition','CameraUpVector'}); % link rotations
                case 'M_deb'
                    dMax = max(sum(inStruct.ir.sig.^2,1));
                    figure()
                    spd1 = subplot(2,3,1);
                    sPwrVisMulti(inStruct.ir.sig,inStruct.ir.op.array_pos_rec,[],dRange,dMax);
                    title('Entire IR')
                    spd2 = subplot(2,3,2);
                    sPwrVisMulti(inStruct.ISM_out_order.ISM_order_0.sig,inStruct.ir.op.array_pos_rec,[],dRange,dMax);
                    title('Direct Sound')
                    spd3 = subplot(2,3,3);
                    sPwrVisMulti(inStruct.ISM_out_order.ISM_order_1.sig,inStruct.ir.op.array_pos_rec,[],dRange,dMax);
                    title('ISM, 1st Order')
                    spd4 = subplot(2,3,4);
                    sPwrVisMulti(inStruct.ISM_out_order.ISM_order_2.sig,inStruct.ir.op.array_pos_rec,[],dRange,dMax);
                    title('ISM, 2nd Order')    
                    spd5 = subplot(2,3,5);
                    sPwrVisMulti(inStruct.ISM_out_order.ISM_order_3.sig,inStruct.ir.op.array_pos_rec,[],dRange,dMax);
                    title('ISM, 3rd Order')
                    spd6 = subplot(2,3,6);
                    sPwrVisMulti(inStruct.ir.sig_late,inStruct.ir.op.array_pos_rec,[],[],dMax);
                    title('Late Reverberation')
                    hlink_deb = linkprop([spd1,spd2,spd3,spd4,spd5,spd6],{'CameraPosition','CameraUpVector'});
                case 'S_lat_vrs'
                    % get the reverb directions
                    twk.legacy_calc1 = inStruct.ir.op.jp_reproduction_numerics;
                    [angles, points] = get_reverb_directions(inStruct.room, inStruct.ir.op.fdn_numDownmixDirections, ...
                    inStruct.ir.op.fdn_legacy_angles, inStruct.ir.op.fdn_twist_directions, twk, false);
 
                    [spherty, DT, K] = sphericity(points);
                    DT.Points = [DT.Points(:,2) -DT.Points(:,1) DT.Points(:,3)];
                    dMax = max(sum(inStruct.ir.sig_late_mc.^2,1));
                    
                    % plot power as usual
                    figure()
                    sPwrVis(inStruct.ir.sig_late,inStruct.ir.op.array_pos_rec,[],dRange)
                    hold on
                    %plot VRS power
                    sPwrVisMulti(inStruct.ir.sig_late_mc,min(vecnorm(inStruct.ir.op.array_pos_rec'))*DT.Points,1,dRange,dMax)

                    
                    %plot reverberation polyeder
                    hp = patch('Vertices', min(vecnorm(inStruct.ir.op.array_pos_rec'))*DT.Points, 'Faces', K);
                    set(hp, 'FaceColor', [0.7 0.7 0.7],'FaceAlpha',0.2, 'EdgeColor', [0.1 0.1 0.1], 'AlignVertexCenters', 'on');
%                     daspect([1 1 1])
%                     view(3);
                    camlight
                    lighting gouraud
                    %title('Late Reverberation')        
                    
            end
        end
    end
end

