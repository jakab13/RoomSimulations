function signalmat = apply_source_dir(database,isd,ism_setup,n,signalmat,op,AvgFilter)


%------------------------------------------------------------------------------
% RAZR engine for Mathwork's MATLAB
%
% Version 0.96.2
%
% Author(s): Henning Steffens
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


% apply_source_dir - Convolution of an input matrix with source DIRs

% choose database and set db-specific parameters

dbase = lower(database);
cfg = get_razr_cfg;
fldname = sprintf('dir_path__%s', dbase);


if ~isfield(cfg, fldname)
   error('Path to DIR databse "%s" not specified in razr_cfg.m.', database);
end

switch dbase
  
    case 'hs_filter'
        %%
        
        fs_database = 44100;
        dir_len = 4096;
        flanklen = round(1e-3*fs_database);
        win = hannwin(flanklen*2);
        flank = repmat(win((flanklen + 1):end), 1, 1);
        start_sample = 1;
        
        %%% Avg Filter
        if(AvgFilter)
            h1 = load('ir_avg_hs_filter.mat','ir_avg_hs_filter');
            
            impulse_response =  1 * h1.ir_avg_hs_filter( start_sample:dir_len +  start_sample-1,:); 
            impulse_response((dir_len - flanklen + 1):dir_len, :) = ...
            impulse_response((dir_len - flanklen + 1):dir_len, :) .* flank;
           
            dir_signalmat = impulse_response;
            conv_signalmat = fftfilt(dir_signalmat, signalmat(:,n));
            signalmat(:,n) = conv_signalmat;
            return;
        end
        
        if (isd.order(1) >= op.ism_order_avg)
                    
             h1 = load('ir_avg_hs_filter.mat','ir_avg_hs_filter');
            
            impulse_response =  1 * h1.ir_avg_hs_filter( start_sample:dir_len +  start_sample-1,:); 
            impulse_response((dir_len - flanklen + 1):dir_len, :) = ...
            impulse_response((dir_len - flanklen + 1):dir_len, :) .* flank;
           
            dir_signalmat = impulse_response;
            conv_signalmat = fftfilt(dir_signalmat, signalmat(:,n));
            signalmat(:,n) = conv_signalmat;
            return;
            
        else
         
        src_azim = isd.src_azim(n);
        
        % SE this code would be needed if we have an azimuth and elev
        src_elev = isd.src_elev(n);
         
        %elev -> approx. calc. % SE after we just converted to degree now
        %back gain :)
        src_az_rad = dg2rd(src_azim);
        src_el_rad = dg2rd(src_elev);
        
        % SE comment, all the following code was to convert azimuth and elev
        % just to one angle re front pointing direction, because hs filter
        % only needs one angle. Here the quick version, it is just the
        % projection on x
        
        [src_dirVec_x, ~, ~] = sph2cart(src_az_rad,src_el_rad,1); % FIXME we are converting back an forth from axis to angles
        
        % angle between vectors 
        % FIXME we could have directly used
        % rel_dir_dot from scale_is_pattern, which we now compute again
        % after sph2cart and cart2sph and two deg2rad and rad2deg.
        
        src_azim = 180*acos(src_dirVec_x)/pi; % yes well again in degrees. 
        src_elev = 0; % just zero this for fun. It is not just anyway in src filter.
        % end SE
        %%%%%%%%%%%%%%%%%%%%%%%%

        bUseOldFilter = 0;
        
        if(bUseOldFilter)
            outmat = hsfilter_src(src_azim, src_elev, 0,op.fs, signalmat(:,n),1);
        else
            outmat = hsfilter_src(src_azim, src_elev, 0,op.fs, op.dir_src_radius, op.dir_theta0, signalmat(:,n),1);
        end
        
        signalmat(:,n) = outmat;
        
        return;
        
        end

      
  
    otherwise
        script_name_params = sprintf('dir_params_%s', database);
        script_name_pick   = sprintf('pick_dir_%s', database);
        if exist(script_name_params, 'file')
            eval(script_name_params);
        else
            error('Script "%s.m" needed for DIR database "%s" but not found. See README.txt for details.', ...
                script_name_params, database);
        end
        if ~exist(script_name_pick, 'file')
            error('Script "%s.m" needed for DIR database "%s" but not found. See README.txt for details.', ...
                script_name_pick, database);
        end
end


end
