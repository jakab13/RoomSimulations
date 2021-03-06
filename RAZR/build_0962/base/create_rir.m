function [ir, ism_setup, fdn_setup, ism_data,op] = create_rir(room, op)
% CREATE_RIR - Synthesize (binaural) room impulse response for specified room
%
% Usage:
%   [ir, ism_setup, fdn_setup, ism_data] = CREATE_RIR(room, [op])
%
% Input:
%   room        Room structure (see help RAZR)
%   op          Options structure (see help RAZR)
%
% Output:
%	ir          IR structure (see help RAZR)
%   ism_setup   ISM setup (see GET_ISM_SETUP)
%   fdn_setup   FDN setup (see GET_FDN_SETUP)
%   ism_data    ISM metadata (see IMAGE_SOURCE_MODEL)
%
% See also: RAZR, GET_DEFAULT_OPTIONS

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



scurr = get_rng_state;
tic_total = tic;

%% input

if nargin < 2
    op = struct;
    if nargin < 1
        error('Not enough input parameters.');
    end
end

op = complement_options(op);

% SE 13.12.2018/13.02.2020 keep the original receiver coords for headphone spatialization
op.array_pos_rec = op.array_pos;

op = rotate_array(room, op); % transform from rec coords to world coords
room = complement_room(room, op);
room_flds = {'srcpos', 'srcdir', 'recpos', 'recdir'};

%% HRTF stuff

if any(strcmp(op.spat_mode, 'hrtf')) || strcmp(op.post_spat_mode, 'hrtf')
    op.hrtf_dbase = get_hrtf_dbase(...
        op.hrtf_database, op.hrtf_options, op.force_loading_hrtf);
else
    op.hrtf_dbase = struct;
end

if isfield(op.hrtf_dbase, 'fs') && op.hrtf_dbase.fs ~= op.fs
    warning(['Specified sampling rate (%g Hz) differs from that of the\n', ...
        'HRTF database "%s" (%g Hz). It has been adjusted.'], ...
        op.fs, op.hrtf_dbase.input_string, op.hrtf_dbase.fs);
    op.fs = op.hrtf_dbase.fs;
end

for n = size(room.srcpos, 1):-1:1
    room_tmp = room;
    for f = 1:length(room_flds)
        room_tmp.(room_flds{f}) = room.(room_flds{f})(n, :);
    end
    
    %% ISM
    
    tic_ism = tic;
    
    ism_setup = get_ism_setup(room_tmp, op);
    
    if op.fdn_smart_ch_mapping
        op.return_ism_data = true;
    end
    
    [ir_tmp, ism_data] = image_source_model(room_tmp, ism_setup, op);
    
    if op.verbosity >= 2
        fprintf('Elapsed time ISM:   %f seconds.\n', toc(tic_ism));
    end
    
    %% FDN
    
    if op.fdn_enabled
        tic_fdn = tic;
        fdn_setup = get_fdn_setup(room_tmp, op);
        
        [FDNinputmat, fdn_setup] = gen_fdn_input(...
            ir_tmp, ism_setup, fdn_setup, ism_data, room_tmp, op);
        
       if op.debugMode
        global bugOut;
        bugOut.FDNinputmat = FDNinputmat;
               
        end
        
        [ir_tmp.sig_late, ir_tmp.sig_late_mc, ir_tmp.sig_late_post_spat] = ...
            feedback_delay_network(FDNinputmat, fdn_setup, op);
        
        ir_tmp = struct_zeropad(ir_tmp, {...
            'sig',               'sig_direct',           'sig_early',           'sig_late', ...
            'sig_post_spat',     'sig_direct_post_spat', 'sig_early_post_spat', 'sig_late_post_spat'}, ...
            true);
        
        if ~strcmp(op.spat_mode{end}, '1stOrdAmb')
            ir_tmp.sig = ir_tmp.sig + ir_tmp.sig_late;
            if op.enable_post_spat
                ir_tmp.sig_post_spat = ir_tmp.sig_post_spat + ir_tmp.sig_late_post_spat;
            end
        else
            warning('op.spat_mode == ''1stOrdAmb'' only implemented for FDN yet.');
        end
        
        if op.verbosity >= 2
            fprintf('Elapsed time FDN:   %f seconds.\n', toc(tic_fdn));
        end
    else
        fdn_setup = struct;
    end
    
    %%
    ir_tmp.fs = op.fs;
    
    if isempty(op.rirname) && isfield(room_tmp, 'name')
        ir_tmp.name = room_tmp.name;
    else
        ir_tmp.name = op.rirname;
    end
    
    if op.return_op
        ir_tmp.op = rmfield(op, 'hrtf_dbase');
    end
    
    % restore random number generator:
    init_rng(scurr);
    
    if op.verbosity >= 1
        fprintf('Elapsed time total: %f seconds.\n', toc(tic_total));
        fprintf('Length of RIR:      %f seconds.\n', size(ir_tmp.sig, 1)/ir_tmp.fs);
    end
    
    ir_tmp.version = get_razr_version;
    
    ir(n) = ir_tmp;
    
end
if op.debugMode
    global bugOut;
    bugOut.room = room;
end
if ~op.return_rir_parts
    fld_suffixes = {'direct', 'early', 'late', ...
        'direct_post_spat', 'early_post_spat', 'late_post_spat'};
    for f = 1:length(fld_suffixes)
        if isfield(ir, ['sig_', fld_suffixes{f}])
            ir = rmfield(ir, ['sig_', fld_suffixes{f}]);
        end
    end
end

if ~op.fdn_return_mc_output && op.fdn_enabled
    ir = rmfield(ir, 'sig_late_mc');
end

ir = rmfield(ir, 'curr_sigmat');

if ~op.return_ism_sigmat
    ir = rmfield(ir, 'early_refl_sigmat');
    ir = rmfield(ir, 'early_refl_sigmat_diffuse');
end

if ~op.enable_post_spat && isfield(ir, 'sig_late_post_spat')
    ir = rmfield(ir, 'sig_late_post_spat');
end


% rm paths from hrtf database:
if any(strcmp(op.spat_mode, 'hrtf'))
    try
        if op.hrtf_dbase.did_addpath
            path(op.hrtf_dbase.oldpath);
        end
    catch
        error('Parameter "did_addpath" does not exist. It must be set in %s', ...
            op.hrtf_dbase.fname_params);
    end
end
