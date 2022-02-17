function op = complement_options(options)
% COMPLEMENT_OPTIONS - Merge input options with default options and set some
% derived options.
%
% Usage:
%   op = COMPLEMENT_OPTIONS(options)

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



op_default = get_default_options;

%% Spatialization
if ~iscell(op_default.spat_mode)
    op_default.spat_mode = {op_default.spat_mode};
end

if isfield(options, 'spat_mode') && ~iscell(options.spat_mode)
    options.spat_mode = {options.spat_mode};
end

%%
if isfield(options, 'ism_only') && all(options.ism_only ~= -1)
    op_default.verbosity = 3;
end

op = overwrite_merge(options, op_default, 1, 1);

if isfield(options, 'ism_only') && all(options.ism_only ~= -1)
    op.ism_order = op.ism_only;
    op.fdn_enabled = false;
end

if isempty(op.fdn_numDownmixDirections)
    op.fdn_numDownmixDirections = op.fdn_numDelays;
end

if op.jp_reproduction_numerics
    op.array_TCelsius = 20;
end

%% Check, if there are fields in input options struct that are not in the default options:
invalid_fields = setdiff(fieldnames(op), fieldnames(op_default));

if ~isempty(invalid_fields)
    [op, invalid_fields] = restore_legacy_opts(op, invalid_fields);
    if ~isempty(invalid_fields)
        plural = length(invalid_fields) > 1;
        plural_str = repmat('s', 1, plural);
        warning('Unknown option%s: %s.', plural_str, strjoin(invalid_fields, ', '));
    end
end

end % eof


function [op, invalid_fields] = restore_legacy_opts(op, invalid_fields)
% If legacy options have been passed to razr, here they are converted to current
% equivalents.

% Sometimes an option has simply been renamed:
renaming = @(oldfield) (op.(oldfield));

% ---- Conversion rules old options --> new options: --------
new_name.tlen = 'len';
newvalue.tlen = @(oldfield) (round(op.(oldfield)*op.fs));

new_name.tlen_max = 'len_max';
newvalue.tlen_max = @(oldfield) (round(op.(oldfield)*op.fs));

new_name.ism_randFactorsInCart = 'ism_jitter_type';
newvalue.ism_randFactorsInCart = ...
    @(oldfield) ([...  % choose string depending on boolean:
    repmat('cart_legacy', 1, sign(double( op.(oldfield)))), ...
    repmat('sph',         1, sign(double(~op.(oldfield))))]);

new_name.ism_ISposRandFactor = 'ism_jitter_factor';
newvalue.ism_ISposRandFactor = renaming;

new_name.ism_enable_diffusion = 'ism_enable_scattering';
newvalue.ism_enable_diffusion = renaming;

new_name.fdn_mfactor = 'fdn_delay_stretch';
newvalue.fdn_mfactor = renaming;

new_name.fdn_hrtf_boxdiags = 'fdn_twist_directions';
newvalue.fdn_hrtf_boxdiags = @(oldfield) (logical(op.(oldfield) - 1));
% -----------------------------------------------------------

old_fields = fieldnames(new_name);
still_invalid = false(length(invalid_fields), 1);

% ---- Do the actual conversion: --------
for n = 1:length(invalid_fields)
    if ismember(invalid_fields{n}, old_fields)
        new_fldname = new_name.(invalid_fields{n});
        
        % Take value, if it is specified as a new option:
        if 1% ~isfield(op, new_fldname)
            op.(new_fldname) = newvalue.(invalid_fields{n})(invalid_fields{n});
        end
        warning(['Option "%s" is now handled by "%s".\nUse the new option.', ...
            ' See also the RAZR changelog.'], invalid_fields{n}, new_fldname);
    else
        still_invalid(n) = true;
    end
end
% ---------------------------------------

invalid_fields = invalid_fields(still_invalid);

end % eof
