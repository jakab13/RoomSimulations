function room = add_abscoeff(room)
% ADD_ABSCOEFF - Add absorption coefficient matrix to a room, based on material
% specifications. Absorption data are retrieved from databases stored in
% get_abscoeff_*.m.
%
% Usage:
%   room = ADD_ABSCOEFF(room)
%
% Input:
%   room    Room structure (see RAZR) which must contain the fields "materials"
%           and "freq". room.freq is the frequency base (in Hz) for absorption
%           coefficients. room.materials can be:
%           - A cell array of 6 material strings in the format "dbase.material",
%             where "dbase" is the name of an absorbing-material database that
%             contains absorption coefficients stored under the name "material".
%             For details and an example, see GET_ROOM_L.M.
%           - A row vector containing absorption coefficients that match the
%             frequencies stored in room.freq. These values will be copied to be
%             equal for all walls.
%           - A column vector containing six absorption coefficients being the
%             same for all frequencies but different for the six walls.
%           - A matrix containing absorption coefficients for six walls and a
%             number of length(freq) frequencies. Walls in rows, frequencies in
%             columns. In that case, this matrix is simply copied to
%             room.abscoeff.
%           - A single number representing an absorption coefficient that will
%             be used for all walls and frequencies.
%           - A cell array, one entry per wall, with each entry being either a 
%             material string, a row vector, or a scalar number.
%
%           For examples, see EXAMPLE_ABSCOEFF.M.
%
% Output:
%   room    Room structure with new field "abscoeff" containing absorption
%           coefficients. Walls in rows, frequencies in columns.
%
% See also: MATERIAL2ABSCOEFF, get_abscoeff_*, EXAMPLE_ABSCOEFF, SMAT

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


%% Input parameters

if ~isfield(room, 'freq')
    error('room must contain field ""freq".');
end
if ~isfield(room, 'materials')
    error('room must contain field ""materials".');
end

num_walls = 6;
num_freq = length(room.freq);

%%
if ischar(room.materials)
    room.materials = {room.materials};
end

if iscell(room.materials)
    if length(room.materials) == 1
        room.materials = repmat(room.materials, 1, num_walls);
    elseif length(room.materials) ~= num_walls
        error('There must be either 1 or %d materials be specified in a room.', ...
            num_walls);
    end
    
    isstr = cellfun(@ischar, room.materials);
    lens  = cellfun(@length, room.materials);
    
    if any(~isstr & lens ~= num_freq & lens ~= 1)
        error('Number of absorption coefficients must be 1 or match number of frequencies.');
    end
    
    % ensure suitable dimensionalities:
    room.materials = room.materials(:);
    room.materials(~isstr) = ...
        cellfun(@(x) x(:)', room.materials(~isstr), 'UniformOutput', false);
    
    issingle = ~isstr & (lens == 1);
    room.materials(issingle) = ...
        cellfun(@(x) repmat(x, 1, num_freq), room.materials(issingle), 'UniformOutput', false);
    
    room.abscoeff = zeros(num_walls, num_freq);
    room.abscoeff(~isstr, :) = cell2mat(room.materials(~isstr));
    room.abscoeff( isstr, :) = material2abscoeff(room.materials(isstr), room.freq);
    
elseif isnumeric(room.materials)
    if numel(room.materials) == 1
        room.abscoeff = repmat(room.materials, num_walls, num_freq);
        
    elseif isrowvec(room.materials)
        % same for all walls, freq dependent:
        if length(room.materials) ~= num_freq
            error('Number of materials must match number of frequencies.');
        end
        room.abscoeff = repmat(room.materials, num_walls, 1);
        
    elseif iscolvec(room.materials)
        % same for all frequencies
        if length(room.materials) ~= num_walls;
            error('Number of materials must be %d.', num_walls);
        end
        room.abscoeff = repmat(room.materials, 1, num_freq);    
    
    else
        room.abscoeff = room.materials;
    end
end

if any(any(room.abscoeff >= 1))
    error('Absorption coefficients must be smaller than 1');
end
