function example_abscoeff(spec)
% EXAMPLE_ABSCOEFF - Example to demonstrate how absorption coefficients can be
% specified.
%
% Usage:
%   EXAMPLE_ABSCOEFF(spec)
%
% Input:
%   spec        Specification switch, see code for details.
%
% See also: example_*, RAZR, MATERIAL2ABSCOEFF


%% Please see EXAMPLE_DEFAULT first for a demonstration of the basic concepts of RAZR.
%%

% Create a room:
room = get_room_L;

% Now, we specify the rooms' wall absorption coefficients in the chosen way:

switch spec
    case 1
        %% Cell array, key strings
        % May be row or column cell array
        % For details and available material-key strings, see MATERIAL2ABSCOEFF,
        % GET_ABSCOEFF_HALL
        room.materials = {...
            'hall:brick', ...
            'hall:draperies', ...
            'hall:draperies', ...
            'hall:brick', ...
            'hall:draperies', ...
            'hall:draperies', ...
            };
        
    case 2
        %% Row vector
        % Same for all walls, frequency dependent (room.freq)
        room.materials = [0.7, 0.6, 0.7, 0.7, 0.5];
        
    case 3
        %% Column vector
        % Same for all frequencies, different for walls
        % Order of walls: [-z;  -y;  -x;   +x;  +y;  +z]
        room.materials  = [0.7; 0.6; 0.7; 0.4; 0.5; 0.1];
        
    case 4
        %% Matrix
        % Different for frequencies and walls
        room.materials = [...
            0.0600    0.1500    0.4000    0.6000    0.6000; ... -z
            0.1000    0.0500    0.0400    0.0700    0.1000; ... -y
            0.3000    0.1000    0.1000    0.1000    0.1000; ... -x
            0.2000    0.2000    0.1000    0.0700    0.0400; ... +x
            0.7000    0.6000    0.7000    0.7000    0.5000; ... +y
            0.0100    0.0200    0.0200    0.0200    0.0300  ... +z
            ];
        
    case 5
        %% Cell array, mixed
        % May be row or column cell array
        room.materials = {...
            'hall:brick', ...
            'hall:draperies', ...
            [0.7, 0.6, 0.7, 0.7, 0.5], ...  % frequency dependent
            'hall:draperies', ...
            0.5, ...                        % frequency independent
            [0.7; 0.6; 0.7; 0.7; 0.5], ...  % will be treated as frequency-dependent, although column vector
            };
        
    case 6
        %% Cell array, numeric values
        % Special case of "mixed", may be row or column cell array
        room.materials = {...
            [0.0600    0.1500    0.4000    0.6000    0.6000]; ... -z
            [0.1000    0.0500    0.0400    0.0700    0.1000]; ... -y
            [0.3000    0.1000    0.1000    0.1000    0.1000]; ... -x
            [0.2000    0.2000    0.1000    0.0700    0.0400]; ... +x
            [0.7000    0.6000    0.7000    0.7000    0.5000]; ... +y
            [0.0100    0.0200    0.0200    0.0200    0.0300]  ... +z
            };
        
        %% Desired reverberation time instead of absorption coefficients
        % Please see GET_ROOM_RT for details
        
end

% Add field "abscoeff" to room:
room = add_abscoeff(room);

% For all cell-array specifications, materials will be plotted as labels in scene:
if iscell(room.materials)
    scene(room);
end
