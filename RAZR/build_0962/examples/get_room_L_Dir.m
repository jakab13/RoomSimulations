function room = get_room_L_Dir
% GET_ROOM_RT - Definintion of an example room: Laboratory.
% In contrast to GET_ROOM_L, here a desired reverberation time 
%
% Usage:
%   room = GET_ROOM_L
%
% Output:
%   room    Room structure (see RAZR)
%
% See also: GET_ROOM_A, GET_ROOM_L, SCENE


% For details on the fields see GET_ROOM_L or razr help.
room.name = 'Laboratory';
room.boxsize = [4.97, 4.12, 3];
% 
% Specify a frequency dependent desired reverberation time instead of wall materials.
% Note: If both materials and t60 are fields of a room, razr will use materials and ignore t60.
room.t60  = [0.42, 0.46, 0.42, 0.37, 0.38];
room.freq = [250   500   1e3   2e3   4e3];

room.srcpos = [2.00, 2.06, 1.50]; 
room.recpos = [1.00, 2.06, 1.50];
room.recdir	= [0,0];
room.srcdir    = [180,0];

