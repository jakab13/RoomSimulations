% Example script for source directivity option. % Two source orientations 
% in two rooms are simulated and automatically played back. 



%% Please see EXAMPLE_DEFAULT first for a demonstration of the basic concepts of RAZR.
%%
clear;
clear global;
cfg = select_razr_cfg;

% Set options
% for receiver:
op.spat_mode = 'shm';
% for source:
op.dir_db = 'hs_filter'; 
op.ism_enableSourceDir = 1;      % Enable source directivity filtering?

% Load example room L:
room = get_room_L_Dir;

% Setup of source and receiver
room.srcpos = [2.00, 2.06, 1.50]; 
room.recpos = [1.00, 2.06, 1.50];
room.recdir	= [0,0];
room.srcdir = [180,0];

% generate ir front room L
ir_front_L = razr(room, op);

% rotate source
room.srcdir = [0,0];

% % generate ir back room L
ir_back_L = razr(room, op);


% Load example room A:
room = get_room_A_Dir;

% Setup of source and receiver
room.srcpos = [6.0, 7.50 , 1.2];
room.recpos = [6.0, 5.0 , 1.2];
room.recdir	= [90, 0];
room.srcdir = [-90,0];

% generate ir front room L
ir_front_A = razr(room, op);

% rotate source
room.srcdir = [90,0];

% % generate ir back room L
ir_back_A = razr(room, op);




%% auralize with test signal

% Room L
out_front_L = apply_rir(ir_front_L);
out_back_L  = apply_rir(ir_back_L);

% Room A
out_front_A = apply_rir(ir_front_A);
out_back_A  = apply_rir(ir_back_A);



%% Play audio

%Front Room L
disp('Playing audio: Front orientation Room L...')
sound(out_front_L{1},ir_front_L.fs);

pause(ceil(size(out_front_L{1},1)/ir_front_L.fs))

%Back Room L
disp('Playing audio: Back orientation Room L...')
sound(out_back_L{1},ir_back_L.fs);

pause(ceil(size(out_back_L{1},1)/ir_back_L.fs))

%Front Room A
disp('Playing audio: Front orientation Room A...')
sound(out_front_A{1},ir_front_A.fs);

pause(ceil(size(out_front_A{1},1)/ir_front_A.fs))

%Back Room A
disp('Playing audio: Back orientation Room A...')
sound(out_back_A{1},ir_back_A.fs);
