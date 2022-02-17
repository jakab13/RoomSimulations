% Example script for defining and auralizing different source distances. 
% Two source distances in two rooms are simulated and automatically played 
% back. 


%% Please see EXAMPLE_DEFAULT first for a demonstration of the basic concepts of RAZR.
%%
clear;
clear global;
cfg = select_razr_cfg;


% Load example room L:
room = get_room_L;

% Setup of source and receiver, distance: 1 m
room.srcpos = [2.00, 2.06, 1.50]; 
room.recpos = [1.00, 2.06, 1.50];
room.recdir	= [0, 0];

% generate ir 1 m room L
ir_one_L = razr(room);

% distance: 3 m
room.srcpos = [4.00, 2.06, 1.50]; 

% % generate ir 3 m room L
ir_three_L = razr(room);


% Load example room A:
room = get_room_A;

% Setup of source and receiver, distance 2.5 m
room.srcpos = [6.0, 7.50 , 1.2];
room.recpos = [6.0, 5.0 , 1.2];
room.recdir	= [90, 0];

% generate ir 2.5 m room A
ir_twofive_A = razr(room);

% distance: 7.5 m
room.srcpos = [6.0, 12.50 , 1.2];

% generate ir 7.5 m room A
ir_sevenfive_A = razr(room);




%% auralize with test signal

% Room L
out_one_L   = apply_rir(ir_one_L);
out_three_L = apply_rir(ir_three_L);

% Room A
out_twofive_A   = apply_rir(ir_twofive_A);
out_sevenfive_A = apply_rir(ir_sevenfive_A);



%% Play audio

disp('Playing audio distance 2.5 m Room A...')
% Distance 2.5 m Room A
sound(out_twofive_A{1},ir_twofive_A.fs);
disp('wait...')
pause(ceil(size(out_twofive_A{1},1)/ir_twofive_A.fs))
disp('Playing audio distance 7.5 m Room A...')
% Distance 7.5 m Room A
sound(out_sevenfive_A{1},ir_sevenfive_A.fs);
disp('wait...')
pause(ceil(size(out_sevenfive_A{1},1)/ir_sevenfive_A.fs))
disp('Playing audio distance 1 m Room L...')
% Distance 1 m Room L
sound(out_one_L{1},ir_one_L.fs);
disp('wait...')
pause(ceil(size(out_one_L{1},1)/ir_one_L.fs))
disp('Playing audio distance 3 m Room L...')
% Distance 3 m Room L
sound(out_three_L{1},ir_three_L.fs);
disp('finish')


