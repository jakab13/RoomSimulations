% Example script for defining and auralizing different source distances. 
% Two source distances in two rooms are simulated and automatically played 
% back. 


%% Please see EXAMPLE_DEFAULT first for a demonstration of the basic concepts of RAZR.
%%
clear;
clear global;
cfg = select_razr_cfg;

filename_list = [
    "bark"
    "bum" 
    "chirp" 
    "glass" 
    "lock" 
    "pinknoise" 
    "plug" 
    "sneeze" 
    "waterdrop" 
    "whisper" 
    "wow"
    ];

% distances = [0, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8];
distances = 0: 0.2: 20;

size_x = 30;
size_y = 30;
size_z = 10;

recpos_x = size_x/2;
recpos_y = 3.0;
recpos_z = 1.5;

filename_room = "_room-" + num2str(size_x) + "-" + num2str(size_y) + "-" + num2str(size_z);

for i = 1:length(filename_list)
    
    filename_core = filename_list(i);

    % Load example room L:
    room = get_room_A;
    room.boxsize   = [size_x, size_y, size_z];

    % Setup of source and receiver, distance 2.5 m
    room.recpos = [recpos_x, recpos_y , recpos_z];

    out_root_folder_name = filename_core + filename_room;
    out_root_folder_path = "/Users/jakabpilaszanovich/Documents/GitHub/distance_perception/experiment/samples/" + out_root_folder_name;
    out_simulated_folder_path = out_root_folder_path + "/simulated";

    if ~exist(out_simulated_folder_path, 'dir')
       mkdir(out_simulated_folder_path)
    end

    for distance = distances
        if distance == 0
            room.materials = [0.99, 0.99, 0.99, 0.99, 0.99];
            room.srcpos = [recpos_x, recpos_y + 0.05, recpos_z];
            filename_full = filename_core + filename_room + "_control.wav";
        else
            room.materials = [0.05, 0.1, 0.13, 0.16, 0.22];
            room.recdir	= [90, 0];
            room.srcpos = [recpos_x, recpos_y + distance, recpos_z];
            filename_full = filename_core + filename_room + "_dist-" + num2str(distance * 100) + ".wav";
        end
        ir = razr(room);
        out = apply_rir(ir, 'src', filename_core);

        out_filename_path = out_simulated_folder_path + "/" + filename_full;
        audiowrite(out_filename_path, out{1}, ir.fs);
        disp("writing: " + filename_full);
    end

    disp("done");
end