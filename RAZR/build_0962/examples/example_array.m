% Example script for creating a multichannel RIR suited for playback in a
% loudspeaker array.
% The example array provided here, was also used in, e.g., [1]
%
% References:
%   [1] Denk, F., Ernst, S. M. A., Ewert, S. D., Kollmeier, B. (2018).
%       Adapting Hearing Devices to the Individual Ear Acoustics: Database and
%       Target Response Correction Functions for Various Device Styles,
%       Trends in Hearing 22, DOI: 10.1177/2331216518779313.
%   [2] Pulkki, V. (1997). Virtual Sound Source Positioning Using Vector Base
%       Amplitude Panning. J. Aud. Eng. Soc., 45(6), 456-466.
%   [3] Laitinen, M., Vilkamo, J., Jussila, K., Politis, A., Pulkki, V. (2014). 
%       Gain normalization in amplitude panning as a function of frequency and 
%       room reverberance. 55th Int. Conference of the AES. Helsinki, Finland.
%
% See also: example_*, RAZR


%% Please see EXAMPLE_DEFAULT and EXAMPLE_OPTIONS first for a demonstration of
%% the basic concepts of RAZR.
%%

clear;

room = get_room_L;
op.return_rir_parts = true;     % return direct sound, early and late BRIR parts as fields of ir
op.return_op = 1;
% Load loudspeaker-array coordinates and assign them to the respective option
% (Positions of speakers relative to receiver, cartesian coordinates):
load('ls_array_vrroom_uol.mat');  % loads "array_pos"
op.array_pos = array_pos;

% Specify spatialization for loudspeaker array:
op.spat_mode = 'array';

% Specify to perform also binaural synthesis of the array signals to enable
% pre-listening via headphones:
op.enable_post_spat = 1;
op.post_spat_mode = 'shm';  % 'shm', 'hrtf', or diotic (for 'hrtf', the database
                            % specified in op.hrtf_database will be used)

% Specify the rendering method:
% 'nearest': "round" virtual directions to nearest speakers, respectively.
% 'vbap': Apply Vector Base Amplitude Panning (after Pulkki [2])
op.array_render = 'nearest';
% To use VBAP, enter 'vbap'. Code is included under /base/external/vbap.
% For up-to-date version, see:
% https://github.com/polarch/Vector-Base-Amplitude-Panning


op.vbap_color_compen = 1;  % Apply coloration compensation for VBAP? [3]

ir = razr(room, op);


% Generate binaural-synthesis playback signal:
if op.enable_post_spat
    out = apply_rir(ir, 'fieldname', 'sig_post_spat');
end

% Generate multichannel playback signal:
out_array = apply_rir(ir);

%% Plot array and power of loudspeaker signals
visWrap(ir,24,1)
return

%% Playback test signal (use ctrl + ENTER to run this section)
soundsc(out{1}, ir.fs);


