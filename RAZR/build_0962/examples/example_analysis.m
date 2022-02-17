% Example script for demonstrating the usage of some RIR analysis tools.
% More functions for analyzing RIRs can be found in the folder ANALYSIS_TOOLS.
%
% See also: example_*, RAZR


%% Please see EXAMPLE_DEFAULT and EXAMPLE_OPTIONS first for a demonstration of the basic concepts
%% of RAZR.
%%

clear;

% This time we take a room for which a desired reverberation time is specified
% (have a look into that function for details):
room = get_room_RT;

% Create a sketch of the room (see help scene for more plotting options):
scene(room);

% Set some options:
op.rirname = room.name;         % Give a name to the BRIR
op.return_rir_parts = true;     % return direct sound, early and late BRIR parts as fields of ir

% BRIR synthesis:
ir = razr(room, op);

%% Analysis

% Plot BRIR time signal on dB scale, use separate colors for direct, early and late part
% (see help plot_ir for more plotting options):
plot_ir(ir, 'delg');

% Plot BRIR frequency response with 3rd-octave smoothing (see help plot_irspec for more options):
plot_irspec(ir, 'smo', '3rd');

% Calculate and plot T30 from BRIR, and plot EDC (see help schroeder_rt for more options):
[T30, freq, edc, tvec, h] = schroeder_rt(...
    ir, 'freq', room.freq, 'plot', true, 'plot_rt', true);

% Compare desired T60 with T30 calculated from synth. RIR:
hold(h.ax_rt, 'on');
h_desired = plot(h.ax_rt, room.freq, room.t60, ...
    'rs-', 'markerfacecolor', 'w', 'linewidth', 1.3);
legend([h.plot_rt, h_desired], {'Synthesized', 'Desired'});

%% Convolve BRIR with dry test signals

% There are multiple ways to chose dry test signals used for convolution in apply_rir. Please see
% help apply_rir for details (and also for more options):
test_snds = {'olsa2', randn(22050, 1), '06044.wav'};
[out, in] = apply_rir(ir, 'src', test_snds, 'normalize', true);

return

%% Listen to the n-th result by the following command:
n = 1;
soundsc(out{n}, ir.fs);

%% Listen to BRIR by the following command:
soundir(ir);
