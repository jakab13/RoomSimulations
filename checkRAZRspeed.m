% script checkRAZRspeed
%
% Script to check how fast RAZR calculates BRIRs
%
% (C) 2021 by JMA, Johannes M. Arend
%             TU Berlin, Audio Communication Group

%% Room parameters

%For now, based on razr lab
clear all;

%Sample rate of FABIAN HRTF dataset
fs = 44100;

%Basic room settings
name = 'Lab';
dimensions = [6.97, 4.12, 3];
freq = [250 500 1e3 2e3 4e3];
temp = 20;
%abscoeff = [.02 .04 .06 .08 .1]; % Concrete rough;
abscoeff = [0.05, 0.1, 0.13, 0.16, 0.22];
srcPos = [2.40, 3.02, 1.80];
recPos1 = [2.40, 0.30, 1.80];
recDir	= [90, 0];
nRecPos = 3;

ismOrder = 3;

%% RAZR

%Room
room_razr.name = name;
room_razr.boxsize = dimensions;
room_razr.freq = freq;
room_razr.materials = abscoeff;
room_razr.srcpos = srcPos; %Static source

%Options
op.verbosity = 1;
op.spat_mode = 'hrtf';
op.hrtf_database = 'fabian.sofa'; %Path set in razr_cfg_user.m
op.fs = fs;
op.array_TCelsius = temp;
op.ism_order = ismOrder;
%Only ISM
op.ism_only = 1;
%Source directivity based on spherical head model
%op.ism_enableSourceDir = 1; 
%op.dir_db = 'hs_filter'; 
%room_razr.srcdir = [0,0];

stepWidth = 0.05;
steps = [0:stepWidth:2]; %2 meter on y axis in 5 cm steps
for s = 1:length(steps)
     
recPos1(2) = recPos1(2)+stepWidth;
%Update room position and direction in loop
room_razr.recpos(1,:) = recPos1;
room_razr.recdir(3,:) = recDir;

%Simulate BRIR
%tic
brir_razr = razr(room_razr,op);
%toc

end
