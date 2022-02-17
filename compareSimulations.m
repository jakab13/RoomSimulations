% script compareSimulations
%
% Script to compare various shoebox room acoustic simulation engines
%
% (C) 2021 by JMA, Johannes M. Arend
%             TU Berlin, Audio Communication Group

%% Room parameters

%For now, based on razr lab
clear all;

%Sample rate of FABIAN HRTF dataset
fs = 44100;

%Test samples
applyHPCF = true;
ns = 2;
s1 = audioread('olsa2_44.wav');
s1 = s1*0.5; %Reduce by 6 dB just to have enough headroom (more important with the drum signal)
s2 = audioread('drums_44.wav');
s2 = s2*0.5;

if applyHPCF
    hpcf = audioread('HpFilter_fabian_HD600.wav');
    s1 = supdeq_convFFT(s1,hpcf);
    s2 = supdeq_convFFT(s2,hpcf);
end
sig = {s1,s2};

%Basic room settings
name = 'Lab';
dimensions = [6.97, 4.12, 3];
freq = [250 500 1e3 2e3 4e3];
temp = 20;
%abscoeff = [.02 .04 .06 .08 .1]; % Concrete rough;
abscoeff = [0.05, 0.1, 0.13, 0.16, 0.22];
srcPos = [4.40, 3.02, 1.80];
recPos1 = [4.40, 0.30, 1.80];
recPos2 = [4.40, 2.30, 1.80];
recPos3 = [5.40, 2.30, 1.80];
recDir	= [90, 0];
nRecPos = 3;

ismOrder = 3;

%% RAZR

%Room
room_razr.name = name;
room_razr.boxsize = dimensions;
room_razr.freq = freq;
room_razr.materials = abscoeff;
room_razr.srcpos = srcPos;
room_razr.recpos(1,:) = recPos1;
room_razr.recpos(2,:) = recPos2;
room_razr.recpos(3,:) = recPos3;
room_razr.recdir(1,:) = recDir;
room_razr.recdir(2,:) = recDir;
room_razr.recdir(3,:) = recDir;

%Options
op.verbosity = 1;
op.spat_mode = 'hrtf';
op.hrtf_database = 'fabian.sofa'; %Path set in razr_cfg_user.m
op.fs = fs;
op.array_TCelsius = temp;
op.ism_order = ismOrder;
%Only ISM
%op.ism_only = 1;
%Source directivity based on spherical head model
%op.ism_enableSourceDir = 1; 
%op.dir_db = 'hs_filter'; 
%room_razr.srcdir = [0,0];

%Plot room
close all;
scene(room_razr)

%Simulate BRIRs
brir_razr = razr(room_razr,op);

%Convolve with anechoic material
nCond = 0;
for s = 1:ns
    for k = 1:nRecPos
       
       nCond = nCond+1;
       
       sig_s = sig{s}; 
       brir = brir_razr(k).sig;
       sig_razr{nCond} = supdeq_convFFT(sig_s,brir);
       
    end
end


%% SofaMyRoom (SMR)

%Room settings
room_smr.dimension = dimensions;
room_smr.humidity = 0.42;
room_smr.temperature = temp; %Also default in RAZR

room_smr.surface.frequency = freq;
room_smr.surface.absorption = repmat(abscoeff,6,1);
room_smr.surface.diffusion = repmat([0.5 0.5 0.5 0.5 0.5],6,1);%0 Diffusion, 0.5 default -- Diffusion takes of course energy from ISM early reflections -- ISM not correct anymore. Only correct with diffusion = 0, but then, ray tracing is wrong and BRIRs look corrupt

%Simulation Options
options_smr.fs                  = fs;                   % sampling frequency in Hz
options_smr.responseduration    = 2;%? 1.25 defult               % duration of impulse response
options_smr.bandsperoctave      = 1;                    % simulation frequency accuracy (1, 2, 3, or 4 bands/octave)
options_smr.referencefrequency  = 125;   %125 default               % reference frequency for frequency octaves
options_smr.airabsorption       = true;                 % apply air absorption?
options_smr.distanceattenuation = true;                 % apply distance attenuation?
options_smr.subsampleaccuracy   = false;                % apply subsample accuracy?
options_smr.highpasscutoff      = 0;                    % 3dB frequency of high-pass filter (0=none)
options_smr.verbose             = true;                 % print status messages?

options_smr.simulatespecular    = true;                 % simulate specular reflections?
options_smr.reflectionorder     = [ ismOrder ismOrder ismOrder ];%Default 10         % maximum specular reflection order (x,y,z)

options_smr.simulatediffuse     = true;                 % simulate diffuse reflections?
options_smr.numberofrays        = 2000; %?                % number of rays in simulation (20*K^2)
options_smr.diffusetimestep     = 0.010;  %?              % time resolution in diffuse energy histogram (seconds)
options_smr.rayenergyfloordB    = -80;      %?            % ray energy threshold (dB, with respect to initial energy)
options_smr.uncorrelatednoise   = true;                 % use uncorrelated poisson arrivals for binaural impulse responses?

options_smr.outputname			= 'output';           	% name of the output file
options_smr.mex_saveaswav       = false;                % enable or disable saving the results of sofamyroom on disk
                                                    % when using MATLAB
%Source
source_smr(1).location           = srcPos;       % location of source (x,y,z; meters)
source_smr(1).orientation        = [ -90 0 0 ];         % orientation of source (yaw,pitch,roll; degrees)
source_smr(1).description        = 'omnidirectional';       % source type

%Receiver
receiver_smr(1).location         = recPos1;        % location of receiver (x,y,z; meters)
receiver_smr(1).orientation      = [ recDir, 0 ];           % orientation of receiver (yaw,pitch,roll; degrees)
receiver_smr(1).description      = 'SOFA C:\Users\JMA\sciebo\TUB\RoomSimulations\hrirs\FABIAN_HRIR_measured_HATO_0.sofa';  % receiver type
receiver_smr(2).location         = recPos2;        % location of receiver (x,y,z; meters)
receiver_smr(2).orientation      = [ recDir, 0 ];           % orientation of receiver (yaw,pitch,roll; degrees)
receiver_smr(2).description      = 'SOFA C:\Users\JMA\sciebo\TUB\RoomSimulations\hrirs\FABIAN_HRIR_measured_HATO_0.sofa';  % receiver type
receiver_smr(3).location         = recPos3;        % location of receiver (x,y,z; meters)
receiver_smr(3).orientation      = [ recDir, 0 ];           % orientation of receiver (yaw,pitch,roll; degrees)
receiver_smr(3).description      = 'SOFA C:\Users\JMA\sciebo\TUB\RoomSimulations\hrirs\FABIAN_HRIR_measured_HATO_0.sofa';  % receiver type

%Write everything in setup struct and simulate
setup_smr.room = room_smr;
setup_smr.options = options_smr;
setup_smr.source = source_smr;
setup_smr.receiver = receiver_smr;

%Plot room
close all;
plotroom(setup_smr)

%Simulate BRIRs
brir_smr = sofamyroom(setup_smr);

%Convolve with anechoic material
nCond = 0;
for s = 1:ns
    for k = 1:nRecPos
       
       nCond = nCond+1;
       
       sig_s = sig{s}; 
       brir = brir_smr{k};
       sig_smr{nCond} = supdeq_convFFT(sig_s,brir);
       
    end
end


%% RAZR - Listen

idRecPos = 1;
idSig = 2;

if idSig == 1
    idListen = idRecPos;
elseif idSig == 2
    idListen = idRecPos+nRecPos;
end
sound(sig_razr{idListen},fs);

%% SofaMyRoom - Listen

idRecPos = 1;
idSig = 2;

if idSig == 1
    idListen = idRecPos;
elseif idSig == 2
    idListen = idRecPos+nRecPos;
end
sound(sig_smr{idListen}*1.75,fs);



