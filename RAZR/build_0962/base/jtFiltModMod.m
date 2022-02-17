% jtFilt - Get filter coefficients for parametric EQs, to approximate specified
% frequency response. 
% Mod: Extrapolate additional command gains for shelves (computationally
% costly)
% Implementation is based on
% Oliver, R. and Jot, J-M.: Efficient Multi-Band Digital Audio Graphic
% Equalizer with Accurate Frequency Response Control, AES Convention Paper
% 9406 (2015)
%
% Usage:
%   [B, A, b, a] = jtFilt(freq, dGain, fs, linGain)
%
% Input:
%   freq            Frequency base in Hz, must have the same column length as
%                   dGain
%   dGain           Desired frequency response, values in dB or linear (specified by input parameter 'linGain')
%   fs              Sampling rate
%   linGain         Set true, if dGain values are linear attenuation factors, false if they are specified in dB
%
% Output:
%   B, A            Filter coefficients of the shelving filters as Matrix (each filter in one row)
% CK development version, 2020-03-31
