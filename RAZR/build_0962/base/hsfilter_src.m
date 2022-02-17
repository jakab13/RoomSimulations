function [output, out_magn] = hsfilter_src(theta, elev, theta_ear, ...
    Fs, dir_src_radius, dir_theta0, input, normalise)
% based on: C. P. Brown and R.O. Duda, "A Structural Model for Binaural 
% Sound Synthesis," IEEE Trans. Speech Audio Processing, vol. 6, no. 5, 
% September 1998.
%
% Usage:
%   [output, out_magn] = hsfilterMod(theta, elev, theta_ear, warpMethod, ...
%       Fs, input, options)
%
% Input:    
%   theta           azimuth angle source
%   elev            elevation angle source (not used right now)
%   theta_ear       azimuth angle of ear re nose in degrees 
%                       (in horizontal plane)
%	Fs              Sampling frequency in Hz
% 	input           Input signal
%  	normalise       normalise to 0 dB at 0° (true false)
%
% Output:   
% 	output          Output of filtered input signal according
%                       to head shadowing
%	out_magn           Magnitude of the output signal

%------------------------------------------------------------------------------
% RAZR engine for Mathwork's MATLAB
%
% Version 0.96.2
%
% Author(s): Henning Steffens
%
% Copyright (c) 2014-2021, Torben Wendt, Steven van de Par, Stephan Ewert,
% University of Oldenburg, Germany.
%
% This work is licensed under the
% Creative Commons Attribution-NonCommercial-NoDerivs 4.0 International
% License (CC BY-NC-ND 4.0).
% To view a copy of this license, visit
% http://creativecommons.org/licenses/by-nc-nd/4.0/ or send a letter to
% Creative Commons, 444 Castro Street, Suite 900, Mountain View, California,
% 94041, USA.
%------------------------------------------------------------------------------

[len, chN] = size(input); % number of input channels

% make sure that theta and elev are row vectors:
theta = theta(:)';

%theta0 = 180;   % (150 in Zölzer's DAFX) 
%a = 0.08;       % radius of head 0.08; %0.04
%c = 340;        % speed of sound
alfa_min = 0.05; % (0.05 in Zölzer's DAFX) %0.1
c = speedOfSound();
theta0 = dir_theta0;
a = dir_src_radius;
w0 = c/a;       % frequency related to radius of head


% head shadow filter
alfa = 1 + alfa_min/2 + (1- alfa_min/2) * cos(theta/theta0*pi);  

B = [(alfa+w0/Fs)/(1+w0/Fs), (-alfa+w0/Fs)/(1+w0/Fs)];  % numerator of TF
A = [1, -(1-w0/Fs)/(1+w0/Fs)];  % denominator of TF

% normalise zero azim
if(normalise)
alfa_norm = 1 + alfa_min/2 + (1- alfa_min/2) * cos(0/theta0*pi);  
B_norm = [(alfa_norm +w0/Fs)/(1+w0/Fs), (-alfa_norm +w0/Fs)/(1+w0/Fs)];  % numerator of TF
A_norm = [1, -(1-w0/Fs)/(1+w0/Fs)];  % denominator of TF
end

% SE remove delays, we assume a point source rotating around its acoustic center 
% % calculate delay in samples
% if (abs(theta) < 90)
%     gdelay = - Fs/w0*(cos(theta*pi/180) - 1);  % the -1 compensates for 
%                                                % negative non-causal group 
%                                                % delay re center of the sphere
% else
%     gdelay = Fs/w0*((abs(theta) - 90)*pi/180 + 1); % here the +1 does the same
% end
% 
% % separate the delay in integer samples and fractional delay
% gdelay_samples = floor(gdelay);
% gdelay = gdelay - gdelay_samples;   % ranges from 0 to 1;
% 
% % derive first order allpass filter coeff ( for fractional delay)
% % ( H(z) = a + z^-1/(1+a*z^-1) )
% a = (1 - gdelay)./(1 + gdelay);



%% 
sig_range = repmat([1, len], length(theta), 1);

% pseudo input for magnitude and output filter below
out_magn_temp = input;
output_temp = input;

% filters
for idx = 1:length(theta)
    % check for mono- or multichannel input
    if chN == 1
        % magnitude filter for mono input
        out_magn_temp(sig_range(idx, 1):sig_range(idx, 2), idx) = filter(B(idx,:), A, input(sig_range(idx, 1):sig_range(idx, 2)));
    elseif length(theta) ~= chN
        error('multichannel signal: number of angles must match number of channels')
    else
        % magnitude filter for multichannel input
        out_magn_temp(sig_range(idx, 1):sig_range(idx, 2),idx) = filter(B(idx,:), A, input(sig_range(idx, 1):sig_range(idx, 2), idx));
    end
    
    % normalise zero azim filter
    if(normalise)
        out_magn_temp(sig_range(idx, 1):sig_range(idx, 2), idx) = filter(A_norm, B_norm(idx,:), out_magn_temp(sig_range(idx, 1):sig_range(idx, 2)));
    end
    
% SE remove delays, we assume a point source rotating around its acoustic center    
%     % integer sample delay
%     output_temp(sig_range(idx, 1):sig_range(idx, 2),idx) = ...
%         filter([zeros(gdelay_samples(idx,:), 1); 1], 1, ...
%         out_magn_temp(sig_range(idx, 1):sig_range(idx, 2), idx));
%     
%     % excess fractional delay via first-order allpass
%     output_temp(sig_range(idx, 1):sig_range(idx, 2),idx) = ...
%         filter([a(idx,:), 1],[1, a(idx,:)], ...
%         output_temp(sig_range(idx, 1):sig_range(idx, 2), idx));
end

out_magn = out_magn_temp;
output = out_magn_temp;

