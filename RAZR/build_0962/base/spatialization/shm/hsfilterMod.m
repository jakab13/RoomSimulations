% based on: C. P. Brown and R.O. Duda, "A Structural Model for Binaural 
% Sound Synthesis," IEEE Trans. Speech Audio Processing, vol. 6, no. 5, 
% September 1998.
%
% Usage:
%   [output, out_magn] = hsfilterMod(theta, elev, theta_ear, warpMethod, ...
%       Fs, input, options)
%
% Input:    
%   theta           azimuth angle re nose in degrees (positive numbers
%                       to right)
%   elev            elevation angle re horizontal plane in degrees
%                       (positive numbers up), range [90 ... -90]
%   theta_ear       azimuth angle of ear re nose in degrees 
%                       (in horizontal plane)
%	warpMethod   	prewarp theta to better account for smaller 
%                       changes in attenuation at small theta
%                       warpMethod = 0 --> do not warp
%                       warpMethod = 1 --> exponent transform warp
%                       warpMethod = 2 --> alternative warping
%	Fs              Sampling frequency in Hz
% 	input           Input signal
%  	limit_range_start_spls
%                   Vector of length size(in, 2) containing time samples
%                       (for each column of 'in') at which the channel-wise
%                       filtering starts. From this sample on, filtering is
%                       applied on the doubled filter length. Using this option
%                       saves computation time. Default: filtering on the whole
%                       signal.
%
% Output:   
% 	output          Output of filtered input signal according
%                       to head shadowing
%	out_magn   		Magnitude of the output signal

%------------------------------------------------------------------------------
% RAZR engine for Mathwork's MATLAB
%
% Version 0.96.2
%
% Author(s): Oliver Buttler
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
