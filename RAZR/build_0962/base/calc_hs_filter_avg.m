% Calculation of average source directivity filter (uncomment to save)

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

clear
close all
clc

% Sphere Integral
radius = 1; %0.08
phi = 360;
theta = 0:2.5:90.0; % 2.5° resolution (widh of rings)
   
phi = dg2rd(phi);
theta = dg2rd(theta);    

fun = @(theta,phi) radius^2 * sin(theta);

surface_integral_half_sphere = zeros(1,size(theta,2)-1);

n = 1;

for i = 2:size(surface_integral_half_sphere,2)+1
     
lwr_phi = 0;
upr_phi = phi;
   
lwr_theta = theta(i-1);
upr_theta = theta(i);
   
surface_integral_half_sphere(1,n) = integral2(fun,lwr_theta,upr_theta,lwr_phi,upr_phi);
n=n+1;
end

 
surface_integral_sphere = sum(surface_integral_half_sphere)*2; % -> 4*pi


%% Power Spectral Density of hs_filter

N = 4096;
fftlen= N;

input = zeros(N,1);
input(1) = 1;
fs = 44100;
freq = (0:N-1).*(fs/N); 
freqs = freq(freq < (fs/2)+1);  
normalise = 1;

dir_src_radius = 0.08; 
dir_theta0 = 180;

degree = 0:2.5:177.5; % rotationally symmetric 


PSD = zeros(size(freq,2),size(degree,2));


for k = 1:size(degree,2) 

    h1 = hsfilter_src(degree(:,k), 0, 0, fs, dir_src_radius, dir_theta0, input,normalise);
    
    H1 = fft(h1,fftlen);
    
    %psd/freqbin
    psdx = (abs(H1).^2);
    %psd/freq bin over degree
    PSD(:,k) = psdx;


end




%% PSD-SPHERE

%ring surfaces whole sphere
surface_integral_full_sphere = [surface_integral_half_sphere(1:end) surface_integral_half_sphere(end:-1:1)]; 

%Weighting of PSD Matrix with ring surfaces  
psd_sphere = surface_integral_full_sphere .* PSD;

%sums of weighted psds per freqbins
psd_sphere_sum = sum(psd_sphere,2);

% divided by surface of whole sphere for average
avg_pf =  psd_sphere_sum / (sum(surface_integral_full_sphere));

% sqrt of ring weighted average psd/freqbin
avg_pf_sqr = sqrt(avg_pf);

%min. phase
avg_pf_sqr_phase = avg_pf_sqr .* exp(-1i .* imag(hilbert(log(avg_pf_sqr))));

%impulse response (real() because of small rounding errors)
ir_avg_hs_filter = real(ifft(avg_pf_sqr_phase));


%% PLOT

figure;
% hold on
semilogx(freqs,20*log10(abs(avg_pf_sqr_phase(1:(N/2)+1))))
% 
 xlim([63 22000]);
 ylim([-35 5]);
 grid on
 title('Avg. Filter (Resolution 2.5 degree)')
 xlabel('Frequency (Hz)')
 ylabel('20*log10*abs(H) / dB')

set(gcf,'color','w');   

%% SAVE
%save('ir_avg_hs_filter.mat','ir_avg_hs_filter');

 
        
        
        
        
        
       
