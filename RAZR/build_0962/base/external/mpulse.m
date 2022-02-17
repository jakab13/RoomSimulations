% mpulse.m - generates spectrally pink-shaped pulse.
%
% Usage: out = mpulse(len, flow, fhigh, ntype, cweight, lnnoise, level, circular, fs)
%
% WARNING TEST VERSION FOR SPECIFIC PURPOSE. DO NOT WITHOUT CONSULTING SE
%
% len      = output length in samples
% flow     = lower cutoff frequency in Hz
% fhigh    = upper cutoff frequency in Hz
% ntype    = 0 : white
%            1 : pink
%            2 : red (brown)
%            3 : Bark-based uniform exciting
%            4 : ERB-based uniform exciting
%            5 : Threshold-equalizing noise, TEN(SPL)
% cweight  = if type is a scalar value, it specifies compensation weights
%            0 : none
%            1 : ISO threshold
%            2 : HDA200 retspl, flow > 125, fhigh < 16000
%            if type is a two-column vector, it specifies an
%            arbitrary frequency weighting
%            column 1 : frequency points
%            column 2 : linear amplitude weighting at frequency points.
%                       Linear interpolation is performed if less
%                       values then frequency bins are specified
% minphase = if 1 use minimum phase
% level    = in dB re full scale
% circular = [1/0] circular or fast using radix-2 FFT
% fs       = sampling rate in Hz
%
% out      = output vector 
%
% (c) Stephan Ewert 1999-2017, Universitaet Oldenburg. All rights reserved.
% Merged all all my 1999 pnoise, cnoise, bpnoise functions

function out = mpulse(len,lcut,hcut,ntype,cweight,minphase,level,circular,fs)


% find nice fft length if not circular
if ( circular == 1 )
    fftpts = len;
else
    fftpts = findnextpow2(len);
end

% freq spacing
spacing = fs/fftpts;

% nyquist freq bin, depending on even or odd number of samples
nyBin = floor(fftpts/2) + 1;

% vector of freq per bin up to Nyquist freq
freqVec = [0:nyBin-1]' * spacing;

lsmooth = [];
hsmooth = [];

if length(lcut) == 1
    lbin = max( round(lcut/spacing) + 1, 2);    % never add a DC, start at bin 2
else
    if lcut(1) >= lcut(2) 
        error('lcut(1) must be smaller than lcut(2)');
    end
    
    lbin = max( round(lcut(1)/spacing) + 1, 2);    % never add a DC, start at bin 2
    
    lbin2 = max( round(lcut(2)/spacing) + 1, 2);
    lsmooth = (0.5*(1 - cos(pi*[0:lbin2-lbin]'/(lbin2-lbin+1)))).^0.5;
end

if length(hcut) == 1
    hbin = min( round(hcut/spacing) + 1, nyBin );
else
    if hcut(1) >= hcut(2) 
        error('hcut(1) must be smaller than hcut(2)');
    end
    
    hbin = min( round(hcut(2)/spacing) + 1, nyBin );
    
    hbin2 = min( round(hcut(1)/spacing) + 1, nyBin );
    hsmooth = (0.5*(1 + cos(pi*[0:hbin-hbin2]'/(hbin-hbin2+1)))).^0.5;
end

% derive freq weighting depending on type
color_weight = zeros(nyBin,1);

switch (ntype )
    case 0
        color_weight = ones(nyBin,1);
    case 1
        color_weight(2:end) = 1./sqrt(freqVec(2:end));
    case 2
        color_weight(2:end) = 1./freqVec(2:end);
    case 3
        % Bark scale, calculate critical bandwidth for each frequency bin
        % according to
        % Zwicker, E., & Terhardt, E. (1980). Analytical expressions for critical-band rate and critical bandwidth as a function of frequency.
        % The Journal of the Acoustical Society of America, 68, 1523.
        cb = 25 + 75 *(1 + 1.4*(freqVec/1000).^2).^(0.69);
        color_weight = 1./sqrt(cb);
    case 4
        % ERB scale
        % according to
        % Glasberg and Moore, JASA 1990
        cb = 24.7 + freqVec/9.265;
        color_weight = 1./sqrt(cb);
    case 5
        % TEN(SPL) based on description in Moore, Glasberg, and Stone (2004),
        % "New Version of the TEN Test With Calibrations in dB HL" Ear & Hearing 25; 478–487

        % ERB scale according to
        % Glasberg and Moore, JASA 1990
        cb = 24.7 + freqVec(lbin:hbin)/9.265;
        color_weight(lbin:hbin) = 1./sqrt(cb);

        % correction by parameter "K" based on Figure 9 of
        % Moore, Glasberg, and Bear, 1997, AES; tabelized by SE
        paramK(:,1) = [50 60 70 100 150 200 300 400 500 700 1000 50000]'; % frequencies
        paramK(:,2) = [13.5 10.6 8.5 4.8 1.5 0.2 -1.2 -1.9 -2.5 -2.8 -3 -3]'; % values

        % derive freq weights based on K
        %paramK_weight = findlinearint(paramK, freqVec(lbin:hbin) );
        paramK_weight = interp1q(paramK(:,1),paramK(:,2),freqVec(lbin:hbin));

        % merge both with K relativ to 1000 Hz (-3 dB)
        color_weight(lbin:hbin) = color_weight(lbin:hbin)./10.^((paramK_weight + 3)/20);
    otherwise
        error('invalid noise type selected');
end

if max(size(cweight)) == 1
    % single number, thus it is cweight selection
    switch ( cweight )
        case 0
            % nothing
        case 1
            % Threshold in dB SPL for tones (ISO226)
            % tabelized by jens-e. appell / 3.95
            %
            % values from 20 Hz to 12500 Hz are taken from ISO 226 (1985-05-01)
            % values at 14000 Hz and 15000 Hz are taken from ISO-Threshold-table
            % in Klaus Bethges thesis.
            % values at 0 and 20000 Hz are not taken from ISO Threshold contour !!
            ISOTable(:,1) = 1000*[0.0    0.02 0.025 0.0315 0.04 0.05 0.063 0.08 0.1  0.125 0.16 0.2  0.25 0.315 0.4 0.5 0.63 0.8 1.0 1.25 1.6 2.0 2.5  3.15 4.0  5.0  6.3 8.0  10.  12.5    14.0 15.0    20.0]';
            ISOTable(:,2) = [ 80    74.3 65.0  56.3   48.4 41.7 35.5  29.8 25.1 20.7  16.8 13.8 11.2 8.9   7.2 6.0 5.0  4.4 4.2 3.8  2.6 1.0 -1.2 -3.6 -3.9 -1.1 6.6 15.3 16.4 11.6    16.0 24.1    70.0]';

            % interpolate freq weights
            %c_weight = findlinearint(ISOTable, freqVec(lbin:hbin) );
            c_weight = interp1q(ISOTable(:,1),ISOTable(:,2), freqVec(lbin:hbin) );

            % merge with color_weigth relativ to 1000 Hz (-4.2 dB)
            color_weight(lbin:hbin) = color_weight(lbin:hbin).*10.^((c_weight-4.2)/20);

        case 2
            % Sennheiser HDA 200 RetSPL
            % tabelized by SE taken from my abstreshtoaudiogram 11/2004

            % combined from ISO/DIS 389-8 and ISO/TR 389-5:1998(E)
            % 12000 is interpolated! by me, values at 0 and 50000 Hz are not
            % valid!
            retspltable = [ ...
                0       30.5; ...
                125     30.5; ...
                250     18.0; ...
                500     11.0; ...
                1000    5.5; ...
                2000    4.5; ...
                3000    2.5; ...
                4000    9.5; ...
                5000    14.0; ...
                6000    17.0; ...
                8000    17.5; ...
                9000    18.5; ...
                10000    22; ...
                11200    23; ...
                12000    26; ...
                12500    28; ...
                14000    36; ...
                16000    56; ...
                50000    56; ...
                ];

            % interpolate freq weights
            %c_weight = findlinearint(retspltable, freqVec(lbin:hbin) );
            c_weight = interp1q(retspltable(:,1),retspltable(:,2), freqVec(lbin:hbin) );

            % merge with color_weigth relativ to 1000 Hz (-5.5 dB)
            color_weight(lbin:hbin) = color_weight(lbin:hbin).*10.^((c_weight-5.5)/20);
        otherwise
            error('invalid cweight selected');
    end
elseif ( size(cweight,1) > 1 & size(cweight,2) == 2 )
    % arbitrary spectral coloration based on frequency/value pairs
    % based on my old cnoise

    % derive freq weights
    %c_weight = findlinearint(cweight, freqVec(lbin:hbin) );
    c_weight = interp1q(cweight(:,1),cweight(:,2), freqVec(lbin:hbin) );
    color_weight(lbin:hbin) = color_weight(lbin:hbin).*c_weight;
else
    error('invalid noise type or colorVec selected');
end

a = zeros(fftpts,1);

a(lbin:hbin) = ones(hbin-lbin+1,1);
fspec = a;

% SE 2017 fix at nyBin
if floor(fftpts/2) < fftpts/2
    % odd
else
    % even
    fspec(nyBin) =  fspec(nyBin)/2;
end


%normalize first than equalize
%fspec = fspec/sqrt(2*sum((abs(fspec)).^2)) * 10^(level/20);

%size(fspec)
%size(color_weight)


% assert 0 weight for DC (should never be applied due to lbin limits, but
% anyhow
color_weight(1) = 0;

fspec(lbin:hbin) = fspec(lbin:hbin) .* color_weight(lbin:hbin);

% smooth edges
if ~isempty(lsmooth)
    fspec(lbin:lbin2) = fspec(lbin:lbin2) .* lsmooth;
end

if ~isempty(hsmooth)
    fspec(hbin2:hbin) = fspec(hbin2:hbin) .* hsmooth;
end





if ( minphase > 0 )
    
    if floor(fftpts/2) < fftpts/2
        % odd
        fspec(floor(fftpts/2)+2:end)=conj(flipud(fspec(2:floor(fftpts/2)+1)));
    else
        % even
        fspec(fftpts/2+2:end)=conj(flipud(fspec(2:fftpts/2)));
        fspec(nyBin) = fspec(nyBin)*2;
    end
    
    %fspec2 = mps(fspec);
    
     % make minimum phase
     % fspec should be smooth at this point, if it reaches zero or has
     % strong changes we get time aliasing, limit to 60 dB dynamic range 
     floorl = max(abs(fspec))/10^(60/20);
     
     fspec = (fspec.^2+floorl^2).^0.5; % use root-of-sum-squares for smooth blending
     
     fspec = abs(fspec).*exp(-1i*imag(hilbert(log(abs(fspec)))));
    %fspec = abs(fspec).*exp(-1i*imag(hilbert(log(abs(fspec)+eps))));
    
    
 
    
%     if floor(fftpts/2) < fftpts/2
%         % odd
%         fspec(floor(fftpts/2)+2:end)=conj(flipud(fspec(2:floor(fftpts/2)+1)));
%     else
%         % even
%         fspec(fftpts/2+2:end)=conj(flipud(fspec(2:fftpts/2)));
%     end
   
    
   
end

out = fftpts*ifft(fspec);
out = 2*real(out(1:len));


% we need to get rid of time aliasing in the response. Take only the first
% quarter of the output, ramp the second quarter and zero the second half
if ( minphase > 0 )
    nyBin = nyBin-1;
    out(nyBin/2+1:nyBin) = out(nyBin/2+1:nyBin).*(0.5*(1 + cos(pi*[0:nyBin/2-1]'/(nyBin/2))));
    out(nyBin+1:end) = out(nyBin+1:end)*0;
end

%normalize level
out = out ./ rms(out) .* 10^(level/20);


% NOT USED, USING MATLAB INTERP1Q
% usage: value = findlinearint(vector, search);
function value = findlinearint(vector, search);

value = 0*search;

for i = 1:length(search)

    match = find( vector(:,1) == search(i) );

    if match
        value(i) = vector(match,2);
    else
        match_l = max(find( vector(:,1) < search(i) ));
        match_g = min(find( vector(:,1) > search(i) ));

        if isempty(match_l)
            value = [];
            error('interpolation failed at lower frequency bound');
            return;
        elseif isempty(match_g)
            value = [];
            error('interpolation failed at upper frequency bound');
            return;
        end

        value(i) = vector(match_l,2) + (vector(match_g,2) - vector(match_l,2))/(vector(match_g,1) - vector(match_l,1))*(search(i) - vector(match_l,1));
    end

end

%SE Matlab nextpow2 is not save 18.03.2014
function np2len = findnextpow2(len);

np2 = nextpow2(len);
np2len = 2^np2;
if np2len < len
    np2len = 2^(1 + np2);
end


% eof
