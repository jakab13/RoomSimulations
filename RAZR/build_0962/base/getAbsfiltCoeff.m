function [B, A, h] = getAbsfiltCoeff(room, fdn_setup, op)
% GETABSFILTCOEFF - Get absorption filter coefficients for FDN
%
% Usage:
%   [b, a, h] = GETABSFILTCOEFF(room, fdn_setup, op)
%
% Input:
%   room        room structure (see RAZR)
%   fdn_setup   Output of GET_FDN_SETUP
%   op          Options structure (see RAZR)
%
% Output:
%   b, a        Filter coefficients as matrix; rows: channels
%   h           Frequency responses

%------------------------------------------------------------------------------
% RAZR engine for Mathwork's MATLAB
%
% Version 0.96.2
%
% Author(s): Torben Wendt
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


%% test call
if op.debugMode
    global bugOut;
    global whereAmI;
end

if nargin == 0
    room.boxsize  = [5, 4, 3];
    room.materials = {...
        'hall:concrete_block_painted', 'hall:draperies', ...
        'hall:carpet_on_conc', 'hall:windowglass', ...
        'hall:concrete_block_painted', 'hall:tile'};
    
    room.freq = octf(250, 4e3);
    room.abscoeff = material2abscoeff(room.materials, room.freq);
    fdn_setup.delays = [900; 1200; 1400];
    
    fdn_setup.b_bp = 1;
    fdn_setup.a_bp = 1;
    
    op.fs = 44100;
    op.filtCreatMeth = 'cq';
    op.plot_absFilters = true;
    op.rt_estim = 'eyring';
    
    [b, a, h] = getAbsfiltCoeff(room, fdn_setup, op);
    return;
end

%%

numFrq  = size(room.freq,2);
numCh   = length(fdn_setup.delays);
b_bplen = length(fdn_setup.b_bp);
a_bplen = length(fdn_setup.a_bp);

if isfield(room, 't60')  % room.materials had their chance in complement_room()
    rt = room.t60;
else
    rt = estimate_rt(room, op.rt_estim);
end

% frequency responses according to Jot and Chaigne (1991) (linear scale):
h = 10.^(-3/op.fs./rt(:)*fdn_setup.delays(:)')';

switch op.filtCreatMeth
    case 'jt'
        for ch = 1:numCh
            [B(:,ch),A(:,ch)] = jtFiltModMod(room.freq,h(ch,:),op.fs,true);
        end
            if op.debugMode
                
                bugOut.AbsFilterDebug.FDN.method = op.filtCreatMeth;
                bugOut.AbsFilterDebug.FDN.cqA{:,ch} = A;
                bugOut.AbsFilterDebug.FDN.cqB{:,ch} = B;
                bugOut.AbsFilterDebug.FDN.h = h;    
                 
            end
    case 'cq'
        % currently still without prealloc:
        %b = zeros(numCh, 2*numFrq-1);
        %a = b;
        bb = cell(numCh,1);
        aa = cell(numCh,1);
        
        for ch = 1:numCh
            [B(:,ch), A(:,ch), bb{ch}, aa{ch}] = composedPEQ(room.freq, h(ch,:), op.fs, true, false);
            % filter optimisation parameter set to false for ease of data
            % management, CK 20-03-30
            if op.debugMode
                
                bugOut.AbsFilterDebug.FDN.method = op.filtCreatMeth;
                bugOut.AbsFilterDebug.FDN.cqA{:,ch} = A;
                bugOut.AbsFilterDebug.FDN.cqB{:,ch} = B;
                bugOut.AbsFilterDebug.FDN.h(ch,:) = h(ch,:);    
                
%                 bugOut.absFilt.cq.cqA{:,ch} = cqA;
%                 bugOut.absFilt.cq.cqB{:,ch} = cqB;
%                 bugOut.absFilt.cq.h(ch,:) = h(ch,:);   
                
            end
        end
        blen = cellfun(@length, bb);
        maxblen = max(blen);
        b = zeros(numCh, maxblen);
        a = zeros(numCh, maxblen);
        
        for ch = 1:numCh
            b(ch,:) = [bb{ch}, zeros(1, maxblen-blen(ch))];
            a(ch,:) = [aa{ch}, zeros(1, maxblen-blen(ch))];
        end
        
        if any(fdn_setup.b_bp ~= 1) || any(fdn_setup.a_bp ~= 1)
            fprintf('Warning in %s: global bandpass in FDN not implemented for method cq\n', mfilename);
        end
        if op.debugMode
                                  
            bugOut.AbsFilterDebug.FDN.cqA{:,ch} = A;
            bugOut.AbsFilterDebug.FDN.cqB{:,ch} = B;
            bugOut.AbsFilterDebug.FDN.h(ch,:) = h(ch,:);    
            
            bugOut.AbsFilterDebug.FDN.aOG = aa';
            bugOut.AbsFilterDebug.FDN.bOG = bb';
            bugOut.AbsFilterDebug.FDN.a = a;
            bugOut.AbsFilterDebug.FDN.b = b;
            
            
%             bugOut.absFilt.cq.cqA = cqA;
%             bugOut.absFilt.cq.cqB = cqB;
%             bugOut.absFilt.cq.aOG = aa';
%             bugOut.absFilt.cq.bOG = bb';
%             bugOut.absFilt.cq.a = a;
%             bugOut.absFilt.cq.b = b;
        end
    case 'cs'
        % prealloc:
        b = zeros(numCh, numFrq*2+1 + b_bplen-1);
        a = zeros(numCh, numFrq*2+1 + a_bplen-1);
        
        for ch = 1:numCh
            [bb, aa] = composedShelving(h(ch,:), room.freq, op.fs, 'peak', 1, 1);
            b(ch,:) = conv(bb, fdn_setup.b_bp);
            a(ch,:) = conv(aa, fdn_setup.a_bp);
            A(1,ch) = {aa};
            B(1,ch) = {bb};
            A(2,ch) = {fdn_setup.a_bp};
            B(2,ch) = {fdn_setup.b_bp};
            if op.debugMode
                
                bugOut.AbsFilterDebug.FDN.method = op.filtCreatMeth;
                bugOut.AbsFilterDebug.FDN.aOG{:,ch} = aa';
                bugOut.AbsFilterDebug.FDN.bOG{:,ch} = bb';
            end
            
        end
        
        if op.debugMode
                
                bugOut.AbsFilterDebug.FDN.method = op.filtCreatMeth;
                bugOut.AbsFilterDebug.FDN.a = a;
                bugOut.AbsFilterDebug.FDN.b = b;
        end
        
    case 'yw'
        % prealloc:
        M = 1;
        b = zeros(numCh, numFrq*M+1 + b_bplen-1);
        a = zeros(numCh, numFrq*M+1 + a_bplen-1);
        
        [bb, aa] = getFiltCoeff_yulewalk(room.freq, h, op.fs, M*length(room.freq), 0);
        
        if op.debugMode
                
                bugOut.AbsFilterDebug.FDN.method = op.filtCreatMeth;
                bugOut.AbsFilterDebug.FDN.aOG = aa;
                bugOut.AbsFilterDebug.FDN.bOG = bb;
        end
        
        for ch = 1:numCh
            b(ch,:) = conv(bb(ch,:), fdn_setup.b_bp);
            a(ch,:) = conv(aa(ch,:), fdn_setup.a_bp);
            A(1,ch) = {aa(ch,:)'};
            B(1,ch) = {bb(ch,:)'};
            A(2,ch) = {fdn_setup.a_bp};
            B(2,ch) = {fdn_setup.b_bp};
        if op.debugMode
        
            bugOut.AbsFilterDebug.FDN.a{ch,:} = a(ch,:);
            bugOut.AbsFilterDebug.FDN.b{ch,:} = b(ch,:);
            bugOut.AbsFilterDebug.FDN.A = A;
        end
            
        end
        
    case 'sh'
        % prealloc:
        b_c = cell(numCh,1);
        a_c = b_c;
        flens = zeros(1,numCh);			% filter lengths
        
        for ch = 1:numCh
            [bb, aa] = shEQ(room.freq, h(ch,:), op.fs, 1);
            b_c{ch} = conv(bb, fdn_setup.b_bp);
            a_c{ch} = conv(aa, fdn_setup.a_bp);
            flens(ch) = length(b_c{ch});
            A(1,ch) = {aa'};
            B(1,ch) = {bb'};
            A(2,ch) = {fdn_setup.a_bp};
            B(2,ch) = {fdn_setup.b_bp};
            if op.debugMode
                
                bugOut.AbsFilterDebug.FDN.method = op.filtCreatMeth;
                bugOut.AbsFilterDebug.FDN.aOG = aa;
                bugOut.AbsFilterDebug.FDN.bOG = bb;
                bugOut.AbsFilterDebug.FDN.a_c{ch} = a_c{ch};
                bugOut.AbsFilterDebug.FDN.b_c{ch} = b_c{ch};
                
                
            end
            
        end
        
        % write coefficients into matrix (pad rows with zeros)
        maxflen = max(flens);
        b = zeros(numCh, maxflen);
        a = b;
        
        for ch = 1:numCh
            b(ch,:) = [b_c{ch}, zeros(1,maxflen-flens(ch))];
            a(ch,:) = [a_c{ch}, zeros(1,maxflen-flens(ch))];
            
             if op.debugMode
  
               bugOut.AbsFilterDebug.FDN.a{ch,:} = a(ch,:);
               bugOut.AbsFilterDebug.FDN.b{ch,:} = b(ch,:);
                
                
            end
            
        end
        
    otherwise
        error('Method not available: %s.', op.filtCreatMeth);
end
if op.debugMode
    bugOut.AbsFilterDebug.FDN.A = A;
    bugOut.AbsFilterDebug.FDN.B = B;
end
% CK 2020-04-08
% Annotation about data structur of filter coefficients:
% As of now, all valid options for op.filtCreatMeth yield an assortment of
% filter coefficients in two cell arrays (one for B, one for A). 
% Along the first dimension, different stages of the filter cascade are
% described, FDN channels resp. walls for reflection filters have separate
% columns.
% For the filtCreatMeths 'cq', 'cs' and 'jt' these single cells already
% correspond to second order structures, for the other methods there might
% be one higher order filter and then something like the global bandpass
% (if existent).
% Depending on the filtering strategy, these cells are now either left alone
% for direct use in 'casc' or later rearrangement into a SOS matrix, merged 
% into one % filter per wall  (still in cell array) for 'conv' or an 
% arrangement of filters with a maximum order of 9 to avoid numerical artefacts
% while keeping decent performance for 'smart'.
if strcmp(op.fltMode,'conv')
    bb = cell(1,numCh);
    bb(:) = {1};
    aa = bb;
    for ch = 1:numCh         % convolve cascade into single filter
        for k = 1:size(A(:,ch),1)
            bb{ch} = conv(bb{ch}, B{k,ch});
            aa{ch} = conv(aa{ch}, A{k,ch});
        end
    end
    A = aa; % return summed values instead of cascade
    B = bb;
elseif strcmp(op.fltMode,'smart')
    clear aaa bbb
   for ch = 1:numCh         % convolve cascade into single filter
    [bbb(:,ch),aaa(:,ch)] = mergeFilters(B(:,ch),A(:,ch),9);
   end
   B = bbb;
   A = aaa;
elseif strcmp(op.fltMode,'casc') || strcmp(op.fltMode,'sos') % don't do anything
else
    error(['unknown filter mode: ' op.fltMode])
end
if op.plot_absFilters
    figure
    han = semilogx(room.freq, 20*log10(h), 's', 'linewidth', 2);
    hold on
    hf = plot_freqrsp(B, A, 'ax', gca, 'fs', op.fs, 'disp_mode', 'sgl');
    
    for n = 1:numCh
        hf.plot(n).Color = han(n).Color;
    end
    
    hold off
    title(sprintf('FDN absorption filters; Synth. method: %s', op.filtCreatMeth));
    %ylim([-10 0])
end
