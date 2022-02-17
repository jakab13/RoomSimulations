function [B, A] = getReflFiltCoeff(freq, absolRefl, b_bp, a_bp, op)
% GETREFLFILTCOEFF - Get reflection filter coefficients
%
% Usage:
%   [bRefl, aRefl] = GETREFLFILTCOEFF(freq, absolRefl, b_bp, a_bp, fs, method, makePlot)
%
% Input:
%   absolRefl           Reflection coefficients as matrix; rows: walls
%   freq                Frequency base
%   b_bp, a_bp          Coefficients for global bandpass filter, applied as last step
%   method              Method for filter synthsis
%                       'cs': composedShelving,
%                       'sh': shEQ,
%                       'cq': compodesPEQ (recommended),
%                       'yw': getFiltCoeff_yulewalk
%   makePlot            If true, make test plot
%
% Output:
%   bRefl, aRefl        Filter coefficients as cell array
%                       (for method 'sh' filter length may vary)

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

if op.debugMode
    global bugOut;
    global whereAmI;
end

numCoeffs = size(absolRefl, 1);
bRefl = cell(numCoeffs, 1);
aRefl = bRefl;

switch op.filtCreatMeth
    case 'jt'
        for w = 1:numCoeffs
            [B(:,w),A(:,w)] = jtFiltModMod(freq,absolRefl(w,:),op.fs,true);
        end
        if op.debugMode && strcmp(whereAmI,'FDNsetup')
            bugOut.ReflFilterDebug.FDN.method = method; 
        elseif op.debugMode && strcmp(whereAmI,'ISMsetup')
            bugOut.ReflFilterDebug.ISM.method = op.filtCreatMeth; 
        end
    case 'cq'
        
        for w = 1:numCoeffs
            [B(:,w), A(:,w), bRefl{w}, aRefl{w}] = composedPEQ(freq, absolRefl(w,:), op.fs, true, false);
            % filter optimisation parameter set to false for ease of data
            % management, CK 20-03-30
            if op.debugMode && strcmp(whereAmI,'FDNsetup')
                
                bugOut.ReflFilterDebug.FDN.method = op.filtCreatMeth;                
                bugOut.ReflFilterDebug.FDN.aOG{w} = aRefl{w};
                bugOut.ReflFilterDebug.FDN.bOG{w} = bRefl{w};
                bugOut.ReflFilterDebug.FDN.h(w,:) = absolRefl(w,:);
%                 bugOut.ReflFilterDebug.FDN.cqA{:,w} = A;
%                 bugOut.ReflFilterDebug.FDN.cqB{:,w} = B;
%                                
%                 bugOut.rflFiltFDN.cq.aOG{w} = aRefl{w};
%                 bugOut.rflFiltFDN.cq.bOG{w} = bRefl{w};
%                 bugOut.rflFiltFDN.cq.h(w,:) = absolRefl(w,:);
%                 bugOut.rflFiltFDN.cq.cqA{:,w} = cqA;
%                 bugOut.rflFiltFDN.cq.cqB{:,w} = cqB;

            elseif op.debugMode && strcmp(whereAmI,'ISMsetup')
                
                bugOut.ReflFilterDebug.ISM.method = op.filtCreatMeth;   
                bugOut.ReflFilterDebug.ISM.aOG{w} = aRefl{w};
                bugOut.ReflFilterDebug.ISM.bOG{w} = bRefl{w};
                bugOut.ReflFilterDebug.ISM.h(w,:) = absolRefl(w,:);
                bugOut.ReflFilterDebug.ISM.cqA{:,w} = A;
                bugOut.ReflFilterDebug.ISM.cqB{:,w} = B;     
                
%                 bugOut.rflFiltISM.cq.aOG{w} = aRefl{w};
%                 bugOut.rflFiltISM.cq.bOG{w} = bRefl{w};
%                 bugOut.rflFiltISM.cq.h(w,:) = absolRefl(w,:);
%                 bugOut.rflFiltISM.cq.cqA{:,w} = cqA;
%                 bugOut.rflFiltISM.cq.cqB{:,w} = cqB;     
            end
            bRefl{w} = conv(bRefl{w}, b_bp);
            aRefl{w} = conv(aRefl{w}, a_bp);
        end
        tmpABP = cell(1,numCoeffs);
        tmpBBP = tmpABP;
        tmpABP(:) = {a_bp};
        tmpBBP(:) = {b_bp};
            B(end+1,:) = tmpBBP;
            A(end+1,:) = tmpABP;
        if op.debugMode && strcmp(whereAmI,'FDNsetup')
            bugOut.ReflFilterDebug.FDN.aBP = aRefl;
            bugOut.ReflFilterDebug.FDN.bBP = bRefl;
%             
%             bugOut.rflFiltFDN.cq.aBP = aRefl;
%             bugOut.rflFiltFDN.cq.bBP = bRefl;
        elseif op.debugMode && strcmp(whereAmI,'ISMsetup')
           
             bugOut.ReflFilterDebug.ISM.aBP = aRefl;
             bugOut.ReflFilterDebug.ISM.bBP = bRefl;

%             bugOut.rflFiltISM.cq.aBP = aRefl;
%             bugOut.rflFiltISM.cq.bBP = bRefl;
        end
    case 'cs'
        for w = 1:numCoeffs
            [bb, aa] = composedShelving(absolRefl(w,:), freq, op.fs, 'peak', 1, 1);
            bRefl{w} = conv(bb, b_bp);
            aRefl{w} = conv(aa, a_bp);
            A(1,w) = {aa};
            B(1,w) = {bb};
            A(2,w) = {a_bp};
            B(2,w) = {b_bp};
            if op.debugMode && strcmp(whereAmI,'FDNsetup')
                bugOut.ReflFilterDebug.FDN.method = op.filtCreatMeth;    
                bugOut.ReflFilterDebug.FDN.aa{:,w} = aa;
                bugOut.ReflFilterDebug.FDN.bb{:,w} = bb;
                bugOut.ReflFilterDebug.FDN.aRefl{w} = aRefl{w};
                bugOut.ReflFilterDebug.FDN.bRefl{w} = bRefl{w};
            elseif op.debugMode && strcmp(whereAmI,'ISMsetup') 
                bugOut.ReflFilterDebug.ISM.method = op.filtCreatMeth;  
                bugOut.ReflFilterDebug.ISM.aa{:,w} = aa;
                bugOut.ReflFilterDebug.ISM.bb{:,w} = bb;
                bugOut.ReflFilterDebug.ISM.aRefl{w} = aRefl{w};
                bugOut.ReflFilterDebug.ISM.bRefl{w} = bRefl{w};
               
            end
                
        end
        
    case 'yw'
        [bb, aa] = getFiltCoeff_yulewalk(freq, absolRefl, op.fs, 1*length(freq), 0);
        for w = 1:numCoeffs
            bRefl{w} = conv(bb(w,:), b_bp);
            aRefl{w} = conv(aa(w,:), a_bp);
            A(1,w) = {aa(w,:)'};
            B(1,w) = {bb(w,:)'};
            A(2,w) = {a_bp};
            B(2,w) = {b_bp};
           if op.debugMode && strcmp(whereAmI,'FDNsetup')
                bugOut.ReflFilterDebug.FDN.method = op.filtCreatMeth;    
                bugOut.ReflFilterDebug.FDN.aa{:,w} = aa;
                bugOut.ReflFilterDebug.FDN.bb{:,w} = bb;
                bugOut.ReflFilterDebug.FDN.aRefl{w} = aRefl{w};
                bugOut.ReflFilterDebug.FDN.bRefl{w} = bRefl{w};
            elseif op.debugMode && strcmp(whereAmI,'ISMsetup') 
                bugOut.ReflFilterDebug.ISM.method = op.filtCreatMeth;  
                bugOut.ReflFilterDebug.ISM.aa{:,w} = aa;
                bugOut.ReflFilterDebug.ISM.bb{:,w} = bb;
                bugOut.ReflFilterDebug.ISM.aRefl{w} = aRefl{w};
                bugOut.ReflFilterDebug.ISM.bRefl{w} = bRefl{w};
            end
            
            
        end
        
    case 'sh'
        for w = 1:numCoeffs
            [bb, aa] = shEQ(freq, absolRefl(w,:), op.fs, 1, 0);
            bRefl{w} = conv(bb, b_bp);
            aRefl{w} = conv(aa, a_bp);
            A(1,w) = {aa'};
            B(1,w) = {bb'};
            A(2,w) = {a_bp};
            B(2,w) = {b_bp};
           if op.debugMode && strcmp(whereAmI,'FDNsetup')
                bugOut.ReflFilterDebug.FDN.method = op.filtCreatMeth;    
                bugOut.ReflFilterDebug.FDN.aa{:,w} = aa;
                bugOut.ReflFilterDebug.FDN.bb{:,w} = bb;
                bugOut.ReflFilterDebug.FDN.aRefl{w} = aRefl{w};
                bugOut.ReflFilterDebug.FDN.bRefl{w} = bRefl{w};
            elseif op.debugMode && strcmp(whereAmI,'ISMsetup') 
                bugOut.ReflFilterDebug.ISM.method = op.filtCreatMeth;  
                bugOut.ReflFilterDebug.ISM.aa{:,w} = aa;
                bugOut.ReflFilterDebug.ISM.bb{:,w} = bb;
                bugOut.ReflFilterDebug.ISM.aRefl{w} = aRefl{w};
                bugOut.ReflFilterDebug.ISM.bRefl{w} = bRefl{w};
               
            end
            
        end
        
    otherwise
        error('Method not available: %s.', method);
end
if op.debugMode && strcmp(whereAmI,'FDNsetup')
    bugOut.ReflFilterDebug.FDN.A = A;
    bugOut.ReflFilterDebug.FDN.B = B;
elseif op.debugMode && strcmp(whereAmI,'ISMsetup')
    bugOut.ReflFilterDebug.ISM.A = A;
    bugOut.ReflFilterDebug.ISM.B = B;
end

if strcmp(op.fltMode,'conv')
    bb = cell(1,numCoeffs);
    bb(:) = {1};
    aa = bb;
    for ch = 1:numCoeffs         % convolve cascade into single filter
        for k = 1:size(A(:,ch),1)
            bb{ch} = conv(bb{ch}, B{k,ch});
            aa{ch} = conv(aa{ch}, A{k,ch});
        end
    end
    A = aa; % return summed values instead of cascade
    B = bb;
elseif strcmp(op.fltMode,'smart')
   for ch = 1:numCoeffs         % convolve cascade into single filter
    [bbb(:,ch),aaa(:,ch)] = mergeFilters(B(:,ch),A(:,ch),9);
   end
   B = bbb;
   A = aaa;
elseif strcmp(op.fltMode,'casc') || strcmp(op.fltMode,'sos') % don't do anything
else
    error(['unknown filter mode: ' op.fltMode])
end

% test plot:
if op.plot_filters_fdn_bin %% ck 20-04-08 hier nochmal eleganter fixen
    figure
    for w = 1:numCoeffs
        p1 = semilogx(freq,20*log10(absolRefl(w,:)),'s','Color',[1/w 0.5 w/numCoeffs]);
        hold on
        p2(w) = plotFrqRsp(bRefl{w}, aRefl{w}, op.fs, 'mode', 'sgl', 'color', get(p1,'Color'));
        hold on
    end
    hold off
    title(sprintf('Reflection factors; Filter synth. method: %s', op.filtCreatMeth), 'interpreter', 'none')
    ylim([-6 0])
    legend show;
    legend('location', 'best')
end
