function [out, in_respl] = apply_rir(ir, varargin)
% APPLY_RIR - Convolve RIR with one or several dry test signals.
%
% Usage:
%	[out, in_respl] = APPLY_RIR(ir)
%	[out, in_respl] = APPLY_RIR(ir, Name, Value)
%
% Input:
%	ir          RIR structure (see RAZR)
%
% Optional Name-Value-pair arguments:
%   src         Dry test signals or cell array of multiple test signals, specified as ...
%               - Soundfile names (including absolute or relative path if the file is not in the
%                 Matlab search path). If your Matlab version contains AUDIOREAD, all file formats
%                 supported by that function are possible, otherwise only wav files are supported.
%               - Key strings for sound samples specified in RAZR_CFG_DEFAULT and RAZR_CFG.
%                 By default, the following keys are supported: 'ppulse', 'olsa1', 'olsa2', 'olsa3'. 
%                 'ppulse' refers to a pink pulse signal.
%                 'olsa' refers to sound examples taken from the Oldenburg Sentence Test for the AFC
%                 Software Package by Stephan Ewert, Daniel Berg, Hoertech gGmbH. See also
%                 ANALYSIS_TOOLS/SAMPLES/OLSA_README.TXT.
%                 Own samples can be used, whose key strings and location on your harddrive have to
%                 be specified in RAZR_CFG.
%               - Sound signal as a matrix. The same sampling rate as of ir will be assumed.
%               Specifications of multiple sounds can be mixed, i.e. the cell array can look like,
%               e.g., {'olsa', rand(44100, 2), 'path/to/soundfile.wav', pink_pulse, ...}.
%               Default: 'olsa1'.
%   samples     [start, end] samples of audiofile(s) (specified in 'src') to read. If the file is
%               resampled (in order to match samplerate of ir), the samples apply to the original
%               samples of the file. Default: [1, Inf] (the whole file is read).
%               For multiple source signals, i.e. 'src' is a cell array, the same samples will
%               be used for all files. source_signals specified as matrices will not be limited.
%   fieldname   Fieldname in ir structure that is used as RIR signal. Default: 'sig'.
%   win         [attack, decay], lengths in samples of Hann flanks to be applied on 'src'.
%               Default: [0, 0].
%   eq          Headphone type key string. Currently supported: 'hd650' or 'none'. Own equalizations
%               can be used. For details, see RAZR_CFG_DEFAULT.
%               Default: Specified in RAZR_CFG_DEFAULT or, potentially overwritten, in RAZR_CFG
%   normalize   If true, normalize output signals to 0.99 (default: false)
%   cconvlen    For circular convolution (using cconv), cconvlen specifies the convolution length in
%               samples. If not specified or set to zero or empty, no circular convolution will be
%               performed. If set to -1, cconvlen will be set to each of the lengths of the dry
%               source signals, respectively.
%   thresh_rms_params = [thresh_db, maxval], where trhesh_db is the threshold for
%               thresh_rms calculation (see THRESH_RMS) and maxval is the maximum value the signal
%               is normalized to (has only an effect, if do_normalize == false). Default behaviour:
%               no application of THRESH_RMS.
%
% Output:
%	out         Cell array of convolved signals (sampling rate same as for ir)
%	in_respl	Dry source signals (resampled to same sampling rate as ir, if original differs)
%
% See also: RAZR_CFG_DEFAULT, THRESH_RMS, SOUNDIR, PINK_PULSE


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



%% input

if isempty(fieldnames(ir))
    return;
end

cfg = get_razr_cfg;
check_thresh = @(x) (length(x) == 2 || isempty(x));

p = inputParser;
addparam = get_addparam_func;
addparam(p, 'src', 'olsa5');
addparam(p, 'samples', [1, Inf]);
addparam(p, 'win', [0, 0]);
addparam(p, 'eq', cfg.default_headphone);
addparam(p, 'normalize', false);
addparam(p, 'cconvlen', []);
addparam(p, 'thresh_rms_params', [], check_thresh);
addparam(p, 'fieldname', 'sig');

parse(p, varargin{:});

if ~isfield(ir, p.Results.fieldname)
    error('Specified field does not exist: %s', p.Results.fieldname);
else
    ir_sig = ir.(p.Results.fieldname);
end

%% apply headphone eq

if ~strcmp(p.Results.eq, 'none')
    eq = load_headphone_eq(p.Results.eq, ir.fs);
    for n = 1:size(ir_sig, 2)
        ir_sig(:, n) = filter(eq, 1, ir_sig(:, n));
    end
end

%%
% convert 'src' to cell:
if iscell(p.Results.src)
    srcsig = p.Results.src;
else
    srcsig = {p.Results.src};
end

numSrc = length(srcsig);
in_respl = cell(numSrc, 1);
out      = cell(numSrc, 1);

if exist('audioread', 'file') || exist('audioread', 'builtin')
    audioread_fcn = @audioread;
else
    audioread_fcn = @wavread;
end

for n = 1:numSrc
    if isnumeric(srcsig{n})
        in = srcsig{n};
        fs = ir.fs;
    else
        fldname = sprintf('sample__%s', srcsig{n});
        if isfield(cfg, fldname)
            audiofilename = cfg.(fldname);
        elseif exist(srcsig{n}, 'file')
            audiofilename = srcsig{n};
        else
            error('Sound sample ID unknown or file not found: %s', srcsig{n});
        end
        [in, fs] = audioread_fcn(audiofilename, p.Results.samples);
    end
    
    % stereo2mono:
    if size(in, 2) == 2
        in = mean(in, 2);
    end
    
   % resample:
   if (size(ir)>1)
   
       if ir(n).fs ~= fs
            in_respl{n} = resample(in, ir(n).fs, fs);% (n) fix for multiple sources
       else
            in_respl{n} = in;
       end
       
   else
        
       if ir.fs ~= fs
            in_respl{n} = resample(in, ir.fs, fs);
       else
            in_respl{n} = in;
       end
   end

%     
    % flank:
    if ~isempty(p.Results.win)
        if length(p.Results.win) == 2
            in_respl{n} = flank_it(in_respl{n}, ...
                p.Results.win(1), p.Results.win(2));
        else
            error('win must be specified as [attacklen, decaylen]');
        end
    end
end

if ~isempty(p.Results.thresh_rms_params)
    in_respl = normalize_snd_samples(...
        in_respl, p.Results.thresh_rms_params(1), p.Results.thresh_rms_params(2));
end

do_circ_conv = ~isempty(p.Results.cconvlen) && p.Results.cconvlen ~= 0;
[len_ir, numCh] = size(ir_sig);

for n = 1:numSrc
    if ~do_circ_conv
        % pad zeros (space for reverb):
        in_respl{n} = [in_respl{n}; zeros(size(ir_sig, 1), 1)];
        len_src = size(in_respl{n}, 1);
        out{n} = zeros(max(len_src, len_ir), numCh);
        
        for ch = 1:numCh
            out{n}(:, ch) = dlfconvComp(in_respl{n}, ir_sig(:, ch));
        end
    else
        if p.Results.cconvlen == -1
            cconvlen = size(in_respl{n}, 1);
        else
            cconvlen = p.Results.cconvlen;
        end
        
        out{n} = zeros(cconvlen, numCh);
        
        for ch = 1:numCh
            out{n}(:, ch) = cconv(in_respl{n}, ir_sig(:, ch), round(cconvlen));
        end
    end
    
    if p.Results.normalize
        out{n} = 0.99*out{n}./(max(max(abs(out{n}))));
    end
end
end


function signal = flank_it(signal, attacklen, decaylen)

if attacklen > length(signal)
    error('Attack flank too long for signal.');
end
if decaylen > length(signal)
    error('Decay flank too long for signal.');
end

attackwin = hannwin(2*attacklen);
decaywin  = hannwin(2*decaylen);
signal(1:attacklen) = signal(1:attacklen).*attackwin(1:attacklen);
signal((end - decaylen + 1):end) = ...
    signal((end - decaylen + 1):end).*decaywin((decaylen + 1):(2*decaylen));

end
