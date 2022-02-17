function addparam_func = get_addparam_func
% GET_ADDPARAM_FUNC - Get function handle for either ADDPARAMETER or
% ADDPARAMVALUE, depending on your current MATLAB version
%
% Usage:
%   addparam_func = GET_ADDPARAM_FUNC
%
% Output:
%   addparam_func       Either @addParameter or @addParamValue
%
% Example:
%   If you use the INPUTPARSER to parse Name-Value inputs of a function, you can
%   use GET_ADDPARAM_FUNC as follows:
%
%   % -----------------------
%   function myfunc(varargin)
%   p = inputParser;
%   addparam = get_addparam_func;
%   addparam(p, 'freq', 1e3);       % add parameter and set default value
%   parse(p, varargin{:});
%   f = p.Results.freq;             % get value of input 'freq'
%   % -----------------------
%
% See also: INPUTPARSER

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


mc = metaclass(inputParser);
method_names = cellfun(@(x){x.Name}, mc.Methods);

if any(strcmp(method_names, 'addParameter'))
    addparam_func = @addParameter;
else
    addparam_func = @addParamValue;
end
