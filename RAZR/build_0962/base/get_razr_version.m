function razr_version = get_razr_version()
% GET_RAZR_VERISON - Returns the current version number as a string


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


file = [get_razr_path,filesep,'version.mat'];
razr_version = load(file,'Razr_Version'); %todo: modify with make_release
razr_version = razr_version.Razr_Version;
end
