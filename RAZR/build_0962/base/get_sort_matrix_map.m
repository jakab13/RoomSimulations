function gainTable = get_sort_matrix_map(ismData, fdnSetup, signalString)
% GET_SORT_MATRIX - Calculates channel mapping of ISM to FDN
% channels.
%
% Usage:
%   S_delta = GET_SORT_MATRIX_MAP(angles, valid_angles)
%
% Input:
%   ismData         for positions image sources centered around receiver
%   fdnSetup        for positions of reverb sources centered around receiver 
%   signalString    to check current signal matrix
%
% Output:
%   gainTable         gain matrix for ism -> fdn

%------------------------------------------------------------------------------
% RAZR engine for Mathwork's MATLAB
%
% Version 0.96.2
%
% Author(s): Michael Schutte, Henning Steffens
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



if (strcmp(signalString,'curr_sigmat'))
    
    ism_relpos = ismData.relpos((ismData.order == max(ismData.order)),:);   
    % L2 normalization optional: output stimmt überein mit ism_relpos
    ism_relpos = ism_relpos ./ norm(ism_relpos,'fro');
    ism_relpos = (-1).*ism_relpos;
     
elseif (strcmp(signalString,'early_refl_sigmat_diffuse'))
    
    ism_relpos = ismData.relpos;
   
end

fdn_relpos = fdnSetup.relpos;
fdn_channels = zeros(size(fdn_relpos,1),1);
ism_channels = zeros(size(ism_relpos,1),3);
squared_norm = 0;
image_source_gains = zeros(size(fdn_channels,1),size(ism_channels,1));

% 
for k = 1:size(ism_channels,1)
    squared_norm = 0;
    for kk = 1:size(fdn_channels,1)
        dot_product = dot(fdn_relpos(kk,:),ism_relpos(k,:));
        if (dot_product > 0)
            squared_norm = squared_norm +(dot_product * dot_product);
            image_source_gains(kk,k) = dot_product;
        end
        
    end
    
    if (squared_norm ~= 0)
        inverse_norm = 1/sqrt(squared_norm);
        for jj = 1:size(fdn_channels,1)
            
            image_source_gains(jj,k) = image_source_gains(jj,k) * inverse_norm;
            
        end
        
    end
end

gainTable = image_source_gains';

end
