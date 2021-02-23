function [hbo, hbr, hbt] = spm_fnirs_calc_hb(X, P)
% Calculate hemoglobin changes using the modified Beer-Lambert law
% FORMAT [hbo, hbr, hbt] = spm_fnirs_calc_hb(X, P)
%
% X         optical density changes 
% P         structure array of model parameters 
%
% hbo     oxy-hemoglobin changes 
% hbr      deoxy-hemoglobin changes 
% hbt      total hemoglobin changes 
%__________________________________________________________________________
%
% References: 
% Modified Beer-Lambert law:
% Delpy, D.T., Cope, M., Van der Zee, P., Arridge, S., Wray, S.,Wyatt, J., 1988. 
% Estimation of optical pathlength through tissue from direct time of flight measurement. 
% Phys. Med. Biol. 33 (12), 1433-1442.
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

%
% $Id$

if isstruct(X), X = X.od; end
        
ns = size(X, 1); nch = size(X, 2); 

a = diag(P.dpf) * P.acoef; 
b = diag(1./((a(2,2)*a(1,1) - a(1,2)*a(2,1)) .* P.d));

hbo = (a(2,2) * X(:,:,1) - a(1,2) * X(:,:,2)) * b; 
hbr = (a(1,1) * X(:,:,2) - a(2,1) * X(:,:,1)) * b;
hbt = hbo + hbr; 
