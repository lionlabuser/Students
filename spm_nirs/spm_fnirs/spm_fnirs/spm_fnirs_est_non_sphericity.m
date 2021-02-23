function [xVi] = spm_fnirs_est_non_sphericity(Y, X, xVi, KY, KX)
% Non-sphericity estimation using ReML 
% FORMAT [xVi] = spm_fnirs_est_non_sphericity(Y, X, xVi, KY, KX) 
% 
% Y         measurements (eg, hemoglobin changes) 
% X         design matrix 
% xVi      structure describing intrinsic non-sphericity 
% KY       filtered measurements 
% KX       filtered design matrix 
% 
% xVi      estimated non-sphericity  
%
% This code is based on spm_est_non_sphericity.m 
% Karl Friston & Guillaume Flandin
% $Id: spm_est_non_sphericity.m 6015 2014-05-23 15:46:19Z guillaume $
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

%
% $Id$

if nargin < 4,
    KY = Y;
    KX = X;
end

ns = size(Y, 1);
nch = size(Y, 2); 

xKXs = spm_sp('Set', KX); 
xKXs.X = full(xKXs.X); 
pKX = spm_sp('x-', xKXs); 

beta = pKX * KY; 
res = KY - KX * beta; 
ResSS = sum(res.^2); 

trRV = spm_SpUtil('trRV', xKXs); 

Cy = 0; % sample covariance
q = spdiags(sqrt(trRV./ResSS'), 0, nch,nch);

Y = Y * q;
Cy = Cy + Y*Y';
Cy = Cy./nch;

[V, h] = spm_reml(Cy, X, xVi.Vi);

V      = V*ns/trace(V);

xVi.h = h; 
xVi.V = sparse(V); 
%xVi.Cy = Cy; 

