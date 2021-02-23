function [SPM] = spm_fnirs_spm(SPM) 
% Estimate GLM parameters for each channel from fNIRS measurements 
% Channel-wise estimates are then interpolated using spm_fnirs_spm_interp.m
% FORMAT [SPM] = spm_fnirs_spm(SPM)
%
% SPM       structure array of GLM parameters 
%
% Note: This script is based on spm_spm.m 
% Karl Friston & Guillaume Flandin
% $Id: spm_spm.m 6015 2014-05-23 15:46:19Z guillaume $
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

%
% $Id$

%--------------------------------------------------------------------------
% Get SPM file 

if ~nargin, % specify file name of SPM.mat
    [SPM sts] = spm_select(1, '^SPM*.*.\.mat$', 'Select a SPM file for fNIRS');
    if ~sts, SPM = []; return; end 
end

if ischar(SPM), load(SPM); end 

fprintf('--------------------------------------------------------- \n'); 
fprintf('Estimate parameters for SPM-%s analysis...\n', SPM.xY.type); 
fprintf('--------------------------------------------------------- \n'); 

%--------------------------------------------------------------------------
% Get design (X) and measurement (Y)
fprintf('Apply a temporal filter to measurements and design matrix... \n'); 

xX = SPM.xX; 
xVi = SPM.xVi; 
ns = SPM.nscan; 
fs = 1/SPM.xY.RT; 

load(SPM.xY.VY); 

Y = reshape(spm_vec(rmfield(Y, 'od')), [P.ns P.nch 3]); 
Y = spm_fnirs_preproc(Y, P); 

% identify channels of interest 
load(P.fname.pos); 
mask = zeros(1, R.ch.nch); 
indx = find(R.ch.mask == 1); 
mask(R.ch.label(indx)) = 1; 

mask = mask .* P.mask; 
chs = find(mask == 1); 

switch SPM.xY.type
    case 'HbO'
        Y = Y(:, chs, 1); 
    case 'HbR'
        Y = Y(:, chs, 2); 
    case 'HbT'
        Y = Y(:, chs, 3); 
end
clear od R; 

KY = spm_fnirs_filter(Y, P, fs); 
KX = spm_fnirs_filter(xX.X, P, fs); 

fprintf('Completed. \n\n'); 

W = speye(ns); 
% estimate intrinsic temporal correlation using AR(1) model 
if ~isfield(xVi, 'V') 
    fprintf('Estimate intrinsic temporal correlation... \n'); 
    [xVi] = spm_fnirs_est_non_sphericity(Y, xX.X, xVi, KY, KX);

    %-Get weight/whitening matrix:  W*W' = inv(V)
    %--------------------------------------------------------------------------
    W          = spm_sqrtm(spm_inv(xVi.V));
    W          = W.*(abs(W) > 1e-6);
    W = sparse(W); 
    
    % update KX, KY
    KY = spm_fnirs_filter(W*Y, P, fs); 
    KX = spm_fnirs_filter(W*xX.X, P, fs); 

    fprintf('Completed. \n\n'); 
end

%-Design space and projector matrix [pseudoinverse] for WLS
%--------------------------------------------------------------------------
fprintf('Estimate GLM parameters... \n');
xX.xKXs        = spm_sp('Set', KX);    % KWX
xX.xKXs.X      = full(xX.xKXs.X);
xX.pKX         = spm_sp('x-',xX.xKXs);                     % Projector

%-Estimate GLM parameters 
%--------------------------------------------------------------------------
beta = NaN(size(xX.X, 2), P.nch); 
res = NaN(ns, P.nch); 

beta(:, chs) = xX.pKX * KY; 
res(:, chs) = KY - xX.xKXs.X * beta(:,chs); 
ResSS = res' * res; 

%-Use non-sphericity to compute [effective] degrees of freedom
%--------------------------------------------------------------------------
V = spm_fnirs_filter(spm_fnirs_filter(W * xVi.V * W', P, fs)', P, fs); % KWVW'K'
[trRV, trRVRV] = spm_SpUtil('trRV', xX.xKXs, V);
xX.trRV        = trRV;                                     % <R'*y'*y*R>
xX.trRVRV      = trRVRV;                                   %-Satterthwaite
xX.erdf        = trRV^2/trRVRV;                            % approximation
xX.Bcov        = xX.pKX*V*xX.pKX';                      % Cov(beta)

%--------------------------------------------------------------------------
% save SPM 
SPM.xVi = xVi;               % Non-sphericity structure
SPM.xX = xX;                 % Design structure 
SPM.xCon = struct([]);      % Contrast structure 
SPM.beta = beta;           % Beta value for each chanel  
SPM.ResSS = ResSS;        % Residual sum of square for each channel 
SPM.xM = mask;      % Mask structure (channels to be analysed)

SPM.Vcorr = fullfile(SPM.swd, 'Vcorr'); 
save(SPM.Vcorr, 'V', spm_get_defaults('mat.format')); clear V; 

SPM.Vres = fullfile(SPM.swd, 'Vres'); 
save(SPM.Vres, 'res', spm_get_defaults('mat.format')); clear res; 
fprintf('Completed. \n\n');

if nargout == 0 
    fprintf('Save SPM.mat... \n'); 
    save(fullfile(SPM.swd, 'SPM.mat'), 'SPM', spm_get_defaults('mat.format')); 
    fprintf('Completed. \n\n'); 
end

fprintf('--------------------------------------------------------- \n'); 
fprintf('%-40s: %30s\n','Completed.',spm('time')); 
fprintf('--------------------------------------------------------- \n'); 



