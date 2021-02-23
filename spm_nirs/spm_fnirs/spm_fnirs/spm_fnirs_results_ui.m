function spm_fnirs_results_ui(SPM)
% Compute specified and thresholded SPM 
% FORMAT spm_fnirs_results_ui(SPM) 
%
% SPM       structure array of GLM parameters 
%
%--------------------------------------------------------------------------
% In this code, 
% (i) estimated GLM parameters are interpolated over surfaces of rendered
% brain, using spm_fnirs_spm_interp(SPM). 
% (ii) contrast parameters and inference SPM are computed, using
% spm_fnirs_contrasts.m which is a function in spm_fnirs_viewer_stat.m 
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

%
% $Id$

%--------------------------------------------------------------------------
% load SPM.mat file 
if ~nargin, 
    [SPM sts] = spm_select(1, '^SPM*.*.\.mat$', 'Select a SPM file for fNIRS');
    if ~sts, SPM = []; return; end
end

if ischar(SPM), load(SPM); end

%--------------------------------------------------------------------------
% Interpolate GLM parameters over surfaces of rendered brain 
if ~isfield(SPM, 'Vbeta') || ~isfield(SPM, 'VResMS') || ~isfield(SPM.xVol, 'VRpv')  
    [SPM] = spm_fnirs_spm_interp(SPM);
end

%--------------------------------------------------------------------------
% Specify contrast vector, compute inference SPM, and then visualize
% thresholded SPM on a surface of rendered brain template 
spm_fnirs_viewer_stat(SPM);
