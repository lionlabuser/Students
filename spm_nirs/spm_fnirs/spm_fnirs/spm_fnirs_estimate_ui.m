function spm_fnirs_estimate_ui
% Estimate GLM parameters from fNIRS measurements
% FORMAT spm_fnirs_estimate_ui 
%
% See spm_fnirs_spm.m for more details. 
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

%
% $Id$

% specify file name of SPM.mat
[SPM sts] = spm_select([1 Inf], '^SPM*.*.\.mat$', 'Select SPM file for fNIRS');
if ~sts, return; end

% estimate parameters for SPM analysis of fNIRS data 
for i = 1:size(SPM, 1)
    spm_fnirs_spm(deblank(SPM(i,:))); 
end
