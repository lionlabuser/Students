function varargout = spm_fnirs_dct(X, fs, cutoff)
% Removes low frequency confounds from X
% FORMAT [Y] = spm_nirs_dct(X, P)
%
% X - data matrix [# samples x # channels]
% fs - sampling frequency 
% cutoff - cutoff periods [s] 
%
% varargout - filtered data 
%
% Note: This script is based on spm_filter.m 
% Karl Friston
% $Id: spm_filter.m 4498 2011-09-23 18:40:43Z guillaume $
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

%
% $Id$

%-Configure filter
%==========================================================================
k = size(X, 1); 
RT = 1/fs; 
n       = fix(2*(k*RT)/cutoff + 1);
X0      = spm_dctmtx(k,n);
X0 = X0(:,2:end);

Y = X - X0 * (X0' * X);

if nargout == 1
    varargout = {Y}; 
elseif nargout == 2
    varargout = {Y, X0}; 
end