function varargout = spm_fnirs_tsmooth(X, fs, Ltype, fwhm)
% Removes high frequency confounds from Y
% FORMAT [Y] = spm_fnirs_tsmooth(X, P)
%
% X - data matrix [# samples x # channels]
% fs - sampling frequency 
% Ltype - Gaussian or HRF 
% fwhm - FWHM of Gaussian kernel 
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
k = size(X,1); 
RT = 1/fs; 

switch Ltype % low-pass filter
    case 'Gaussian'
        sigma   = fwhm/RT;
        h       = round(4*sigma);
        h       = exp(-[-h:h].^2/(2*sigma^2));
        n       = length(h);
        d       = [1:n] - (n + 1)/2;
        if      n == 1, h = 1; end
    case 'HRF'
        h = spm_hrf(RT);
        h = [h; zeros(size(h))];
        g = abs(fft(h));
        h = real(ifft(g));
        h = fftshift(h)';
        n = length(h);
        d = [1:n] - n/2 -1;
end
KL = spdiags(ones(k,1)*h, d, k,k);
KL = spdiags(1./sum(KL')',0,k,k)*KL;

Y = KL * X;

if nargout == 1
    varargout = {Y}; 
elseif nargout == 2
    varargout = {Y, KL}; 
end
