function [fY, KL] = spm_fnirs_filter(Y, P, fs)
% Apply temporal filters to hemoglobin changes 
% FORMAT [fY, KL] = spm_fnirs_filter_hb(Y,P,fs)
%
% Y     matrix of hemoglobin changes [# samples x (# channels x 3)] 
% P     structure array of filter parameters (P.K)
% fs    sampling frequency 
%
% fY    matrix of filtered data 
% KL    temporal smoothing matrix
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

%
% $Id$

if nargin < 3, fs = P.fs; end 

dim = size(Y); 
if ndims(Y) == 3, Y = reshape(Y, [dim(1) dim(2)*dim(3)]); end

n = size(Y, 2); 

%--------------------------------------------------------------------------
% i. detrending based on DCT 
if strcmpi(P.K.H.type, 'DCT') 
    Y = spm_fnirs_dct(Y, fs, P.K.H.cutoff);
end

%--------------------------------------------------------------------------
% ii. temporal smoothing 
switch P.K.L.type 
    case 'Gaussian' 
        [Y, KL] = spm_fnirs_tsmooth(Y, fs, P.K.L.type, P.K.L.fwhm); 
    case 'HRF'
        [Y, KL] = spm_fnirs_tsmooth(Y, fs, P.K.L.type); 
    case 'no'
        KL = []; 
end

%--------------------------------------------------------------------------
% iii. temporal smoothing 
fY = reshape(Y, dim); 

