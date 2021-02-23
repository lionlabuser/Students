% Produce individual contrast images containing the experimental effects of
% interest on a 2D regular grid (canonical scalp surface). 
%
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

%
% $Id$

spm('Defaults', 'EEG')
[fname, sts] = spm_select(Inf, '^con*.*.\.mat$', 'Select contrast beta files');
nsubj = size(fname, 1);

fwhm = 14; % full width at half maximum of Gaussian kernel
for i = 1:nsubj 
    fname_c = deblank(fname(i,:)); 
    load(fname_c); 
    [N, interpY] = spm_fnirs_2dtopo(S, fname_c); 
    
    indx_m = find(isnan(interpY) == 0); 
    spm_smooth(N.dat.fname,N.dat.fname, [fwhm fwhm fwhm]); 
    
    vol = spm_vol(N.dat.fname); 
    siData = spm_read_vols(vol); 
    sY = zeros(size(siData)); 
    sY(indx_m) = siData(indx_m); 
    
    vol = spm_write_vol(vol,sY); 
end

    
    
