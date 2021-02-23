function [P] = spm_fnirs_specify_params(P)
% Specify parameters of the modified Beer-Lambert law 
% FORMAT [P] = spm_fnirs_specify_params(P)
%
% P - information about data structure and 
%      parameters for the modified Beer-Lambert law 
%
%--------------------------------------------------------------------------
% P.wav      light wavelengths
% P.fs            sampling frequency
% P.nch          number of channels
% P.ns            number of samples
% P.mask       mask of measurements
% P.fname     names of raw files and converted file (NIRS.mat) 
%
% P.d             distance between source and detector 
% P.acoef      molar absorption coefficients [1/(mM*cm)] 
% P.dpf         differential pathlength factor [unitless] 
% P.base       baseline period [scan] 
%__________________________________________________________________________
%
% References: 
% 1. Modified Beer-Lambert law:
% Delpy, D.T., Cope, M., Van der Zee, P., Arridge, S., Wray, S.,Wyatt, J., 1988. 
% Estimation of optical pathlength through tissue from direct time of flight measurement. 
% Phys. Med. Biol. 33 (12), 1433-1442.
%
% 2. Molar absorption coefficients: 
% Biomedical Optics Research Laboratory, UCL
% http://www.ucl.ac.uk/medphys/research/borl/intro/spectra
%
% 3. Differential pathlength factors: 
% Duncan, A., Meek, J.H., Clemence, M., Elwell, C.E., Fallon, P., Tyszczuk,
% L., Cope, M., Delpy, D.T., 1996
% Measurement of cranial optical path length as a function of age using
% phase resolved near infrared spectroscopy. 
% Pediatric research 39, 889-894.
%
% Scholkmann, F., Wolf, M., 2013. 
% General equation for the differential pathlength factor of the frontal
% head depending on wavelength and age. 
% J. Biomed. Opt. 18, 105004. 
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

%
% $Id$

if ~nargin, P = []; end

%--------------------------------------------------------------------------
% specify data information
spm_input('Specify information of fNIRS data:...', 1, 'd');

%--------------------------------------------------------------------------
str = 'sampling frequency [Hz]: ';
if ~isfield(P, 'fs')
    P.fs = spm_input(str, '+1', 'r', '', 1);
else 
    spm_input([str num2str(P.fs)], '+1', 'd'); 
end

%--------------------------------------------------------------------------
str = 'wavelengths [nm]: ';
if ~isfield(P, 'wav') 
    P.wav = spm_input(str, '+1', 'r', '', 2)';
else
    spm_input([str num2str(P.wav)], '+1', 'd'); 
end

%--------------------------------------------------------------------------
str = ['number of channels: ' num2str(P.nch)];
spm_input(str, '+1', 'd');

%--------------------------------------------------------------------------
str = 'number of measurements: '; 
if ~isfield(P, 'mask'), P.mask = ones(1, P.nch); end 
n = length(find(P.mask ~= 0)); 
spm_input([str num2str(n)], '+1', 'd'); 

%--------------------------------------------------------------------------
str = 'period of baseline [s]';
if ~isfield(P, 'base') 
    base = [0 P.ns./P.fs];
    base = round((spm_input(str, '+1', 'r', base)).*P.fs);
else 
    base = P.base; 
end

try
    if base(1) < 1, base(1) = 1; end
    if base(2) > P.ns, base(2) = P.ns; end
end
P.base = base;

%--------------------------------------------------------------------------
% specify model parameters
spm_input('Specify parameters for the modified Beer-Lambert law:...', 1, 'd');
str = 'age of subject: '; 
if ~isfield(P, 'age') 
    P.age = spm_input(str, '+1', 'n', ''); 
else 
    spm_input([str num2str(P.age)], '+1', 'd'); 
end 

str = 'distance (source-detector) [cm]: ';
if ~isfield(P, 'd') 
    P.d = spm_input(str, '+1', 'r', '');
else 
    spm_input([str num2str(P.d)], '+1', 'd'); 
end

%--------------------------------------------------------------------------
spm_input('Molar absorption coefficient [1/(mM cm)]:...', '+1', 'd');
str1 = 'HbO HbR at wavelength 1: ';
str2 = 'HbO HbR at wavelength 2: ';

if ~isfield(P, 'acoef')
    load('acoef_MCope.mat');
    wdiff = abs(acoef(:,1) - ones(size(acoef,1),1) * P.wav(1));
    indx_w1 = find(wdiff == min(wdiff));
    
    wdiff = abs(acoef(:,1) - ones(size(acoef,1),1) * P.wav(2));
    indx_w2 = find(wdiff == min(wdiff));
    
    P.acoef(1,:) = spm_input(str1, '+1', 'r', mean(acoef(indx_w1,2:3),1));
    P.acoef(2,:) = spm_input(str2, '+1', 'r', mean(acoef(indx_w2,2:3),1));
else
    spm_input([str1 num2str(P.acoef(1,:))], '+1', 'd'); 
    spm_input([str2 num2str(P.acoef(2,:))], '+1', 'd'); 
end

%--------------------------------------------------------------------------
str = 'DPF at wave 1, wave 2: ';
if ~isfield(P, 'dpf')
    dpf = 223.3 + 0.05624*(P.age^0.8493) + (-5.723) * (10^(-7)) .* (P.wav.^3) + 0.001245 .* (P.wav.^2) + (-0.9025) .* P.wav; 
    % note: dpf equation is robust between 690 and 832 nm wavelengths
    P.dpf = spm_input(str, '+1', 'r', dpf); 
else 
    spm_input([str num2str(P.dpf)], '+1', 'd'); 
end
