function [y, P] = spm_fnirs_read_artinis(F)
% Read optical density data acquired using Artinis Oxymon MKIII system, and
% convert them to variables compatible with SPM-fNIRS toolbox 
% FORMAT [y,P] = spm_fnirs_read_artinis(F) 
% 
% F - Name of Artinis Oxymon MKIII file (.oxy3) 
%
% y - Optical density changes [# channels x # samples x # waves] 
%
% P - information about raw data structure
%--------------------------------------------------------------------------
% P.wav        light wavelengths
% P.fs            sampling frequency
% P.nch          number of channels
% P.ns            number of samples
% P.mask       mask of measurements 
% P.fname     names of raw files and converted file (NIRS.mat) 
% P.base       baseline period [s] 
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

%
% $Id$

%---------------------------------------------------------
if ~nargin
    [F, sts] = spm_select(1, '^*.*\.oxy3$', 'Select Artinis Oxymon MKIII data files');
    if ~sts, y = []; return; end
end

if iscell(F), F = cell2mat(F); end

[pathstr, name, ext] = fileparts(F); 
%---------------------------------------------------------
% read optical density data (*.oxy3) 
hdr = read_artinis_oxy3(F); 
od = read_artinis_oxy3(F, hdr); 
auxidx = cellfun(@isempty, strfind(hdr.chantype, 'nirs'));
od(auxidx, :) = [];
%---------------------------------------------------------
% convert the data to variables compatible with SPM-fNIRS 
P.nch = size(od, 1) ./ 2; 
P.fs = hdr.Fs; 
P.ns = hdr.nSamples; 
P.mask = ones(1, P.nch); 

indx_w1 = 1:2:2*P.nch; 
indx_w2 = 2:2:2*P.nch; 

y = zeros(P.ns, P.nch, 2); 
y(:,:,1) = od(indx_w1,:)'; 
y(:,:,2) = od(indx_w2,:)'; 

P.wav = round([median(hdr.opto.wavelength(1:2:end)) median(hdr.opto.wavelength(2:2:end))]); 
P.fname.raw.Y{1,1} = F;
P.fname.raw.type = 'Optical Density'; 

% channel configuration (a pair of source and detector) 
ch_sd = []; 
[i, s_idx] = find(hdr.opto.transceiver<0);
[i, d_idx] = find(hdr.opto.transceiver>0);
s_idx = unique(s_idx);
d_idx = unique(d_idx);
for i = 1:size(hdr.opto.transceiver, 1)
       
    indx_d = find(hdr.opto.transceiver(i, :) >0);
    indx_s = find(hdr.opto.transceiver(i, :) <0);
    indx_s = find(ismember(s_idx, indx_s));
    indx_d = find(ismember(d_idx, indx_d));
    
    k = size(ch_sd, 1)+1;
    if k == 1 | ~ismember(ch_sd(:, [2 3]), [indx_s indx_d], 'rows')
      ch_sd(k, 1) = k;
      ch_sd(k, 2) = indx_s;
      ch_sd(k, 3) = indx_d;
    end
end

if max(ch_sd(:, 2))*2 ~= numel(hdr.opto.wavelength)
  error('Only exactly 2 wavelengths per transmitter are allowed - please contact Artinis support (support@artinis.com)');
end

%--------------------------------------------------------------------------
% write channel configuration as text file
% note: this file can be used in spatil preprocessing
% of spm-fnirs toolbox
%--------------------------------------------------------------------------
sfname = fullfile(pathstr, 'ch_config.txt');
fid_w = fopen(sfname, 'w');
fprintf(fid_w, '%s\r\n', 'Ch, Source, Detector');
for i = 1:P.nch, fprintf(fid_w, '%s\r\n', [num2str(ch_sd(i,1)) ', ' num2str(ch_sd(i,2)) ', ' num2str(ch_sd(i,3))]); end
fclose(fid_w);


if nargout == 0
    % save converted data as spm-fnirs format ('NIRS.mat')
    P.fname.nirs = fullfile(pathstr, 'NIRS.mat');
    save(P.fname.nirs, 'y', 'P', spm_get_defaults('mat.format')); 
end 

