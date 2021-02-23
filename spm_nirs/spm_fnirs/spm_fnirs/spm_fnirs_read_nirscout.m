function [y, P] = spm_fnirs_read_nirscout(F)
% Read light intensity data acquired using NIRScout system,
% and convert them to variables compatible with SPM-fNIRS toolbox
% FORMAT [y,P] = spm_fnirs_read_nirscout(F)
%
% F - Name of NIRScout files (*.wl1, *.wl2, *.hdr)
%
% y - Light intensity [# samples x # channels x # waves]
%
% P - information about raw data structure
%--------------------------------------------------------------------------
% P.fs                   sampling frequency
% P.nch                number of channels
% P.ns                  number of samples
% P.ch_mask        mask of channels
% P.fname           names of raw files and converted file (NIRS.mat) 
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

%
% $Id$

if ~nargin
    [F, sts] = spm_select([1 Inf], 'any', 'Select NIRx NIRScout data files (wl1, wl2, hdr)');
    if ~sts, y = []; return; end
end

if iscell(F), F = cell2mat(F); end

nfiles = size(F, 1);
% read *.wl1 file
for i = 1:nfiles
    fname = deblank(F(i,:)); 
    
    [pathstr,name,ext] = fileparts(fname);
    if strcmpi(ext, '.wl1')
        % read measurements at wavelength 1
        [y(:,:,1)] = spm_fnirs_read_txt(fname);
    elseif strcmpi(ext, '.wl2')
        % read measurements at wavelength 2
        [y(:,:,2)] = spm_fnirs_read_txt(fname);
    elseif strcmpi(ext, '.hdr')
        % read header information
        fid = fopen(fname);
        strcomp = {};
        strcomp{1} = 'Sources';
        strcomp{2} = 'Detectors';
        strcomp{3} = 'SamplingRate';
        strcomp{4} = 'Events';
        strcomp{5} = 'S-D-Key';
        strcomp{6} = 'S-D-Mask';
        
        var = {}; event = [];
        n = 1;
        
        while 1
            tline = fgetl(fid);
            if ~ischar(tline); break, end;
            flag = regexpi(tline, strcomp{n});
            if isempty(flag) == 0
                if n < 4 % # source, # dector, sampling rate
                    var{n} = str2num(tline(1,flag+length(strcomp{n})+1:end));
                    n = n +1;
                    
                elseif n == 4 % event
                    flag2 = [];
                    tline = fgetl(fid);
                    count = 1;
                    while isempty(flag2) == 1
                        tmp = str2num(tline);
                        event(count,:) = tmp(1:2);
                        count = count + 1;
                        tline = fgetl(fid);
                        flag2 = strfind(tline, '#');
                    end
                    
                    marker = unique(event(:,2));
                    for j = 1:length(marker)
                        indx = find(event(:,2) == marker(j));
                        onset = event(indx,1); 
                        onsets{1, j} = onset(:)'; 
                        names{1,j} = sprintf('Event #%d', marker(j)); 
                    end
                    save(fullfile(pathstr, 'multiple_conditions.mat'), 'onsets', 'names', spm_get_defaults('mat.format')); 
                    n = n + 1;
                    
                elseif n == 5
                    % channel configuration (a pair of source and detector)
                    indx_se = strfind(tline, '"');
                    indx_comma = strfind(tline, ',');
                    indx_colon = strfind(tline, ':');
                    
                    str_ch = {};
                    nch = length(indx_colon);
                    str_ch{1} = tline(1,indx_se(1)+1:indx_colon(1)-1);
                    str_ch{nch} = tline(1,indx_comma(end-1)+1:indx_colon(end)-1);
                    for j = 2:nch-1, str_ch{j} = tline(1,indx_comma(j-1)+1:indx_colon(j)-1); end
                    
                    ch_sd = [];
                    for j = 1:nch
                        str_chconf = str_ch{j};
                        indx = strfind(str_chconf, '-');
                        num_s = str2num(str_chconf(1:indx-1));
                        num_d = str2num(str_chconf(indx+1:end));
                        ch_sd(j, 1) = num_s; % source
                        ch_sd(j, 2) = num_d; % detector
                    end
                    n = n + 1;
                    
                elseif n == 6
                    % Mask of channels of interest
                    flag2 = [];
                    sdmask = [];
                    tline = fgetl(fid);
                    while isempty(flag2) == 1
                        sdmask = [sdmask; str2num(tline)];
                        tline = fgetl(fid);
                        flag2 = strfind(tline, '#');
                    end
                    try
                        sdmask = reshape(sdmask', [1 var{1}*var{2}]);
                    catch
                        sdmask = sdmask(:)';
                    end
                    break;
                end
            end
        end
        fclose(fid);
        
        P.fs = var{3}; 
    end
    P.fname.raw.Y{i,1} = fname; 
end

chs = find(sdmask == 1);  
y = y(:, chs, :); 
ch_sd = ch_sd(chs, :); 

P.nch = size(y, 2); 
P.ns = size(y, 1); 
P.mask = ones(1, P.nch); 
P.fname.raw.type = 'Light Intensity'; 

%--------------------------------------------------------------------------
% write channel configuration as text file
% note: this file can be used in spatial preprocessing
% of spm-fnirs toolbox
%--------------------------------------------------------------------------
sfname = fullfile(pathstr, 'ch_config.txt');
if ~exist(sfname, 'file')
    fid_w = fopen(sfname, 'w');
    fprintf(fid_w, '%s\r\n', 'Ch, Source, Detector');
    for j = 1:P.nch, fprintf(fid_w, '%s\r\n', [num2str(j) ', ' num2str(ch_sd(j,1)) ', ' num2str(ch_sd(j,2))]); end
    fclose(fid_w);
end

if nargout == 0
    % save converted data as spm-fnirs format ('NIRS.mat')
    P.fname.nirs = fullfile(pathstr, 'NIRS.mat');
    save(P.fname.nirs, 'y', 'P', spm_get_defaults('mat.format')); 
end 

