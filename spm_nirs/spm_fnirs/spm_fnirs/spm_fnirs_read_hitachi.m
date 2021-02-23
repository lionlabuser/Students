function [y, P] = spm_fnirs_read_hitachi(F)
% Read light intensity data acquired using Hitachi ETG-4000 system,
% and convert them to variables compatible with SPM-fNIRS toolbox
% FORMAT [y,P] = spm_fnirs_read_nirscout(F)
%
% F - Name of Hitachi ETG-4000 files (*.csv)
% multiple files can be selected using cell structure 
%
% y - Light intensity [# samples x #  channels x # waves]
%
% P - information about raw data structure
%--------------------------------------------------------------------------
% P.wav      light wavelengths
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
    [F, sts] = spm_select([1 2], '^*.*\.csv$', 'Select Hitachi ETG 4000 data files');
    if ~sts, y = []; return; end
end

if iscell(F), F = cell2mat(F); end

%---------------------------------------------------------
% read measurements at wavelength 1
fid = fopen(F(1,:));
while 1
    tline = fgetl(fid);
    indx = find(tline == ','); tline(indx) = ' ';

    if strncmpi(tline, 'Age', 3) % age
        age = tline(indx(1)+1:end); 
        age(regexpi(age, 'y')) = []; 
        age = str2num(age); 
        if ~isempty(age), P.age = age; end 
    end 
    
    if ~isempty(strfind(tline, 'Wave[nm]')) % wavelength
        P.wav = str2num(tline(indx(1)+1:end));
        
    elseif ~isempty(strfind(tline, 'Sampling Period[s]')) % sampling frequency
        P.fs = 1./str2num(tline(indx(1)+1:end));
        
    elseif ~isempty(strfind(tline, 'Data')); % sample number
        tline = fgetl(fid); indx = find(tline == ','); tline(indx) = ' '; 
        
        indx_s = strfind(tline, 'Mark'); indx_s = indx_s(1);
        indx_e = find(indx == indx_s -1) + 1; % event
        
        indx_s = strfind(tline, 'PreScan');
        indx_b = find(indx == indx_s -1) + 1; % baseline
        
        nch1 = length(strfind(tline, 'CH'))./2;
        indx_w1 = 1:2:(2*nch1-1);
        indx_w2 = indx_w1+1;
        
        ns = 1;
        while 1
            tline = fgetl(fid);
            if ~ischar(tline); break; end;
            indx = find(tline == ','); tline(indx) = ' ';
            
            mes = extract(tline, indx, 2, nch1*2+1);
            y(ns,:,1) = mes(indx_w1);
            y(ns,:,2) = mes(indx_w2);
            
            event(ns, 1) = extract(tline, indx, indx_e,[]);
            base(ns, 1) = extract(tline, indx, indx_b,[]);
            
            ns = ns + 1;
        end
        fclose(fid);
        break;
    end
end
P.fname.raw.Y{1,1} = F(1,:); 
P.fname.raw.type = 'Light Intensity'; 

indx = find(base == 1); 
P.base = [min(indx) max(indx)]; 
P.mask(1, 1:nch1) = 1; 

[pathstr, name, ext] = fileparts(F(1,:)); 

%---------------------------------------------------------
% This routine was written for reading events from data. The information
% can be used in model specification of spm-fnirs toolbox 
%---------------------------------------------------------
% marker = unique(event(:,1));
% indx_z = find(marker == 0); 
% if ~isempty(indx_z), marker(indx_z) = []; end 
% for j = 1:length(marker)
%     onsets{1, j} = find(event == marker(j)); 
%     names{1, j} = sprintf('Event #%d', marker(j));
% end
% save(fullfile(pathstr, 'multiple_conditions.mat'), 'onsets', 'names', spm_get_defaults('mat.format'));
%---------------------------------------------------------

%---------------------------------------------------------
% read measurements at wavelength 2
if size(F, 1) == 2
    fid = fopen(F(2,:));
    while 1
        tline = fgetl(fid);
        if ~isempty(strfind(tline, 'Data')); % sample number
            tline = fgetl(fid); indx = find(tline == ','); tline(indx) = ' ';

            nch2 = length(strfind(tline, 'CH'))./2;
            indx_ch2 = (nch1+1):(nch1+nch2); 
            
            indx_w1 = 1:2:(2*nch2 - 1);
            indx_w2 = indx_w1+1;
            
            ns = 1;
            while 1
                tline = fgetl(fid);
                if ~ischar(tline); break; end;
                indx = find(tline == ','); tline(indx) = ' ';
                
                mes = extract(tline, indx, 2, nch2.*2+1);
                y(ns, indx_ch2, 1) = mes(indx_w1);
                y(ns, indx_ch2, 2) = mes(indx_w2);
      
                ns = ns + 1;
            end
            fclose(fid);
            break;
        end
    end
    P.mask(indx_ch2) = 1;
    P.fname.raw.Y{2,1} = F(2,:);
end

P.nch = size(y, 2);
P.ns = size(y,1); 

if nargout == 0
    % save converted data as spm-fnirs format ('NIRS.mat')
    P.fname.nirs = fullfile(pathstr, 'NIRS.mat');
    save(P.fname.nirs, 'y', 'P', spm_get_defaults('mat.format')); 
end 

function data = extract(sline, indx_comma, start_c, end_c)
if isempty(end_c), end_c = start_c; end

if end_c > length(indx_comma)
    data = str2num(sline(indx_comma(start_c-1)+1:end));
else
    data = str2num(sline(indx_comma(start_c-1)+1:indx_comma(end_c)-1));
end
