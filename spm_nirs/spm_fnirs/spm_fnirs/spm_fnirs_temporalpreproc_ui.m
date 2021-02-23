function [P] = spm_fnirs_temporalpreproc_ui(F)
% Apply temporal filters to fNIRS data 
% FORMAT [P] = spm_fnirs_temporalpreproc_ui(F)
%
% F     'NIRS.mat' file to be analysed 
%
% P     structure array of filter parameters (P.K) 
%
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

%
% $Id$

fprintf('--------------------------------------------------------- \n'); 
fprintf('Apply temporal filters to fNIRS data...\n'); 
fprintf('--------------------------------------------------------- \n'); 

if ~nargin, % specify file name of NIRS.mat 
    [F sts] = spm_select(1, '^NIRS*.*.\.mat$', 'Select fNIRS file for temporal preprocessing');
    if ~sts, P = []; return; end 
end

spm_input('Apply temporal preprocessing to fNIRS data:', 1, 'd'); 

%--------------------------------------------------------------------------
fprintf('Read and display time series of fNIRS data...\n'); 

load(F); K = [];

% identify channels of interest 
if isfield(P.fname, 'pos') 
    load(P.fname.pos);
    mask = zeros(1, R.ch.nch);
    indx = find(R.ch.mask == 1);
    mask(R.ch.label(indx)) = 1; clear R;
else 
    mask = ones(1, P.nch); 
end

mask = mask .* P.mask;
ch_roi = find(mask ~= 0); 

% display time series of fNIRS data
spm_fnirs_viewer_timeseries(Y, P, [], ch_roi);

%==================================
% specify parameters for temporal preprocessing
%==================================
fprintf('Specify parameters for temporal preprocessing...\n'); 

%--------------------------------------------------------------------------
str = 'Motion artifact correction? ';
if ~isfield(K, 'M')
    K.M.type = spm_input(str, '+1', 'MARA|no');
    if strcmpi(K.M.type, 'MARA')
        str = 'specify parameters using a file?';
        answer = spm_input(str, '+1', 'yes|no');
        switch answer
            case 'yes'
                [fname_p sts] = spm_select(1, 'mat', 'Select a file of motion correction parameters');
                load(fname_p);
            case 'no'
                str = 'channels to be analysed:';
                gui_s = spm_input(str, '+1', 'b', {'All', 'Selected'}, [0 1], 0); 
                if ~gui_s % all channels 
                    chs = ch_roi; 
                else 
                    chs = spm_input('channels: ', '+1', 'n', ch_roi); 
                end
                
                str = 'moving window length [sec]';
                L = spm_input(str, '+1', 'r', '1', 1);
                
                str = 'threshold factor-motion detection';
                th = spm_input(str, '+1', 'r', '3', 1);
                
                str = 'smoothing factor-motion artifact';
                alpha = spm_input(str, '+1', 'r', '5', 1);
                
                sfname = fullfile(spm_file(P.fname.nirs, 'path'), 'motion_params.mat');
                save(sfname, 'chs', spm_get_defaults('mat.format'));
        end
        K.M.chs = chs;
        K.M.L = L;
        K.M.th = th;
        K.M.alpha = alpha;
    end
else
    spm_input([str K.M.type], '+1', 'd');
end
    
%--------------------------------------------------------------------------
str = 'Physiological noise removal? ';
if ~isfield(K, 'C')
    K.C.type = spm_input(str, '+2', 'Band-stop filter|no');
    if strcmpi(K.C.type, 'Band-stop filter')
        K.C.cutoff = spm_input('stopband frequencies Hz [start end]:', '+1', 'r', '[0.12 0.35; 0.7 1.5]');
    end
else
    spm_input([str K.C.type], '+1', 'd');
end

%--------------------------------------------------------------------------
str = sprintf('Change sampling rate from %3.2fHz?', P.fs);
if ~isfield(K, 'D')
    K.D.type = spm_input(str, '+2', 'yes|no');
else
    spm_input([str K.D.type], '+1','d');
end

if strcmpi(K.D.type, 'yes')
    fs = spm_input('new sampling rate [Hz]:', '+1', 'r', 1);
    K.D.nfs = fs;
else
    fs = P.fs; 
end

%--------------------------------------------------------------------------
str = 'Detrending? ';
if ~isfield(K, 'H') 
    K.H.type = spm_input(str, 1, 'DCT|no');
    if strcmpi(K.H.type, 'DCT')
        K.H.cutoff = spm_input('cut-off period [sec]:', '+1', 'r', '128');
    end
else
    spm_input([str K.H.type], '+1','d'); 
end
    
%--------------------------------------------------------------------------
str = 'Temporal smoothing? ';
if ~isfield(K, 'L')
    K.L.type = spm_input(str, '+2', 'Gaussian|HRF|no');
    if strcmpi(K.L.type, 'Gaussian')
        K.L.fwhm = spm_input('FWHM [sec]','+1','r','4',1);
    end
else 
    spm_input([str K.L.type], '+1','d'); 
end

fprintf('Completed. \n\n');

%==================================
% apply filters to hemoglobin changes
%================================== 
fprintf('Apply a temporal filter to hemoglobin changes... \n'); 

load(F); % mask of measurements might be changed using spm_fnirs_viewer_timeseries 
P.K = K; % update structure P 

% change structure array to matrix 
y = spm_vec(rmfield(Y, 'od')); 
y = reshape(y, [P.ns P.nch 3]); 

[fy, P] = spm_fnirs_preproc(y, P); 
[fy] = spm_fnirs_filter(fy, P, fs); 

fprintf('Completed. \n\n');

%--------------------------------------------------------------------------
% save the results (parameters)
if nargout ==0, 
    fprintf('Update header of NIRS.mat file...\n'); 
    save(P.fname.nirs, 'P', '-append', spm_get_defaults('mat.format')); 
    fprintf('Completed. \n\n');
end 

%--------------------------------------------------------------------------
fprintf('Display results...\n'); 
spm_fnirs_viewer_timeseries(y, P, fy, ch_roi);

fprintf('--------------------------------------------------------- \n'); 
fprintf('%-40s: %30s\n','Completed.', spm('time')); 
fprintf('--------------------------------------------------------- \n'); 

