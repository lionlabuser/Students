function [Y, P] = spm_fnirs_convert_ui(F, P)
% Convert fNIRS data to optical density and hemoglobin changes 
% FORMAT [Y, P] = spm_fnirs_convert_ui(F, P)
%
% F - name of fNIRS files (.txt; .mat) or matrix of measurments (y)
% P - information about data structure and 
%      parameters for the modified Beer-Lambert law 
%
% Y - structure array of fNIRS data 
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
%--------------------------------------------------------------------------
% Input F can be 
% (i) names of text files which include measurements at two wavelengths 
% (ii) a name of mat file which includes measurement matrix (y)
% (iii) measurement matrix (y), where y is [# samples x # channels]. 
% 
%--------------------------------------------------------------------------
% Examples: 
% >> spm_fnirs_convert_ui; 
% all input variables will be specified using GUI. 
%
% >> F{1,1} = 'meas1.txt'; F{2,1} = 'meas2.txt'; 
% >> spm_fnirs_convert_ui(F); 
% Model parameters, P will be specified using GUI 
%
% >> F = 'NIRS.mat'; 
% >> spm_fnirs_convert_ui(F); 
%
% if meaurements and parameters are saved as matrix y, and structure P, 
% >> spm_fnirs_convert_ui(y, P); 
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

%
% $Id$

fprintf('--------------------------------------------------------- \n'); 
fprintf('Convert fNIRS data to optical density and hemoglobin changes...\n'); 
fprintf('--------------------------------------------------------- \n'); 

%---------------------------------------------------------
if ~nargin 
    [F, sts] = spm_select([1 2], 'any', 'Select fNIRS file(s) to be used for SPM analysis');
    if ~sts, Y =[]; return; end 
elseif nargin == 1
    P = []; 
end

%---------------------------------------------------------
spm_input('Read raw data and convert it to hemoglobin changes:...', 1, 'd');
fprintf('Read raw data...\n'); 

y = []; % matrix of measurements 
if ischar(F) % F: file names 
    if iscell(F), F = cell2mat(F); end 
    
    fname = deblank(F(1,:)); 
    [sdir, name, ext] = fileparts(fname); 
    
    if strcmpi(ext, '.txt') || strcmpi(ext, '.csv') % raw data (text format) 
        for i =1:2, 
            fname = deblank(F(i,:));
            [y(:,:,i), P.ns, P.nch] = spm_fnirs_read_txt(fname); 
            P.fname.raw.Y{i,1} = fname; 
        end
    elseif strcmpi(ext, '.mat') % mat file of Y matrix 
        load(fname); P.fname.raw.Y = fname; 
    end 
else
    y = F; sdir = []; 
end

fprintf('Completed. \n\n');

%---------------------------------------------------------
% specify type of measurements 
str = 'Measurement type? '; 
if ~isfield(P.fname.raw, 'type')
    [P.fname.raw.type, sts] = spm_input(str, '+1', 'Light Intensity|Optical Density');
else
    spm_input([str P.fname.raw.type], '+1', 'd');
end

switch P.fname.raw.type 
    case 'Light Intensity' 
        fprintf('Calculate optical density changes...\n'); 
        [Y.od P] = spm_fnirs_calc_od(y, P);
        fprintf('Completed. \n\n');
        
    case 'Optical Density'
        Y.od = y; 
        P.base = []; 
end

%---------------------------------------------------------
% specify parameters for measurements Y
fprintf('Specify parameters for the modified Beer-Lambert law...\n');
P = spm_fnirs_specify_params(P);
fprintf('Completed. \n\n');

%---------------------------------------------------------
% calculate hemoglobin changes using optical density
fprintf('Calculate hemoglobin concentration changes...\n');
[Y.hbo, Y.hbr, Y.hbt] = spm_fnirs_calc_hb(Y, P);
fprintf('Completed. \n\n');

if nargout == 0
    fprintf('Save NIRS.mat... \n'); 
    if isempty(sdir), 
        [sdir, sts] = spm_select(1, 'dir', 'Select a directory where NIRS.mat is saved.');
        if ~sts, Y = []; return; end 
    end
    P.fname.nirs = fullfile(sdir, 'NIRS.mat');
    save(P.fname.nirs, 'Y', 'P', spm_get_defaults('mat.format')); 
    fprintf('Completed. \n\n'); 
end 

fprintf('--------------------------------------------------------- \n'); 
fprintf('%-40s: %30s\n','Completed.',spm('time')); 
fprintf('--------------------------------------------------------- \n'); 
