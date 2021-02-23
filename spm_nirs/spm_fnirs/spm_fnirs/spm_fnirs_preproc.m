function [fY, P] = spm_fnirs_preproc(Y, P)
% Apply temporal filters to fNIRS time series (eg, hemoglobin changes)
% FORMAT [fY, P] = spm_fnirs_preproc(Y,P)
%
% Y     matrix of fNIRS time series [# samples x # channels] 
%        or [# samples x # channels x # hemoglobins (eg, HbO, HbR, HbT)
%
% P     structure array of filter parameters (P.K)
%
% fY    matrix of filtered data
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

%
% $Id$

dim = size(Y);
if ndims(Y) == 3,
    Y = reshape(Y, [dim(1) dim(2) * dim(3)]);
else
    dim(3) = 1;
end

n = size(Y, 2);

%--------------------------------------------------------------------------
% i. motion artifact correction
if strcmpi(P.K.M.type, 'MARA')
    M = P.K.M;
    sfname = fullfile(spm_file(P.fname.nirs, 'path'), 'motion_params.mat'); 
    indx_m = []; % indices of measurements to be corrected
    for i = 1:dim(3), indx_m = [indx_m M.chs+dim(2)*(i-1)]; end
    
    % L: moving window length
    if ~iscell(M.L) % L: scalar
        if isscalar(M.L)
            L = NaN(1, n);
            L(indx_m) = M.L;
            L = mat2cell(L, 1, dim(2) * ones(1, dim(3)));
            M.L = L;
            save(sfname, 'L', '-append', spm_get_defaults('mat.format')); 
        else
            fprintf('Error: parameter L should be scalar or cell array.\n');
        end
    end
    L = round(cell2mat(M.L) .* P.fs);
    
    % alpha: smoothing factor-motion artifact
    if ~iscell(M.alpha) % L: scalar
        if isscalar(M.alpha)
            alpha = NaN(1, n);
            alpha(indx_m) = M.alpha;
            alpha = mat2cell(alpha, 1, dim(2) * ones(1, dim(3))); 
            M.alpha = alpha; 
            save(sfname, 'alpha', '-append', spm_get_defaults('mat.format'));
        else
            fprintf('Error: parameter alpha should be scalar or cell array.\n');
        end
    end
    alpha = cell2mat(M.alpha);
    
    % threshold for motion detection
    mstd_y = NaN(3, n);
    
    spm_input('Remove motion artifact from Hb changes:', 1, 'd');
    nd = size(indx_m, 2);
    spm_progress_bar('Init', nd, 'Std estimation', 'Total number of data');
    
    for i = 1:nd 
        std_y = spm_fnirs_MovStd(Y(:, indx_m(i)), round(L(indx_m(i))./2)); 
        indx_n = find(isnan(std_y) == 1); std_y(indx_n) = [];
        mstd_y(1, indx_m(i)) = min(std_y);
        mstd_y(2, indx_m(i)) = mean(std_y);
        mstd_y(3, indx_m(i)) = max(std_y);
        spm_progress_bar('Set', i);
    end
    spm_progress_bar('Clear');
    
    if ~iscell(M.th)
        if isscalar(M.th)
            th = M.th * mstd_y(2,:);
            th = mat2cell(th, 1, dim(2) * ones(1, dim(3)));
            M.th = th; 
            save(sfname, 'th', '-append', spm_get_defaults('mat.format'));
        else
            fprintf('Error: parameter th should be scalar or cell array.\n');
        end
    end
    th = cell2mat(M.th);
    
    % apply MARA method
     spm_progress_bar('Init', nd, 'motion artifact removal', 'Total number of data');
    for i = 1:nd 
        if th(indx_m(i)) < mstd_y(3, indx_m(i)) && th(indx_m(i)) > mstd_y(1, indx_m(i))
            Y(:, indx_m(i)) = spm_fnirs_MARA(Y(:, indx_m(i)), P.fs, th(indx_m(i)), L(indx_m(i)) , alpha(indx_m(i)));
        end
        spm_progress_bar('Set', i); 
    end
    spm_progress_bar('Clear'); 
    
    % update structure array for MARA
    P.K.M = M;
end

%--------------------------------------------------------------------------
% ii. physiological noise removal (band-stop filter)
if strcmpi(P.K.C.type, 'Band-stop filter')
    addpath(fullfile(spm('Dir'), 'external', 'fieldtrip', 'preproc'));
    
    forder = 5;
    ftype = 'but';
    nf = size(P.K.C.cutoff, 1);
    fdir = 'twopass';
    
    for j = 1:nf
        try
            Y = ft_preproc_bandstopfilter(Y', P.fs, P.K.C.cutoff(j,:), forder, ftype, fdir)';
        catch
            Y = ft_preproc_bandstopfilter(Y', P.fs, P.K.C.cutoff(j,:), [], 'fir', fdir)';
        end
    end
end

%--------------------------------------------------------------------------
% iii. change sampling rate
if strcmpi(P.K.D.type, 'yes')
    p = round(10*P.K.D.nfs);
    q = round(10*P.fs);
    Y = resample(Y, p, q);
    
    P.K.D.nfs = p/q*P.fs;
    P.K.D.ns = size(Y, 1);
    
    dim(1) = P.K.D.ns;
end

fY = reshape(Y, dim);

