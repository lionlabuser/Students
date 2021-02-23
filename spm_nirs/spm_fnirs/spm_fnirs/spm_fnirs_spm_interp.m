function [SPM] = spm_fnirs_spm_interp(SPM)
% Interpolate GLM parameters over surfaces of rendered brain 
% FORMAT [SPM] = spm_fnirs_spm_interp(SPM) 
%
% SPM       GLM parameters for measurements of fNIRS channels 
%
%--------------------------------------------------------------------------
% In this code, 
% interpolated GLM parameters are saved as *.mat files, 
% whose names are saved in fields of SPM structure: 
%
% Vbeta.mat         Interpolated betas 
% VResMS.mat      Interpolated variance scaled by tr(RV)
% VRpv.mat          Resels per voxel 
% mask.mat           Mask for interpolated data 
%
% SPM.Vbeta         Filehandle - Beta 
% SPM.VResMS      Filehandle - ResMS 
% SPM.xVol.VRpv  Filehandle - Rpv (resels per voxel) 
% SPM.VM             Filehandle - mask 
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% 
% $Id$

fprintf('--------------------------------------------------------- \n'); 
fprintf('Interpolate GLM parameters over brain surfaces... \n'); 
fprintf('--------------------------------------------------------- \n'); 

%--------------------------------------------------------------------------
% load SPM file 
if ~nargin, 
    [SPM sts] = spm_select(1, '^SPM*.*.\.mat$', 'Select a SPM file for fNIRS');
    if ~sts, SPM = []; return; end
end

if ischar(SPM), load(SPM); end

%--------------------------------------------------------------------------
% load GLM parameters for measurements of fNIRS channels 
fprintf('Load channel-wise parameters... \n'); 

beta = SPM.beta;
ResSS = SPM.ResSS; 
ns = SPM.nscan;
erdf = spm_SpUtil('trRV',SPM.xX.xKXs);  % error df

load(SPM.Vres); 
nSres = min(ns, spm_get_defaults('stats.maxres'));
indx_res = round(linspace(1, ns, nSres))';
res = res(indx_res, :);

% load channel positions
load(SPM.xY.VY, 'P');
load(P.fname.pos); clear P;

xy = R.ch.xy; rend = R.rend; chs = R.ch.label; clear R; 

mask = SPM.xM; 
ch_nroi = find(mask == 0 | mask == -1); 
for i = 1:6, 
    xy{i}(:, ch_nroi) = 0; 
    s{i} = size(rend{i}.ren); 
end
    
% identify voxels to be interpolated
[xy_i] = spm_fnirs_xy_interp(xy, rend);
clear rend; 

fprintf('Completed. \n\n'); 

%--------------------------------------------------------------------------
% interpolate GLM parameters over surfaces of rendered brain

fprintf('Interpolate parameters & estimate smoothness... \n'); 
% initialize parameters
nviews = 6;
Vbeta = cell(1, nviews);
VResMS = cell(1, nviews); 
VM = cell(1, nviews); % mask

imethod = 'natural'; % interpolation method

FWHM = cell(1, nviews); % estimated FWHM
VRpv = cell(1, nviews); % resels per voxel
Rc = cell(1, nviews); % vector of resel counts

% interpolate parameters over each view of brain surface (eg, dorsal view) 
for i = 1:nviews
    if ~isempty(xy_i{i})
        indx = find(sum(xy{i}) ~= 0);
        xy_c = xy{i}(:, indx); 
        
        indx_c = chs(indx);
        
        nc = size(beta, 1); % number of contrasts
        nvox = s{i}(1) * s{i}(2);
        
        nch_v = size(xy_c, 2); % number of channels on a view of brain
        nvox_i = size(xy_i, 2); % number of interpolated voxels 
        
        % estimate interpolation matrix 
        emtx = eye(nch_v); 
        kern = zeros(nch_v, nvox); 
        for j = 1:nch_v 
            kern(j, sub2ind(s{i}, xy_i{i}(1,:), xy_i{i}(2,:))) = griddata(xy_c(1,:), xy_c(2,:), emtx(:, j), xy_i{i}(1,:), xy_i{i}(2,:), imethod);
        end 
        indx_nan = find(isnan(kern) == 1); 
        kern(indx_nan) = 0; kern = sparse(kern); 
        
        % compute interpolated beta 
        Vbeta{i} = sparse(beta(:, indx_c) * kern); 
        
        % get mask 
        mask = sparse(s{i}(1), s{i}(2));
        indx_m = find(Vbeta{i}(1,:) ~= 0);
        mask(indx_m) = 1;
        VM{i} = mask;
        S{i} = nnz(mask); % number of non-zero elements
        clear mask;
        
        % compute interpolated residual SSQ 
        [y x] = meshgrid(indx_c, indx_c); 
        indx_s = sub2ind(size(ResSS), x, y); 
        iResSS = sparse(s{i}(1), s{i}(2));
        %iResSS(indx_m) = diag(kern(:, indx_m)' * ResSS(indx_s) * kern(:, indx_m)); 
        
        max_size = 10000; % 763 MB
        num_block = ceil(length(indx_m)/max_size);
        for block = 1:num_block
            ind1 = (block-1) * max_size + 1;
            ind2 = min(block*max_size,length(indx_m));
            bindx_m = indx_m( ind1:ind2 );
            iResSS(bindx_m) = diag(kern(:,bindx_m)' * ResSS(indx_s) * kern(:,bindx_m));
        end
        
        VResMS{i} = iResSS./ SPM.xX.trRV; 
        
        % interpolate standardised residuals
        VResI = NaN(nSres, nvox);
        for j = 1:nSres
            iRes = reshape(res(j, indx_c) * kern, s{i}); 
            VResI(j,indx_m) = iRes(indx_m)./sqrt(iResSS(indx_m)/erdf); 
        end 
        clear iRes kern; 
        
        % estimate smoothness based on residual images 
        [FWHM{i}, VRpv{i}, Rc{i}] = spm_fnirs_est_smoothness(VResI, full(VM{i}), [ns erdf]);
    end
end

fprintf('Completed. \n\n'); 
%--------------------------------------------------------------------------
% save estimated (interpolated) parameters
SPM.Vbeta = fullfile(SPM.swd, 'Vbeta');
save(SPM.Vbeta, 'Vbeta', spm_get_defaults('mat.format'));

SPM.VResMS = fullfile(SPM.swd, 'VResMS');
save(SPM.VResMS, 'VResMS', spm_get_defaults('mat.format'));

SPM.xVol.VRpv = fullfile(SPM.swd, 'VRpv');
save(SPM.xVol.VRpv, 'VRpv', spm_get_defaults('mat.format'));

SPM.VM = fullfile(SPM.swd, 'mask');
save(SPM.VM, 'VM', spm_get_defaults('mat.format'));

SPM.xVol.FWHM = FWHM; % Smoothness data
SPM.xVol.R = Rc; % Resel counts
SPM.xVol.S = S; % non-zero voxels

delete(deblank([SPM.Vres '.mat'])); 
SPM = rmfield(SPM, 'Vres'); 

fprintf('Save SPM.mat... \n');
save(fullfile(SPM.swd, 'SPM.mat'), 'SPM', spm_get_defaults('mat.format'));
fprintf('Completed. \n\n');

fprintf('--------------------------------------------------------- \n'); 
fprintf('%-40s: %30s\n','Completed.',spm('time')); 
fprintf('--------------------------------------------------------- \n'); 

