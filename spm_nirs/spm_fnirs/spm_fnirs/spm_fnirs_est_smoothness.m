function [FWHM,d,R] = spm_fnirs_est_smoothness(V,VM,ndf)
% Estimate smoothness based on residual images 
% FORMAT [FWHM,VRpv,R] = spm_est_smoothness(V,VM,[ndf])
%
% V     -  standardized residual images [# scans x # voxels] 
% VM    - mapped mask image 
% ndf   - A 2-vector, [n df], the original n & dof of the linear model
%
% FWHM  - estimated FWHM in all image directions
% VRpv  - handle of Resels per Voxel image
% R     - vector of resel counts
%--------------------------------------------------------------------------
% note: this script is same as spm_est_smoothness.m 
% Stefan Kiebel, Tom Nichols
% $Id: spm_est_smoothness.m 6157 2014-09-05 18:17:54Z guillaume $
%
% Input formats of the function are changed from 
% filename of NIfTI-1 (V and VM)  
% to matrices of V and VM 
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

%
% $Id$

%-Assign input arguments
%--------------------------------------------------------------------------
dim = size(VM); dim(3) = 1; 

ns = size(V, 1); % number of scans 
V = reshape(V', [dim ns]); 

n_full = ndf(1);
edf    = ndf(2);

%-Dimensionality of image
%--------------------------------------------------------------------------
D        = 3 - sum(dim(1:3) == 1);
if D == 0
    FWHM = [Inf Inf Inf];
    R    = [0 0 0];
    return;
end

%-Find voxels within mask
%--------------------------------------------------------------------------
d = VM; 
[Ix,Iy,Iz] = ndgrid(1:dim(1),1:dim(2),1:dim(3));
Ix = Ix(d~=0); Iy = Iy(d~=0); Iz = Iz(d~=0);

%-Compute covariance of derivatives
%--------------------------------------------------------------------------
str   = 'Spatial non-sphericity (over scans)';
fprintf('%-40s: %30s',str,'...estimating derivatives');                 %-#
spm_progress_bar('Init',100,'smoothness estimation','');

L     = zeros(size(Ix,1),D,D);
ssq   = zeros(size(Ix,1),1);
for i = 1:ns
    [d, dx, dy, dz] = spm_sample_vol(V(:,:,i), Ix,Iy,Iz,1); 
    
    % sum of squares
    %----------------------------------------------------------------------
    ssq  = ssq + d.^2;
    
    % covariance of finite differences
    %----------------------------------------------------------------------
    if D >= 1
        L(:,1,1) = L(:,1,1) + dx.*dx;
    end
    if D >= 2
        L(:,1,2) = L(:,1,2) + dx.*dy;
        L(:,2,2) = L(:,2,2) + dy.*dy;
    end
    if D >= 3
        L(:,1,3) = L(:,1,3) + dx.*dz;
        L(:,2,3) = L(:,2,3) + dy.*dz;
        L(:,3,3) = L(:,3,3) + dz.*dz;
    end
    
    spm_progress_bar('Set',100*i/length(V));
    
end
spm_progress_bar('Clear')

L  = L/size(V,3);
L  = L*(n_full/edf);  % Scale


%-Evaluate determinant (and xyz components for FWHM)
%--------------------------------------------------------------------------
if D == 1
    resel_xyz = L;
    resel_img = L;
end
if D == 2
    resel_xyz = [L(:,1,1) L(:,2,2)];
    resel_img = L(:,1,1).*L(:,2,2) - ...
                L(:,1,2).*L(:,1,2);
end
if D == 3
    resel_xyz = [L(:,1,1) L(:,2,2)  L(:,3,3)];
    resel_img = L(:,1,1).*L(:,2,2).*L(:,3,3) + ...
                L(:,1,2).*L(:,2,3).*L(:,1,3)*2 - ...
                L(:,1,1).*L(:,2,3).*L(:,2,3) - ...
                L(:,1,2).*L(:,1,2).*L(:,3,3) - ...
                L(:,1,3).*L(:,2,2).*L(:,1,3);
end    
resel_img(resel_img<0) = 0;
% Convert det(Lambda) and diag(Lambda) to units of resels
resel_img = sqrt(resel_img/(4*log(2))^D);
resel_xyz = sqrt(resel_xyz/(4*log(2)));


%-Optional mask-weighted smoothing of RPV image
%--------------------------------------------------------------------------
if spm_get_defaults('stats.rft.nonstat')
    fwhm_vox = 3;
else
    fwhm_vox = 0;
end

%-Write Resels Per Voxel image
%--------------------------------------------------------------------------
fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...writing resels/voxel image');%-#

for i = 1:dim(3)
    d = NaN(dim(1:2));
    I = find(Iz == i);
    if ~isempty(I)
        d(sub2ind(dim(1:2), Ix(I), Iy(I))) = resel_img(I);
    end
end


%-(unbiased) RESEL estimator and Global equivalent FWHM
% where we desire FWHM with components proportional to 1./mean(resel_xyz),
% but scaled so prod(1./FWHM) agrees with (the unbiased) mean(resel_img).
%--------------------------------------------------------------------------
i     = isnan(ssq) | ssq < sqrt(eps);
resel_img = mean(resel_img(~i,:));
resel_xyz = mean(resel_xyz(~i,:));

RESEL = resel_img^(1/D)*(resel_xyz/prod(resel_xyz)^(1/D));
FWHM  = full(sparse(1,1:D,1./RESEL,1,3));
FWHM(isnan(FWHM)) = 0;
FWHM(~FWHM) = 1;

%-resel counts for search volume (defined by mask)
%--------------------------------------------------------------------------
% R0   = spm_resels_vol(VM,[1 1 1])';
% R    = R0.*(resel.^([0:3]/3));
% OR
R      = spm_resels_vol(VM,FWHM)';

fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done');               %-#
