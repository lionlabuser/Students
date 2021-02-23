function [xy_i] = spm_fnirs_xy_interp(xy_s, rend)
% Identify voxels of rendered brain template to be interpolated 
% FORMAT [xy_i] = spm_fnirs_xy_interp(xy_s, rend) 
%
% xy_s      channel positions on the rendered brain surfaces 
% rend      structure array of rendered brain surface 
%
% xy_i      indices of pixels to be interpolated 
%
%--------------------------------------------------------------------------
% note: 
% data format of xy_s and rend are consistent with 
% R.ch.xy, and R.rend, which are outputs of spm_fnirs_spatialpreproc_ui.m 
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% 
% $Id$

%--------------------------------------------------------------------------
% calculate average of pixels between channels 
% (average of distance between channels)

Cel = xy_s{2}; 
indx = find(sum(Cel) ~= 0); 
Cel = Cel(:, indx); 

nch = size(Cel,2);
for i = 1:nch 
    indx = 1:nch; indx(i) = []; 
    mat_dist = Cel(:,i) * ones(1, nch-1) - Cel(:, indx);
    dist(i) = min(sqrt(sum(mat_dist.^2,1))); 
end
% dist_th = mean(dist); 
dist_th = max(dist); 

%--------------------------------------------------------------------------
% find pixels to be interpolated 

nviews = 6; % views of rendered brain 
for i = 1:nviews
    img = rend{i}.ren; 
    indx = find(sum(xy_s{i}) ~= 0); 
    Cel = xy_s{i}(:, indx); 
    
    nch = size(Cel, 2); 
    xy = []; 
    if  nch > 2
        [x y] = find(abs(img - 0.04) > 0.01); % pixels of brain mask 
        
        nvox = length(x); 
        for j = 1:nch 
            dist = Cel(:, j)*ones(1, nvox) - [x y]'; 
            dist = sqrt(sum(dist.^2,1)); 
            indx_r = find(dist < dist_th); 
            xy_ch = [x(indx_r) y(indx_r)]'; 
            xy = [xy xy_ch]; 
        end
        xy = unique(xy', 'rows')'; 
    end
    xy_i{i} = xy; 
end 
