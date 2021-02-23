function [vert_interp, xyz_interp] = spm_fnirs_roivert(xyz_ch, r)
% Find vertices within regions of interest 
%
% FORMAT [vert_interp, xyz_interp] = spm_nirs_roivert(xyz_ch, r)
%
% xyz_ch        MNI coordinates of channels [nch x 3]
% r                  radius of roi of each channel 
%
% xyz_interp  MNI coordinates of vertices within ROI
%
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging
 
% $Id$
if nargin < 2, r = 30; end 

%-----------------------------------
% load vertices from 3D scalp mesh 
fsurf = 'scalp_2562.surf.gii'; 
Mh = gifti(fullfile(spm('Dir'), 'canonical', fsurf)); % read scalp mesh 
nch = size(xyz_ch, 1); 

% find vertices closest to fNIRS sensor positions, and to be interpolated
xyz_vch = zeros(nch, 3); 

xyz_interp = []; 
vert_interp = []; 
for i = 1:nch,
    [xyz, vert] = spm_fnirs_vert_mindist(xyz_ch(i,:), Mh.vertices, r); 
    xyz_vch(i,:) = xyz(1,:); 
    vert_ch(i,:) = vert(1);
    
    xyz_interp = [xyz_interp; xyz];
    vert_interp = [vert_interp; vert]; 
end
[vert_interp indx] = unique(vert_interp);
xyz_interp = xyz_interp(indx, :); 
