 function [g, interpY] = spm_fnirs_3dscalp(S, fname, vert_interp)
% Project 3D MNI positions of fNIRS sensors onto 3D scalp surface
% perform spherical spline interpolation, and save them as GIfTI file format.
%
% FORMAT [g, interpY] = spm_fnirs_3dscalp(S, vert_interp)
%
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging
 
% $Id$

if isfield(S.ch, 'xyzH')
    xyz_ch = S.ch.xyzH'; 
else 
    xyz_ch = S.ch.xyzC'; 
end

Y = S.cbeta(S.ch.label)'; 
indx = find(isnan(Y) == 1); 
Y(indx) = []; 
xyz_ch(indx, :) = []; 

[sdir sname ext] = fileparts(fname);
fname = fullfile(sdir, [sname '.gii']); 

%-----------------------------------
% load vertices from 3D scalp mesh 
fsurf = 'scalp_2562.surf.gii'; 
Mh = gifti(fullfile(spm('Dir'), 'canonical', fsurf)); % read scalp mesh 
nch = size(xyz_ch, 1); 

if nargin < 3 
    % find vertices closest to fNIRS sensor positions, and to be interpolated
    xyz_vch = zeros(nch, 3);
    r = 30; % unit: mm
    [vert_interp, xyz_interp] = spm_fnirs_roivert(xyz_ch, r);
else
    xyz_interp = Mh.vertices(vert_interp,:); 
end

for i = 1:nch, xyz_vch(i,:) = spm_fnirs_vert_mindist(xyz_ch(i,:), Mh.vertices); end

% perform spline spherical interpolation 
cur_dir = cd;
fieldtrip_dir = fullfile(spm('Dir'), 'external','fieldtrip', 'private');
cd(fieldtrip_dir);

W =sphericalSplineInterpolate(xyz_vch', xyz_interp', 1e-5, 4);
cd(cur_dir);

nvert = size(Mh.vertices, 1); 
interpY = zeros(nvert, 1); 

interpY(vert_interp, 1) = W * Y; 

% save data as GIFTI format 
Mh.cdata = double(interpY); 

% spm_mesh_render(Mh); % for visualisation 
g = gifti(Mh);
g.private.metadata(1).name = 'SurfaceID';
g.private.metadata(1).value = fsurf;
save(g, fname, 'ExternalFileBinary');

% copy canonical gifti file (scalp_2562.surf.gii) to the directory 
listing = struct2cell(dir(sdir)); 
fnames = cell2mat(listing(1,:)); 

if isempty(strfind(fnames, 'scalp_2562.surf.gii'))
    copyfile(fullfile(spm('Dir'), 'canonical', 'scalp_2562.surf.gii'), sdir);
end

    

