function [N, interpY] = spm_fnirs_2dtopo(S, fname, n)
% Project 3D MNI positions of fNIRS sensors onto 2D topography, and
% interpolate channel-wise contrasts, and then save results as NIfTI format
% FORMAT [N] = spm_fnirs_2dtopo(S, fname, n) 
%
% S             structure array of contrast x beta values for each channel 
% fname     name of NIfTI file to be saved 
% n             spatial resolution of 2D topography 
%
% N            structure array of NIfTI file 
%
%--------------------------------------------------------------------------
% note: 
% S.cbeta   contrast x beta 
% S.ch         channel positions in MNI space 
% see spm_fnirs_contrasts.m for more details. 
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging
 
% $Id$

% default spatial resolutions of 2D topography 
if nargin < 3, n = 64; end 

if isfield(S.ch, 'xyzH')
    xyz_ch = S.ch.xyzH'; 
else 
    xyz_ch = S.ch.xyzC'; 
end

Y = S.cbeta(S.ch.label); 

indx = find(isnan(Y) == 1); 
Y(indx) = []; 
xyz_ch(indx, :) = []; 

[sdir sname ext] = fileparts(fname);
fname = fullfile(sdir, [sname '.nii']); 

%-----------------------------------
% load NIRS sensor positions (3D MNI coordinates) 
mdir = fileparts(which('spm_fnirs')); 
load(fullfile(mdir, 'gen_conimg', 'eeg_128sens.mat')); 

nch = size(xyz_ch, 1); 
nch_eeg = size(sens.chanpos, 1); 

sens.chanpos = [sens.chanpos; xyz_ch]; 
sens.elecpos = [sens.elecpos; xyz_ch];
    
for i = 1:nch 
    sens.chantype{nch_eeg+i,1} = 'eeg';
    sens.chanunit{nch_eeg+i,1} = 'V';
    sens.label{nch_eeg+i,1} = ['CH' num2str(i)];
end

%-----------------------------------
% project 3D positions onto 2D topography map
[xy, label] = spm_eeg_project3D(sens, 'EEG');

% Locate channels and generate mask for converting M/EEG data into images
[Cel, x, y] = spm_fnirs_locate_channels(xy, n);

%-----------------------------------
% perform 2D interpolation on scattered data 
interpY = NaN(n,n);

Cel(1:nch_eeg,:) = [];

% calculate average pixels between channels
nch = size(Cel,1);
for i = 1:nch 
    indx = 1:nch; indx(i) = []; 
    mat_dist = ones(nch-1,1) * Cel(i,:) - Cel(indx,:);
    dist_t = sqrt(sum(mat_dist.^2,2)); 
    dist(i) = min(dist_t); 
end
dist = mean(dist); 

% find pixels to be interpolated 
xy_interp = []; 
for i = 1: nch 
    [xy_ch, ~] = spm_fnirs_vert_mindist(Cel(i,:), [x y], dist); 
    xy_interp = [xy_interp; xy_ch]; 
end
xy_interp = unique(xy_interp, 'rows'); 

interpY = NaN(n,n); 
interpY(sub2ind([n n], xy_interp(:,1), xy_interp(:,2))) = griddata(Cel(:,1), Cel(:,2), double(Y), xy_interp(:,1), xy_interp(:,2),'linear');

%-----------------------------------
% save data as NIFTI format 
% this script is based on spm_eeg_convert2images.m 
% 'scalp' option
N     = nifti;
N.mat = eye(4);
N.mat_intent = 'Aligned';

C     = [68  100];  % origin
V     = [136 172]/n; % voxel size

avflag = [0 1 1];

N.mat = [...
    V(1)  0     0               -C(1);...
    0     V(2)  0               -C(2);...
    0     0     1               1;...
    0     0     0               1];
N.mat(3,4) = N.mat(3,4) - N.mat(3,3);

dat = file_array('', [n n], 'FLOAT32-LE');

cdat       = dat;
cdat.fname = fname;
cdat.dim   = dat.dim; 
N.dat      = cdat;
create(N);

N.dat(:,:)     = interpY;  

