function [xyz_v, indx_v]  = spm_fnirs_vert_mindist(pos, vertices, r)
% Find a vertex closest to sensor position 
%
% FORMAT [xyz_v, indx_v] = spm_fnirs_vert_mindist(pos, vertices, r)
% 
% pos         sensor position [nch x 3] 
% vertices   xyz coordintes of all vertices in the mesh brain [nvert x 3] 
% r             radius [mm, or pixels]
%
% xyz_v     xyz coordinates of vertices closest to sensor positions
% indx_v    indices of vertices closest to sensor positions 
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

%
% $Id$

if nargin == 2, r = 0; end

nvert = size(vertices, 1); 
dist = ones(nvert, 1) * pos - vertices; 
dist = sqrt(sum(dist.^2, 2)); 

[dist, indx] = sort(dist); 

if r == 0
    indx_v = indx(1);
else
    indx_v = indx(find(dist < r));
end

xyz_v = vertices(indx_v, :); 
