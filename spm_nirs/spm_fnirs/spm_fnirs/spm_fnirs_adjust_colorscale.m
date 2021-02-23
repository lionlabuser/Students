function [simg] = spm_fnirs_adjust_colorscale(oimg, min_c, max_c) 
% Adjust color scale of image to user-defined range 
% FORMAT simage = spm_fnirs_adjust_colorscale(oimg, min_c, max_c) 
%
% oimg                   original image 
% [min_c max_c]     lowest and highest values in the color scale 
%
% simg                    image after adjustment of color scale 
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% 
% $Id$

indx_over = find(oimg > min_c | oimg == min_c); 
oimg = oimg(indx_over); 
simg = 63./(max_c - min_c) .* (oimg - min_c) + 65; 

    