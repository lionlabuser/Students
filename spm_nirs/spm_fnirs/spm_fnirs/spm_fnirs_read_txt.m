function [y, ns, nch] = spm_fnirs_read_txt(fname)
% Read data from text file (*.txt or *.csv) 
% FORMAT [y, ns, nch] = spm_fnirs_read_txt(fname)
%
% fname - name of fNIRS files (.txt; .csv) 
% 
% y - data matrix [# samples x # channels] 
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

%
% $Id$


% read csv or txt format 
y = []; ns = 0; 

fid = fopen(fname); 
while 1 
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    ns = ns + 1;
    
    indx_c = find(tline == ',');
    tline(indx_c) = ' ';
    y(ns, :) = str2num(tline);
end
fclose(fid); 

nch = size(y, 2); 

