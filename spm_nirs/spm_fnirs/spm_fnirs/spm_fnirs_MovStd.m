function [y1] = MovStd_spm_fnirs(x,k)
% Function to calculate the moving standard deviation.
%
% INPUT
% x:    input signal
% k:    half of the centered window length (w = 1 + 2*k)
%__________________________________________________________________________
% Dr. Felix Scholkmann, Biomedical Optics Research Laboratory (BORL), 
% Universtiy Hospital Zurich, University of Zurich, Zurich, Switzerland
% Felix.Scholkmann@usz.ch
% Version 1: 30 September 2008. This version: 29 May 2015
%_________________________________________________________________________

% (1) Calculate the MSD
y1 = NaN(length(x),1); % preallocate y1
for i = k+1:length(x)-k
    y1(i) = std(x(i-k:i+k));
end