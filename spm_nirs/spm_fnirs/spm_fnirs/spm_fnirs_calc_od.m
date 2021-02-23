function [y, P] = spm_fnirs_calc_od(X, P)
% Calculate optical density changes from light intensity 
% FORMAT [y,P] = spm_fnirs_calc_od(X, P)
%
% X   light intensity 
% P   structure array of model parameters 
%
% y   optical density changes 
%
%__________________________________________________________________________
%
% References: 
% Chapter 2 in Mark Cope, PhD thesis (1991): 
% The application of near infrared spectroscopy to non invasive monitoring
% of cerebral oxygenation in the newborn infant. 
% Department of Medical Physics and Engineering, UCL, 342. 
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

%
% $Id$

if nargin == 1 || ~isfield(P, 'base')
    P.base = [1 size(X, 1)];
end

% calculate optical density changes using light intensity 
ns = size(X, 1); nch = size(X, 2); 
X = reshape(abs(X), [ns nch*2]); 

[indx] = find(X == 0); X(indx) = realmin('single'); 
y = -log(X) + log(ones(ns, 1) * mean(X(P.base(1):P.base(2), :), 1)); 

y = reshape(y, [ns nch 2]); 

