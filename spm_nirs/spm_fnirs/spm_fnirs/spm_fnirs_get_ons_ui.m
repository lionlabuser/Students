function [U] = spm_fnirs_get_ons_ui
% This function was modified for SPM-fNIRS

% Returns input [designed effects] structures
% FORMAT [U] = spm_get_ons(SPM,s)
%
% SPM   - SPM structure (see spm_fMRI_design.m)
% s     - session number
%
% U     - (1 x n)   struct array of (n) trial-specific structures
%
%   U(i).name   - cell of names for each input or cause
%   U(i).u      - inputs or stimulus function matrix
%   U(i).dt     - time bin (seconds)
%   U(i).ons    - onsets    (in SPM.xBF.UNITS)
%   U(i).dur    - durations (in SPM.xBF.UNITS)
%   U(i).P      - parameter struct.
%
%       U(i).P(p).name - parameter name
%       U(i).P(p).P    - parameter vector
%       U(i).P(p).h    - order of polynomial expansion
%       U(i).P(p).i    - sub-indices of u pertaining to P
%__________________________________________________________________________
%
% SLICE TIMING
%
% With longs TRs you may want to shift the regressors so that they are
% aligned to a particular slice. This is effected by resetting the
% values of defaults.stats.fmri.t and defaults.stats.fmri.t0 in
% spm_defaults.
% defaults.stats.fmri.t is the number of time-bins per scan used when
% building regressors. Onsets are defined in temporal units of scans
% starting at 0.
% defaults.stats.fmri.t0 is the first time-bin at which the regressors are
% resampled to coincide with data acquisition. If defaults.stats.fmri.t0
% is set to 1 then the regressors will be appropriate for the first slice.
% If you want to temporally realign the regressors so that they match
% responses in the middle slice then make defaults.stats.fmri.t0 =
% defaults.stats.fmri.t/2 (assuming there is a negligible gap between
% volume acquisitions).
% Default values are defaults.stats.fmri.t=16 and defaults.stats.fmri.t0=1.
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_get_ons.m 4083 2010-10-08 10:31:55Z guillaume $

%-prompt string
%----------------------------------------------------------------------
s = 1;
str = sprintf('Session %d: trial specification',s);
spm_input(str,1,'d')

U   = {};
v   = spm_input('number of conditions/trials',2,'w1');

% get trials
%--------------------------------------------------------------------------
for i = 1:v
    
    % get names
    %----------------------------------------------------------------------
    try
        Uname     = U(i).name(1);
    catch
        str       = sprintf('name for condition/trial %d ?',i);
        Uname     = {spm_input(str,3,'s',sprintf('trial %d',i))};
        U(i).name = Uname;
    end
    
    % get main [trial] effects
    %======================================================================
    
    % onsets
    %----------------------------------------------------------------------
    try
        ons = U(i).ons;
        ons = ons(:);
    catch
        ons = [];
    end
    if isempty(ons)
        str      = ['vector of onsets - ' Uname{1}];
        ons      = spm_input(str,4,'r',' ',[Inf 1]);
        U(i).ons = ons(:);
    end
    
    % durations
    %----------------------------------------------------------------------
    try
        dur = U(i).dur;
        dur = dur(:);
    catch
        dur = [];
    end
    if isempty(dur)
        str = 'duration[s] (events = 0)';
        while 1
            dur = spm_input(str,5,'r',' ',[Inf 1]);
            if length(dur) == 1
                dur    = dur*ones(size(ons));
            end
            if length(dur) == length(ons), break, end
            str = sprintf('enter a scalar or [%d] vector',...
                length(ons));
        end
        U(i).dur = dur;
    end
end
