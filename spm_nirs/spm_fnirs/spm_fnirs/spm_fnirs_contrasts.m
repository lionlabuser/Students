function SPM = spm_fnirs_contrasts(SPM,Ic)
% Compute and store contrast parameters and inference SPM{.}
% FORMAT SPM = spm_fnirs_contrasts(SPM,Ic)
%
% SPM  - SPM data structure
% Ic   - indices of xCon to compute
%
%--------------------------------------------------------------------------
% note: 
% i. this script is written, based on spm_contrasts.m 
% Karl Friston, Will Penny & Guillaume Flandin
% $Id: spm_contrasts.m 5786 2013-12-06 18:25:00Z guillaume $
%
% ii. contrast x beta values for channels are saved as con_*.mat file. 
% This file will be converted to volumetric NIfTI-1 or surface-based GIfTI
% file formats, and then used for SPM group analysis. 
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

%
% $Id$

% Temporary copy of the SPM variable, to avoid saving it in SPM.mat unless
% it has changed (faster, read-only access)
%--------------------------------------------------------------------------
tmpSPM = SPM;

%-Change to results directory
%--------------------------------------------------------------------------
swd = SPM.swd; 
cd(swd);

%-Get contrast definitions (if available)
%--------------------------------------------------------------------------
try
    xCon = SPM.xCon;
catch
    xCon = [];
end

%-Set all contrasts by default
%--------------------------------------------------------------------------
if nargin < 2
    Ic   = 1:length(xCon);
end
Ic(Ic == 0) = [];

%-Load interpolated beta and variance
%--------------------------------------------------------------------------
load(SPM.Vbeta);
load(SPM.VResMS);

%-size of rendered brain image 
nviews = 6;
for i = 1:nviews
    s{i} = size(VResMS{i});
end

%-Compute & store contrast parameters, contrast/ESS images, & SPM images
%==========================================================================
spm('Pointer','Watch')
for i = 1:length(Ic)
    
    %-Canonicalise contrast structure with required fields
    %----------------------------------------------------------------------
    ic  = Ic(i);
    if isempty(xCon(ic).eidf)
        X1o           = spm_FcUtil('X1o',xCon(ic),SPM.xX.xKXs);
        load(SPM.Vcorr); 
        [trMV,trMVMV] = spm_SpUtil('trMV',X1o,V); 
        clear V; 
        xCon(ic).eidf = trMV^2/trMVMV;
    end
    
    %-Write inference SPM/PPM
    %======================================================================
    if isempty(xCon(ic).Vspm) || xCon(ic).STAT == 'P'
        
        Z = cell(1, nviews);
        switch(xCon(ic).STAT)
            case 'T'                                 %-Compute SPM{t} image
                %----------------------------------------------------------
                
                % compute interpolated t-stat 
                for k = 1:nviews
                    if sum(s{k}) ~= 0
                        cB = reshape(xCon(ic).c' * Vbeta{k}, s{k});
                        
                        Vc  = xCon(ic).c'*SPM.xX.Bcov*xCon(ic).c;
                        SE  = sqrt(VResMS{k}*Vc);
                        Z{k} = cB./SE;
                    end
                end
                
                %-----------------------------------------------------------
                % calculate contrast value for each channel, 
                % write it as con*.mat to be used in SPM group analysis 
                S = []; 
                S.cbeta = xCon(ic).c' * SPM.beta;
                
                % load channel positions 
                load(SPM.xY.VY, 'P');
                load(P.fname.pos); clear P;
    
                S.ch = R.ch;
                fname = fullfile(swd, sprintf('con_%04d.mat', ic));
                save(fname, 'S', spm_get_defaults('mat.format'));
                %-----------------------------------------------------------
                
            case 'F'                                 %-Compute SPM{F} image
                %-----------------------------------------------------------
            
                % compute ESS as sum of squared weighted beta 
                h  = spm_FcUtil('Hsqr',xCon(ic),SPM.xX.xKXs);
                for k = 1:nviews
                    if sum(s{k}) ~= 0 
                        ss = reshape(sum((h*Vbeta{k}).^2, 1), s{k});
                        MVM = ss./ trMV;
                        Z{k} = MVM./VResMS{k};
                    end
                end
        end 
        
        fname = fullfile(swd, sprintf('spm%c_%04d.mat', xCon(ic).STAT, ic));
        xCon(ic).Vspm = fname;
        
        save(fname, 'Z', spm_get_defaults('mat.format'));

    end
end

spm('Pointer','Arrow')

% place xCon back in SPM
%--------------------------------------------------------------------------
SPM.xCon = xCon;

% Check if SPM has changed. Save only if it has.
%--------------------------------------------------------------------------
if spm_check_version('matlab','8.0') >= 0, my_isequaln = @isequaln;
else my_isequaln = @isequalwithequalnans; end
if ~my_isequaln(tmpSPM,SPM)
    save(fullfile(swd, 'SPM.mat'), 'SPM', spm_get_defaults('mat.format')); 
end

