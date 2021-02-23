function spm_fnirs_specify1st_ui(F)
% Specify general linear model (GLM) of fNIRS for the 1st level analysis
% FORMAT [SPM] = spm_fnirs_specify1st_ui(F)
%
% F         'NIRS.mat' file to be analysed
%
% Parameters for 1st level SPM analysis of HbO, HbR, and HbT will be saved
% as SPM.mat file in a directory of ..\HbO, \HbR, \and \HbT, respectively.
%
%--------------------------------------------------------------------------
% This script was written, based on the following scripts:
% (i) spm_fMRI_design.m
% Karl Friston
% $Id: spm_fMRI_design.m 5183 2013-01-10 15:30:20Z ged $
%
% (ii) spm_fmri_spm_ui.m
% Karl Friston
% $Id: spm_fmri_spm_ui.m 6088 2014-07-03 17:57:09Z guillaume $
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

%
% $Id$

SPM = [];

fprintf('--------------------------------------------------------- \n');
fprintf('Specify GLM-fNIRS for the 1st level analysis...\n');
fprintf('--------------------------------------------------------- \n');

%--------------------------------------------------------------------------
% specify NIRS.mat file

spm_input('Specify parameters for 1st level GLM analysis of fNIRS:', 1, 'd');

if ~nargin, % specify file name of NIRS.mat
    [F sts] = spm_select(1, '^NIRS*.*.\.mat$', 'Select a fNIRS file (NIRS.mat)');
    if ~sts, return; end
end

if ischar(F), load(F); end

%--------------------------------------------------------------------------
% specify directory to save SPM.mat file
[sdir, sts] = spm_select(1, 'dir', 'Select a SPM directory');

%--------------------------------------------------------------------------
% specify chromophores (HbO, HbR, HbT) to be analysed

Finter = spm_figure('GetWin','Interactive');
bcolor = get(Finter,'color');
WS     = spm('WinScale');
dx     = 60;

spm_input('Specify 1st level analysis for fNIRS data:...  ',1,'d');
spm_input('Hb signal to be analysed? ', '+1','d');

% radio button for hb changes
str{1,1} = 'HbO';
str{2,1} = 'HbR';
str{3,1} = 'HbT';
for i = 1:3
    h1(i) = uicontrol(Finter,...
        'Position',[193+dx*(i-1) 325 020 020].*WS,...
        'BackgroundColor',bcolor,...
        'Style','radiobutton');
    
    h2(i)  = uicontrol(Finter,'String',str{i,1},...
        'Style','text',...
        'FontSize',10,...
        'BackgroundColor',bcolor,...
        'HorizontalAlignment','right',...
        'Position',[180+dx*(i-1) 345 030 020].*WS);
end

h3 = uicontrol(Finter,'String','done','Position', [300 300 060 020].*WS,...
    'Callback', 'uiresume(gcbf)');

uiwait(Finter);

str_nirs = [];
a = cell2mat(get(h1, 'Value')); indx = find(a == 1);
str_nirs = str(indx); % chromophores to be analysed

delete(findobj(get(Finter,'Children'),'flat'));

fprintf('Hemoglobin signal to be analysed: %s \n\n', cell2mat(str_nirs'));
%--------------------------------------------------------------------------
fprintf('Specify experimental design...\n');

str = 'specify design in';
SPM.xBF.UNITS = spm_input(str,'+1','scans|secs');

str = 'specify conditions using a file?';
answer = spm_input(str, '+1', 'yes|no');

switch answer
    case 'yes'
        [fname_u sts] = spm_select(1, 'mat', 'Select a file of multiple condition');
        load(fname_u); % names = cell{1,n}, onsets, durations
        for i = 1:size(names, 2)
            U(i).name = names(1,i);
            U(i).ons = onsets{1,i}(:);
            U(i).dur = durations{1,i}(:);
        end
    case 'no'
        U = spm_fnirs_get_ons_ui;
end

for i = 1:numel(U)
    U(i).P.name = 'none';
    U(i).P.h = 0;
end
SPM.Sess.U = U;

% in case of resampling
switch P.K.D.type
    case 'yes'
        fs = P.K.D.nfs;
        ns = P.K.D.ns;
    case 'no'
        fs = P.fs;
        ns = P.ns;
end

SPM.xY.RT = 1/fs;
SPM.nscan = ns;

try
    SPM.xBF.T  = spm_get_defaults('stats.fmri.t');
catch
    SPM.xBF.T  = 16;
end

try
    SPM.xBF.T0 = spm_get_defaults('stats.fmri.t0');
catch
    SPM.xBF.T0 = 8;
end
SPM.xBF.dt = SPM.xY.RT/SPM.xBF.T;
SPM.xBF.Volterra = 1; % 1st Volterra expansion
SPM.Sess.U = spm_get_ons(SPM,1);

fprintf('Completed. \n\n');

%--------------------------------------------------------------------------
% covariates - C
fprintf('Specify other regressors...\n');

spm_input('Other regressors (eg, systemic confounds)',1,'d');
C     = [];
c     = spm_input('user specified','+1','w1',0);
while size(C,2) < c
    str = sprintf('regressor %i',size(C,2) + 1);
    C  = [C spm_input(str,2,'e',[],[ns Inf])];
end

Cname = cell(1,size(C,2));
for i = 1:size(C,2)
    str      = sprintf('regressor %i',i);
    Cname{i} = spm_input('name of','+0','s',str);
end

SPM.Sess.C.C = C;
SPM.Sess.C.name = Cname;

fprintf('Completed. \n\n');
%--------------------------------------------------------------------------
%-Intrinsic autocorrelations (Vi) for non-sphericity ReML estimation

str   = 'Correct for serial correlations?';
if strcmp(P.K.L.type, 'no')
    cVi   = {'none','AR(1)'};
    cVi   = spm_input(str,'+1','b',cVi);
else
    cVi = 'none';
    spm_input([str ' ' cVi], '+1', 'd');
end
fprintf('Correct for serial correlation: %s \n\n', cVi);

% create Vi structure
switch cVi
    case 'none'
        SPM.xVi.V  = speye(sum(ns));
        cVi = 'i.i.d';
    case 'AR(1)'  % assume AR(0.2) in xVi.Vi
        SPM.xVi.Vi = spm_Ce(ns,0.2);
        cVi = 'AR(0.2)';
end

SPM.xVi.form = cVi;

%--------------------------------------------------------------------------
%  generate design matrix using specified parameters
SPM.xY.VY = P.fname.nirs; % file name of Y

SPM0 = SPM;
save_SPM = 0;
for i = 1:size(str_nirs, 1)
    fprintf('Generate design matrix for %s...\n', str_nirs{i,1});
    
    [SPM] = spm_fMRI_design(SPM0, save_SPM);
    SPM.xY.type = str_nirs{i,1}; % hemoglobin type (eg, HbO or HbR)
    
    if strcmpi(str_nirs{i,1}, 'HbR')
        if strcmpi(SPM.xBF.name, 'hrf') || strcmpi(SPM.xBF.name, 'hrf (with time derivative)') || strcmpi(SPM.xBF.name, 'hrf (with time and dispersion derivatives)')
            SPM.xX.X(:,1:end-1) = SPM.xX.X(:,1:end-1).*(-1);
        end
    end
    
    fprintf('Completed. \n\n');
    
    %-Design description - for saving and display
    %--------------------------------------------------------------------------
    ntr = length(SPM.Sess.U);
    Bstr = sprintf('[%s] %s', str_nirs{i,1}, SPM.xBF.name);
    Hstr = P.K.H.type; 
    if strcmpi(Hstr, 'DCT'), Hstr = sprintf('%s, Cutoff: %d {s}', Hstr, P.K.H.cutoff); end 
    Lstr = P.K.L.type;
    if strcmpi(Lstr, 'Gaussian'), Lstr = sprintf('%s, FWHM %d', Lstr, P.K.L.fwhm); end
    
    SPM.xsDes = struct(...
        'Basis_functions',      Bstr,...
        'Number_of_sessions',   sprintf('%d',1),...
        'Trials_per_session',   sprintf('%-3d',ntr),...
        'Interscan_interval',   sprintf('%0.2f {s}',SPM.xY.RT),...
        'High_pass_Filter',     Hstr,...
        'Low_pass_Filter', Lstr);
    
    spm_DesRep('DesMtx',SPM.xX,[],SPM.xsDes);
    
    swd = fullfile(sdir, str_nirs{i,1});
    if ~isdir(swd), mkdir(swd), end
    
    SPM.swd = swd;
    
    fprintf('Save SPM.mat... \n');
    save(fullfile(SPM.swd, 'SPM.mat'), 'SPM', spm_get_defaults('mat.format'));
    fprintf('Completed. \n\n');
    
    %--------------------------------------------------------------------------
    % Delete files related to previous SPM.mat file (if exist)
    fname = spm_select('FPList', SPM.swd, '^con_.*\.mat$');
    if ~isempty(fname),
        for i = 1:size(fname, 1), delete(deblank(fname(i,:))); end
    end
    
    fname = spm_select('FPList', SPM.swd, '^spmT_.*\.mat$');
    if ~isempty(fname),
        for i = 1:size(fname, 1), delete(deblank(fname(i,:))); end
    end
    
    fname = spm_select('FPList', SPM.swd, '^spmF_.*\.mat$');
    if ~isempty(fname),
        for i = 1:size(fname, 1), delete(deblank(fname(i,:))); end
    end
end

fprintf('--------------------------------------------------------- \n');
fprintf('%-40s: %30s\n','Completed.',spm('time'));
fprintf('--------------------------------------------------------- \n');

