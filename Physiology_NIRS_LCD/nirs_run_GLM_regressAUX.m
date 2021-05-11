%___________________________________________________________________
% Copyright (C) 2019 LION Lab, Centre de recherche CHU Sainte-Justine
% www.lionlab.umontreal.ca
%___________________________________________________________________
function out = nirs_run_GLM_regressAUX(job)
% 1) Load nirs.mat
% 2) for each block, load the related AUX data (need to be
%    previously filtered + downsampled!). 
% 3) Compute a regression using the AUX data as the regressors (includes 
%    as well a constant to perform correctly the regression). 
% 4) Update the SelectedFactors.mat file (or create a new one) with the 
%    betas for each regressor + each valid block. Correct the NIRS data
%    using the input regressors (job.Covariables).
%
% INPUT:
%       job: structure with fieldname...
%           - NIRSmat, containing at least one path+file of a NIRSmat file
%                       eg. job.NIRSmat={'C:...\NIRS.mat'};
%           - covariables, containing names of AUXregressors. Need to be 
%           linked in AUX <filtered AUX> with the script 
%           <nirs_run_filterAUX.m>. If more than one variables, separate 
%           by a comma (without space). The script will automatically add 
%           a constant (of 1).
%                       eg. job.covariables='Sat,Resp';

load(job.NIRSmat{1});
%use last step of preprocessing
lst = length(NIRS.Dt.fir.pp);
rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
NC = NIRS.Cf.H.C.N; %number channels
fs = NIRS.Cf.dev.fs; %sampling rate
if ~contains(job.covariables,'Constant','IgnoreCase',true)
    cov.labels=['Constant,' job.covariables];
    cov.labels =split(cov.labels,','); %get names of covariables to regress
end

if isfield(job,'nbtimeminimum')
    nbtimeminimum=job.nbtimeminimum;
else
nbtimeminimum=10/40; %the % time non artifacted to be considered as a good channel
end
%%%%%%%%%%%
%GET AUXILIARIES
if ~isfield(NIRS.Dt,'AUX') %verify if theres an AUX
    error('No auxiliary attached to the NIRS file %s', (job.NIRSmat{1}))
end

for iAUX = 1:numel(NIRS.Dt.AUX) %Find the filtered AUX
    newAUX=numel(NIRS.Dt.AUX)+1;
    if contains(NIRS.Dt.AUX(iAUX).label,'AUX') || contains(NIRS.Dt.AUX(iAUX).label,'EEG')
        if contains(NIRS.Dt.AUX(iAUX).label,'fil')
            idAUX=iAUX;
            break
        end
    end
end
if ~exist('idAUX','var') %verify if a filtered AUX could be found
    error('No filtered AUX make sure to run filterAUX before this script')
end
goodAUX=NIRS.Dt.AUX(idAUX);

if ~isfield(goodAUX.pp(end),'sync_timesec') %check if there were synchronisation -take the last field
    error('No segmentation have been made -- ensure that aux synchronisation are ok');
end

%check if SelectedFactors file exist or create one
[nirsPATH,~,~] = fileparts(job.NIRSmat{1});
if exist(fullfile(nirsPATH,'SelectedFactors.mat'),'file')
    load(fullfile(nirsPATH,'SelectedFactors.mat'))
    Prow=length(PARCOMP)+1;
else
    PARCOMP=struct;
    Prow=1;
end

%% 
for f = 1:size(rDtp,1) %for each block
    nirsdata = fopen_NIR(rDtp{f},NC)'; %load nir data for the block
    [cPATH,cFILE,~] = fileparts(rDtp{f});
    tstart=1;
    tstop=size(nirsdata,1);
    badch = NIRS.Cf.H.C.ok(:,f);
    if badch==0 %if all channels are bad, switch to next block!
        continue
    end
    
    %%check yellow + identify channels that have a ratio of bad/total
    %%(time) that is worse than the limit (nbtimeminimum)
    mrks = [];
    ind = [];
    noise =  logical(zeros(size(nirsdata)));
    badchY=ones([size(nirsdata,2) 1]);
    [ind_dur_ch] = read_vmrk_find(fullfile(cPATH,[cFILE '.vmrk']),{'bad_step'}); %load yellow
    if ~isempty(ind_dur_ch) %if yellow
        badind = find((ind_dur_ch(:,1)+ind_dur_ch(:,2))>size(noise,1)); %find if there is longer yellow than the data duration
        if ~isempty(badind)
            disp(['Warning file ' fullfile(cPATH,[cFILE '.vmrk']) ' marker : ' num2str(badind') ' are out of range in the data file'])
            ind_dur_ch(badind,2)=size(noise,2)- ind_dur_ch(badind,1);
        end
        for Chan = 1:size(noise,2) %for each channel
            mrks = find(ind_dur_ch(:,3) == Chan); %extract its noise
            ind = ind_dur_ch(mrks,1); %extract the timestart of each bad segment
            indf = ind + ind_dur_ch(mrks,2) - 1; %calculate the timestop of each bad segment
            if ~isempty(ind)
                for i = 1:numel(ind)
                    noise(ind(i):indf(i),Chan) = 1; %mark each noise segment in the noise variable with 1
                end
            end
            if sum(ind_dur_ch([mrks],2))/size(noise,1) > nbtimeminimum %if the sum of the total noise > threshold - exclude channel
                badchY(Chan)=0;
            end
        end
    end
    
    % for NIRS bad channels, replace the data with NaN values
    idbad = find(badch==0|badchY==0); %remove exclude channel from regressor
    if ~isempty(idbad)
        nirsdata(:,idbad)=nan;
        disp(['Channels ' idbad ' are excluded from the regression'])
    end
    
    % load auxiliaries (regressors)
    try
        atstart=goodAUX.pp(end).sync_timesec{f};
    catch
        atstart = 0;
    end
    atstop = atstart+NIRS.Dt.fir.sizebloc{f}*1/fs;
    [aPATH,aFILE,aEXT]=fileparts(goodAUX.pp(end).p{f}); %current path file and extension
    [adata,ainfoBV,~,~] = fopen_EEG(goodAUX.pp(end).p{f}, atstart, atstop); %load AUX
    
    %get aux numbers that fits the covariable
    for cc=1:length(cov.labels) %create constant variable
        if strcmpi('Constant',cov.labels{cc})
            cov.data(:,cc)=ones([size(adata,1) 1]);
            cov.ConstantID=cc;
        else
            for ich=1:numel(ainfoBV.name_ele) %find the AUX selected as covariables
                if strcmpi(ainfoBV.name_ele{ich},cov.labels{cc})
                    cov.data(:,cc)=adata(:,ich); %create matrice data (X) of covariables that will be regressed
                end
            end
        end
    end
    if ~(size(cov.data,1)==size(nirsdata,1))
        error('Time axis of NIRS data and AUX data (covariables) don''t have the same size')
    end
    %iduse = find(sum(~isnan(score),2)==size(score,2)& ~isnan(MATall(:,i,j)));
    %       X = score(iduse,:);
    %      y = MATall(iduse,i,j);
    %     if ~isempty(iduse)
    %R2 statistic, the F-statistic and its p-value, and an estimate of the error variance.
    goodCH=[];
    for CH = 1:NC %for each channel
        if ~any(CH==idbad) %check if good
            goodCH=[goodCH CH];
            [tmp.b,tmp.bint,~,~,tmp.stats] = regress(nirsdata(:,CH),cov.data); %regress NIRS data with all AUX
            tmpbeta(:,CH) = tmp.b;
            tmpbetainterval{CH} = tmp.bint;
            %tmpresiduals{CH}=tmp.r;
            tmpR2(CH)=tmp.stats(1);
            tmpFstat(CH)=tmp.stats(2);
            tmpSig(CH)=tmp.stats(3);
            tmpErrorVariance(CH)=tmp.stats(4);
        else
            tmpbeta(:,CH) = NaN([size(cov.data,2) 1]);
            tmpbetainterval{CH} = NaN;
            %tmpresiduals{CH}= NaN;
            tmpR2(CH)= NaN;
            tmpFstat(CH)= NaN;
            tmpSig(CH)= NaN;
            tmpErrorVariance(CH)= NaN;
        end
    end
    tmpXcorr=nirsdata;
    for cc=1:length(cov.labels) %for each AUX
        tmpXm{cc} = tmpbeta(cc,:).* cov.data(:,cc); %multiply the beta of each AUX with the original AUX data
        if ~(cc==cov.ConstantID) %substract all regressors (AUX), except the constant
            tmpXcorr = tmpXcorr - tmpXm{cc};
        end
    end
    
    
    %write SELECTED FACTORS new info
    PARCOMP(Prow).file= f;
    PARCOMP(Prow).filestr =  sprintf('Bl%02.0f', f);
    PARCOMP(Prow).label= ['ExtPhysio_' [cov.labels{~(1:length(cov.labels)==cov.ConstantID)}] '_' PARCOMP(Prow).filestr];
    PARCOMP(Prow).type = 'GLM';
    PARCOMP(Prow).beta = tmpbeta; %beta
    PARCOMP(Prow).std = tmpErrorVariance;
    PARCOMP(Prow).AUX.data = cov.data; %données AUX
    PARCOMP(Prow).AUX.indXm = tmpXm; %données AUX * betas
    PARCOMP(Prow).AUX.label = cov.labels;
    PARCOMP(Prow).AUX.origin = goodAUX.pp(end).p{f}; %path AUX
    PARCOMP(Prow).AUX.betaInterval = tmpbetainterval;
    PARCOMP(Prow).AUX.R2 = tmpR2;
    PARCOMP(Prow).AUX.Fstat = tmpFstat;
    PARCOMP(Prow).AUX.RegressionPvalue = tmpSig; %Pvalue
    PARCOMP(Prow).data = nirsdata; %données NIRS
    PARCOMP(Prow).Xm = nirsdata - tmpXcorr; %résiduels entre les données originales et corrigées
    PARCOMP(Prow).dataCORR = tmpXcorr; %données NIRS régressées
    PARCOMP(Prow).indt = [tstart :tstop]; %indice de temps.
    PARCOMP(Prow).listgood =  goodCH;
    PARCOMP(Prow).module  = lst;
    PARCOMP(Prow).modulestr = NIRS.Dt.fir.pp(lst).pre;
    PARCOMP(Prow).ComponentToKeep = 1;
    PARCOMP(Prow).idreg = 1;
    PARCOMP(Prow).topo = tmpbeta(PARCOMP(Prow).ComponentToKeep,:);
    Prow = Prow + 1;
    clear tmp* cov.data goodCH nirsdata idbad
    
    disp(['block ' num2str(f) ' done'] );
end


fprintf('Update NIRSmat COMPLETED ...%s\n*\n**\n***\n',job.NIRSmat{1})
save(job.NIRSmat{1},'NIRS');
save(fullfile(nirsPATH,'SelectedFactors.mat'),'PARCOMP');


out.NIRSmat = job.NIRSmat;
end
