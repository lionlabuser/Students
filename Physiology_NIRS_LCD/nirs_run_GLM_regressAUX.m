%___________________________________________________________________
% Copyright (C) 2021 LCD for lionlab, Centre de recherche CHU Sainte-Justine
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

nbchminimum =  0.37; %verification for each block: minimum % of good channels to keep the block
%min 20 channels out of 54

if isfield(job,'VIEWplot')
    VIEWplot=job.VIEWplot;
else
    VIEWplot=0;
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
IDrows=[]; %Id of rows for the current script in the PARCOMP - will be use later to create a graph output regarding the betas

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
        nirsdata(:,idbad) = nan;
        disp(['Channels ' num2str(idbad') ' are excluded from the regression'])
    end
    
    % verification of minimum number of channels
    currentNchan=size(nirsdata(1,1:NC/2),2);
    if currentNchan < nbchminimum*NC/2
        fprintf('Not enough valid channels (%d HbO channels) to continue with the spatial filtering.\nBlock %s is skipped.\n',...
            currentNchan,f);
        continue
    end
    
    % load auxiliaries (regressors)
    try
        atstart=goodAUX.pp(end).sync_timesec{f};
    catch
        atstart = 0;
        disp('Warning! Synchronization time not found. Automatically set to 0')
    end
    atstop = atstart+NIRS.Dt.fir.sizebloc{f}*1/fs;
    [aPATH,aFILE,aEXT] = fileparts(goodAUX.pp(end).p{f}); %current path file and extension
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
    
    %if the block is bad (no good channels), then skip to next one
    if isempty(goodCH)
        continue
    end

    %Regression
    tmpXcorr = nirsdata;
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
    
    IDrows=[IDrows Prow];
    Prow = Prow + 1;
    clear tmp* cov.data goodCH nirsdata idbad
    
    disp(['block ' num2str(f) ' done'] );
end

%% figure of beta distribution
figg=figure('units','normalized','outerposition',[0 0 1 1]);
figg=tiledlayout(2,length(cov.labels)-1,'TileSpacing','Compact','Padding','Compact');
maintitle='Distribution of AUX beta values for physiology regression on HbO channels, sorted by channels and blocks';
smalltitle='Median + interquartile range';
%get betas for hbo
for f=1:length(IDrows)
    for cc=1:length(cov.labels)
        figbetas.(cov.labels{cc})(f,:)=PARCOMP(IDrows(f)).beta(cc,1:(NC/2));
    end
end
for cc=1:length(cov.labels)
    if ~(cc==cov.ConstantID) %except the constant
        nexttile;
        boxchart(figbetas.(cov.labels{cc}),'markerstyle','.');
        ylabel('Beta values'); xlabel('Channels'); title(cov.labels{cc});yline(0,'--');
        nexttile;
        boxchart(figbetas.(cov.labels{cc})','markerstyle','.');
        ylabel('Beta values'); xlabel('Blocks'); title(cov.labels{cc});yline(0,'--');
        xticklabels([PARCOMP(IDrows).file]);
    end
end
title(figg,maintitle,'fontweight','bold')
subtitle(figg,smalltitle)
saveas(figg,fullfile(nirsPATH,'BetaPhysio_distributionplot.fig'));
saveas(figg,fullfile(nirsPATH,'BetaPhysio_distributionplot.png'));
close
clear figg figbetas

%% figure of corrected data
if VIEWplot
    figg2=figure('units','normalized','outerposition',[0 0 1 1]);
    figg2=tiledlayout(5,3,'TileSpacing','Compact','Padding','Compact');
    sizefig=ceil(length(IDrows)/5);
    yy=[];
    for f=1:length(IDrows)
        if any(f==6:5:70)
            %adjust the Y limits so that all are the same!
            for ir=1:15
                nexttile(ir);
                xlim([0 size(PARCOMP(IDrows(1)).data,1)])
                ylim([min(yy) max(yy)])
            end
            %save fig
            saveas(figg2,fullfile(nirsPATH,['PhysioCorr_b' num2str(f-5) '-' num2str(f-1) '.png']))
            close; clear figg2
            figg2=figure('units','normalized','outerposition',[0 0 1 1]);
            figg2=tiledlayout(5,3,'TileSpacing','Compact','Padding','Compact');
            yy=[];
        end

        nexttile;
        plot(PARCOMP(IDrows(f)).data(:,1:(NC/2))-mean(PARCOMP(IDrows(f)).data(1:38,1:(NC/2)),'omitnan')); yy=[yy ylim];
        ylabel(['Block ' num2str(PARCOMP(IDrows(f)).file)],'fontweight','bold','FontSize',12);
        if any(f==1:5:70); title('HBO initial data'); end

        nexttile;
        plot(PARCOMP(IDrows(f)).Xm(:,1:(NC/2)));  yy=[yy ylim];   hold on
        plot(mean(PARCOMP(IDrows(f)).Xm(:,1:(NC/2)),2,'omitnan'),'Color','k','LineWidth',2);
        yy=[yy ylim];
        if any(f==1:5:70); title('Physio component'); end

        nexttile;
        plot(PARCOMP(IDrows(f)).dataCORR(:,1:(NC/2))-mean(PARCOMP(IDrows(f)).dataCORR(1:38,1:(NC/2)),'omitnan'));  yy=[yy ylim];
        if any(f==1:5:70); title('Corrected data'); end
    end
    %adjust the Y limits so that all are the same!
    for ir=1:(3*(f-5*floor(f/5)))
        nexttile(ir);
        xlim([0 size(PARCOMP(IDrows(1)).data,1)])
        ylim([min(yy) max(yy)])
    end
    %save fig
    saveas(figg2,fullfile(nirsPATH,['PhysioCorr_b' num2str(5*floor(f/5)+1) '-' num2str(f) '.png']))
    close; clear figg2
end

%% save nirs mat and PARCOMP

fprintf('Update NIRSmat COMPLETED ...%s\n*\n**\n***\n',job.NIRSmat{1})
save(job.NIRSmat{1},'NIRS');
save(fullfile(nirsPATH,'SelectedFactors.mat'),'PARCOMP');


out.NIRSmat = job.NIRSmat;
end
