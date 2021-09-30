function  outtxt=DataQualityCheck_CardiacCrossCorrelation(participant,folder,targetNIRS)

%Method by
% Luca Pollonini, Cristen Olds, Homer Abaya, Heather Bortfeld,
% Michael S. Beauchamp, John S. Oghalai (2014).
% Auditory cortex activation to natural speech and simulated cochlear
% implant speech measured with functional near-infrared spectroscopy,
% Hearing Research, Volume 309, Pages 84-93
% https://doi.org/10.1016/j.heares.2013.11.007

%PRE-REQUIS =
% 1) create a subfolder for Cardiac Cross correlation
%
% 2) filtering the data to get cardiac peak (e.g. infant or
%       newborn = good to filter between 1 and 3 Hx, as the peak is around 2Hz

%OUTPUTS =
% 1) save the correlation coefficient matrice within the specified <folder>
% 2) save a text file with a summary of bach channels with their respective bad blocks
% 3) update NIRS.mat file contained in the <targetNIRS> folder

%% Values to adapt
%participant='BB016';
%folder='C:\Users\laura\AnalysesNIRS\MART0m\BB016\Segment\NormV2\CardiacCrossCorr\';
%targetNIRS= 'C:\Users\laura\AnalysesNIRS\MART0m\BB016\Segment\NormV2\'; %target NIRS.mat folder to automatically mark the bad channels
NC=108; %number of channels
threshold=.75; %correlation threshold (here = same as Pollonini 2014)
CHANGE=1; %overwrite nirs.mat file for bad channels

%% Beginning of automatic script
if ~strcmp(folder(end), filesep)
    folder=[folder filesep];
end
if exist('targetNIRS','var')
    if~strcmp(targetNIRS(end), filesep)
        targetNIRS=[targetNIRS filesep];
    end
end

nirsfile=dir([folder 'f*.nir']); %getting the names of .nir files from FILTERED preprocessing step

for f=1:length(nirsfile)
    nir=fullfile(nirsfile(f).folder,nirsfile(f).name);
    
    d = fopen_NIR(nir,NC);
    
    for ch=1:(NC/2)
        if sum(isnan(d([ch ch+NC/2],:)),'all')>0 %at least one of the channel is NAN
            R(ch,f)=0;
        else
            wv1=(d(ch,:)-mean(d(ch,:)))/std(d(ch,:));
            wv2=(d(ch+NC/2,:)-mean(d(ch+NC/2,:)))/std(d(ch+NC/2,:));
            r=corrcoef(wv1,wv2);
            R(ch,f)=r(1,2);
            clear r wv*
        end
        
    end
    clear d
end

%save R matrice
R2=reshape(R<threshold,size(R));
save([folder participant '_WVcorrcoef.mat'],'R','R2')


%% change to bad channels in NIRS.mat file
%if targetfolder is not specified (if the variable does not exist)
% = automatically uses the folder in which the Cardiac subfolder is in
if CHANGE ==1
    if ~exist('targetNIRS','var')
        tmp=split(folder,filesep);
        targetNIRS=[fullfile(tmp{1:end-2}) filesep];
        clear tmp
    end
    load([targetNIRS 'NIRS.mat'],'NIRS');
    if size([R2;R2]) == size(NIRS.Cf.H.C.ok)
        NIRS.Cf.H.C.ok([R2;R2]) =0 ;
    else
        error('Size of coefficient matrice does not match size of badchannel identification in NIRS.mat file')
    end
    save([targetNIRS 'NIRS.mat'],'NIRS');
end

%% create a summary text file

summary=['***********************************\n'...
    'QUALITY CHECK based on zero-lag cardiac correlation between wavelengths\n'...
    '***********************************\n' ...
    'Date: ' date ...
    '\nParticipant ID: ' participant ...
    '\nFolder ID: '  replace(folder,filesep,[filesep filesep])...
    '\nNIRS.mat file in which those channels are identified as bad: ' replace(targetNIRS,filesep,[filesep filesep]) 'NIRS.mat'...
    '\nThreshold for bad correlation: r=' num2str(threshold)...
    '\n_____________________\nChannel with bad zero-lag correlation: \n'];

badch=find(sum(R2,2));
if isempty(badch)
    summary=[summary 'NO BAD CHANNELS'];
else
    for b=1:length(badch)
        badblock=find(R2(badch(b),:));
        summary=[summary ...
            'CH #' num2str(badch(b)) ' - '...
            num2str(length(badblock)) 'bad blocks out of ' num2str(length(nirsfile)) ...
            '\n (blocks #' num2str(badblock) ')\n\n'];
    end
end

%save into .txt file
outtxt=[folder participant '_CardiacQualityCheck.txt'];
fileID=fopen(outtxt,'w');
fprintf(fileID,summary);
fclose(fileID);

end
