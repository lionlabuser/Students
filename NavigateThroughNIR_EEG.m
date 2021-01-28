%infos and code taken from CurrentSetPath2018\spm12\toolbox\LIONirs\nirs_run_E_extractcomponent.m
%from lines 893 (GLM part of the function)

NIRSdirectory='C:\Data\data_NIRS\ELAN\ANALYSED\Martine_0m\BB048\Segment'; 
            %CHANGE HERE according to your NIRS.mat directory

%download NIRS.mat structure            
load([ NIRSdirectory '\NIRS.mat'])

%number of channels (hbo AND hbr)
NC=NIRS.Cf.H.C.N;

%sampling rate : Frequence echantillonage Hz
fsNIRS=NIRS.Cf.dev.fs;

%% upload .nir FILES
% FUNCTION fopen_NIR = to open .nir files
%  Convention: .nir files written in multiplexed format (channels x time),
%  as a long vector, and require specification of number of channels
%  while .nirs files are written in matrix format (time x channels)
%  output of fopen_NIR will be in multiplexed format

%.nir files: name and locations.
nirfiles=NIRS.Dt.fir.pp(end).p; %end = from the latest step of NIRS preprocessing

d1 = fopen_NIR(nirfiles{1},NC); %load e.g. the first .nir data

%time axis from 0 to end based on sampling rate and number of points from
%the .nir file
timeSECaxe = (1/fsNIRS):(1/fsNIRS):(size(d1,1)*1/fsNIRS);

%% process to upload .dat file (EEG)
% 1. Make sure synchronisation is complete (check if you have a
% <sync_timesec> variable into fields of NIRS.Dt.AUX
% 2. AUX1files & AUX1time_start are related = therefore the first data of
% each need to go together
% 3. Also need to have a end time (because we'll want to import only a
% segment of the EEG file <AUX1time_stop>
% 4. Use f_openEEG to open .dat
%       FUNCTION fopen(filename) >> open whole bloc
%       fopen(filename tstart, tstop) >> open a segment only
%       Read EEG.dat format analyser generic data export
%       Input file > EEG generic data export .dat, .vmrk, .vhdr
%       Output > data will be in vectorized (time x channel matrix) event if they
%                are save in multiplex you could use the output of this function directly 
%                in vectorised data matrix for the whole file
%              > info (vhdr information)
% 5. Need sampling rate of EEG data (500Hz) + ratio related to NIRSfs (64)
% 6. Use downsample with the data and the ratio >> create a new dataset
% (AUXdatanew>

%labels of AUX attached to the NIRS.mat
AUXlabel={NIRS.Dt.AUX.label};

AUX1files=NIRS.Dt.AUX(1).pp(end).p; %AUX files (or file) related to the .nir file
%(if multiple blocks in .nir: it is sometimes the same repeated several times)

numberAUX1files=length(AUX1files); %here you have the total number of files 
%(should be equal to your number of blocks)

AUX1time_start=cell2mat(NIRS.Dt.AUX(1).pp(end).sync_timesec); %the syncrhonisation time
%in seconds from the AUX file (eg. eeg) that has been previously
%synchronized from the segmentation step.

AUX1time_stop=AUX1time_start+timeSECaxe(end); %to calculate the <end time>
%for the EEg data to be downloaded (based on the length of the bock

%FUNCTION fopen_EEG here!!
[AUXdata,AUXinfo,~,~] = fopen_EEG(AUX1files{1}, AUX1time_start(1), AUX1time_stop(1)); 
        %1 = represents the first "block" of the NIRS data that has a synchronized EEg data

AUXinfo.name_ele %name of each channels in the AUX file

fsAUX =1/(AUXinfoBV.SamplingInterval/1000000); %Frequence echantillonage Hz EEG

[~,ratio] = rat(fsNIRS/fsAUX,0.0001); %ratio to know how to downsample
%the EEG channels so that it fits the length of NIRS data

AUXdatanew=downsample(AUXdata,ratio);
 %one line for EACH aux CHANNEL (ex. saturometre \ respiration)

 %check if same size as NIRS data
 if size(AUXdatanew,1)==size(d1,1)
     fprintf('NIRS and EEG data have the same length: %d data points\n',size(AUXdatanew,1))
 else
     fprintf('NIRS data have %d data points\nEEG data DOESN''T have the same length: %d data points\n ',...
         size(d1,1),size(AUXdatanew,1))
 end
 
 
%normalize amplitude of EEG data for each EEG channel 
%(into Zscore: minus mean, divided %by STD)
for c=1:size(AUXdatanew,2) %for each channel
    m=mean(AUXdatanew(:,c),'omitnan');
    st=std(AUXdatanew(:,c),0,'omitnan');
    AUXdatanew(:,c)=(AUXdatanew(:,c)-m)/st; %update data
    
    
   REGRESSORS.(AUXinfo.name_ele{c})=AUXdatanew(:,c); %save each AUX 
        % channel independently in a REGRESSORS structure
end
    

%% OBTAIN SPATIAL MEAN (average of all channels, one for HbO, one for HbR)

HBOchan=1:NC/2; %channel numbers for hbo
HBRchan=1+NC/2 : NC; %channel numbers for hbr

%save each GlobalAVG in a REGRESSORS structure
REGRESSORS.globalAVG_HbO=mean(d1(:,HBOchan),2,'omitnan');
REGRESSORS.globalAVG_HbR=mean(d1(:,HBRchan),2,'omitnan');