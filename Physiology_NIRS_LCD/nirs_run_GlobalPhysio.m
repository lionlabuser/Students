function out = nirs_run_GlobalPhysio(job)
% job:
%   - job.NIRSmat = {'path +  NIRS.mat'}
%   - job.physzone = {'path + filename + .zone'} %%NOTE REGARDING THE ZONE
%   FILE: IT CAN CONTAIN EITHER ALL HBO AND HBR CHANNELS OR ONLY HBO
%   CHANNELS. BOTH WORK IN THE CURRENT SCRIPT. BUT IN JULIE'S SCRIPTS, CAREFUL TO ONLY PUT HBO
%   CHANNELS IN THE ZONE FILES
%   - job.multimodalPATH={'path that ends with global'}
%   - job.trig = [2 3 4] trigger or triggers related to the task. if
%       resting state, then enter 0.
%   - job.globalavg = [1 or 0] If we want to extract the global average.
%   - job.globalpca = [1 or 2 or 3 or 0] Normal PCA
%       Refers to the number of components that will be extracted 
%       (max the number of channels. BUT CAREFUL IT CREATED A VERY LARGE 
%       SET OF DATA... not really useful!) If no PCA, then enter 0.
%   - job.spatialpca = [1 or 0] If we want to extract the global component
%       from spatial smoothing PCA
%
%   The current function doesnt take for account of the yellow
%   identification INE SEPARATING COMPONENTS. If yellow identification accounts for block separation
%   (e.g. for resting-state data), you should first segment your data
%   in blocks or copy-paste this script and adapt it (can be based on
%   nirs_run_E_extractcomponent.m).
%
%   The current script is inspired from nirs_run_E_extractcomponent.m -
%   physiology section. Extract global average / global pca from the data
%   and create .dat file for each block (like an AUXILIARY). Update NIRS.mat information
%   accordingly regarding AUX files (link those as auxiliaries, + make sure it is synchronized with the block). 
% Created 2021-02-23 by LCD.
% 
% UPDATE 2021-03-04
% 1) Add spatial filtered PCA, based on the article from Noah et al. 2021 
%     (https://doi.org/10.1117/1.NPh.8.1.015004)
% 	& Zhang et al. 2016 (https://doi.org/10.1117/1.NPh.3.1.015004)
%     The function is at the end of the current script, and each step 
%     is detailed in commentary.
% 
% 2) Write/Update SelectedFactors.mat (inspired from the GLM part of 
%     the nirs_run_E_extractcomponent.m script). Important to know what 
%     each column means in this file:
%     - data: original data (time x channels)
%     - Xm: the global component (time x channels): each channel has 
%         its "own" global component. It represents the singular shape 
%         of the component multiplied by the beta ot the channel
%     - dataCORR: data (minus) Xm (substraction of the global component)
%         THEREFORE IT DOES THE SUBSTRACTION! BUT YOU WILL FIND THE DATA IN
%         THE SELECTEDFACTORS MATFILE. 
%     - Each block x type of global decomposition = represents a line 

for filenb=1:size(job.NIRSmat,1) %Loop over all subjects
    
    load(job.NIRSmat{filenb});
    
    lst = length(NIRS.Dt.fir.pp); %last step of nirs processing
    rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
    NC = NIRS.Cf.H.C.N; %number of channels
    fs = NIRS.Cf.dev.fs; %NIRS sampling rate
    wl=NIRS.Cf.H.C.wl';
    %nbchminimum = 5/100; %0.05;             %en pourcentage
    if isfield(job,'nbtimeminimum')
        nbtimeminimum=job.nbtimeminimum;
    else
    nbtimeminimum=10/40; %the % time non artifacted to be considered as a good channel
    end
    load(job.physzone{filenb}   ,'-mat');
    
    if ~contains(job.multimodalPATH{filenb}(end),filesep)
        job.multimodalPATH{filenb}(end+1)=filesep;
    end
    if job.globalavg || job.globalpca
    if ~exist(job.multimodalPATH{filenb},'dir')
        mkdir(job.multimodalPATH{filenb});
    end
    end
    
    %check SelectedFactors file
    [nirsPATH,~,~] = fileparts(job.NIRSmat{filenb});
    if exist(fullfile(nirsPATH,'SelectedFactors.mat'),'file')
        load(fullfile(nirsPATH,'SelectedFactors.mat'))
        Prow=length(PARCOMP)+1;
    else
        PARCOMP=struct;
        Prow=1;
    end
    
    %Update NIRS.mat file for AUX
    if job.globalavg
        AUXavg=numel(NIRS.Dt.AUX)+1;
        NIRS.Dt.AUX(AUXavg).label='GlobalAVG';
    end
    if job.globalpca
        AUXpca=numel(NIRS.Dt.AUX)+1;
        NIRS.Dt.AUX(AUXpca).label='GlobalPCA';
    end
    
    %%Determine regressor zone and channel zones
    %(first column regresor channel to average, second row column ch to apply)
    RegZone=[];
    ChanZone=[];
    for iz =1:numel(zone.label)
        tmp = upper(zone.label{iz});
        if numel(tmp)>9
            if strcmp(tmp(1:9),'REGRESSOR')
                RegZone = [RegZone,iz];
            end
        end
    end
    if numel(RegZone)>1
        error('More than 1 regressor zone. Script made for global Physio.\n')
    end
%     for iR = 1:numel(RegZone)
%         tmp = upper(zone.label{RegZone(iR)});
%         zoneidentification = tmp(10:end);
%         for iz = 1:numel(zone.label)
%             tmpzone = upper(zone.label{iz});
%             if strcmp( strtrim(zoneidentification), strtrim(tmpzone))
%                 ChanZone = [ChanZone,iz];
%             end
%         end
%     end
%     
%     if numel(ChanZone) ~= numel(RegZone)
%         msgbox('Regressor Zone and Zone must be in equal number in the zone list')
%         return
%     end
    
    for f = 1:size(rDtp,1)
        d = fopen_NIR(rDtp{f,1},NC)';
        [cPATH,cFILE,~] = fileparts(rDtp{f});
        tstart=1;
        tstop=size(d,1);
        badch = NIRS.Cf.H.C.ok(:,f);
        if badch==0 %if all channels are bad, switch to next block!
            continue
        end
        
        %%check yellow
        mrks = [];
        ind = [];
        noise =  logical(zeros(size(d)));
        badchY=ones([size(d,2) 1]);
        [ind_dur_ch] = read_vmrk_find(fullfile(cPATH,[cFILE '.vmrk']),{'bad_step'});
        if ~isempty(ind_dur_ch)
            badind = find((ind_dur_ch(:,1)+ind_dur_ch(:,2))>size(noise,1));
            if ~isempty(badind)
                disp(['Warning file ' fullfile(cPATH,[cFILE '.vmrk']) ' marker : ' num2str(badind') ' are out of range in the data file'])
                ind_dur_ch(badind,2)=size(noise,2)- ind_dur_ch(badind,1);
            end
            for Chan = 1:size(noise,2)
                mrks = find(ind_dur_ch(:,3)==Chan);
                ind = ind_dur_ch(mrks,1);
                indf = ind + ind_dur_ch(mrks,2) - 1;
                if ~isempty(ind)
                    for i = 1:numel(ind)
                        noise(ind(i):indf(i),Chan) = 1;
                    end
                end
                if sum(ind_dur_ch([mrks],2))/size(noise,1) > nbtimeminimum
                    badchY(Chan)=0;
                end
            end
        end
        
        %% DATA PREPARATION
        chlistRegressor = zone.plotLst{RegZone};
        
        %FIRST WAVELENGTH ~ HBO
        %Create a Channel regressor list with only the HbO channels!
        for tt=1:length(chlistRegressor)
            if wl(chlistRegressor(tt))==1 %check up if only HbO
                goodhbo(tt)=tt;
            end
        end
        CHhbo =chlistRegressor(goodhbo);
        
        idbad = find(badch( CHhbo )==0|badchY( CHhbo )==0); %remove exclude channel from regressor
        if ~isempty(idbad)
            CHhbo (idbad) = [];
            if isempty(CHhbo)
                disp(['No good channel in the regressor zone please verify your zone ', zone.label{RegZone}])
                continue
            end
        end
        data = d(tstart:tstop,CHhbo);
        
        %SECOND WAVELENGTH ~ HBR
        CHhbr =  CHhbo + size(d,2)/2;
        data2 = d(tstart:tstop,CHhbr);
        
        %BASELINE CORRECTION
        if ~job.trig %IF RESTING: baseline correction using global mean across time
            data=data-mean(data,'omitnan');
            data2=data2-mean(data2,'omitnan');
        else %IF TASK-BASED: baseline correction using mean of prestim time
            for tt=1:length(job.trig)
              
                [t1,t2]=find(NIRS.Dt.fir.aux5{f}==job.trig(tt));
                if t1
                    trigpos=NIRS.Dt.fir.aux5{f}(t1,t2+1);
                    break
                end
            end
            data=data-mean(data(tstart:trigpos,:),'omitnan');
            data2=data2-mean(data2(tstart:trigpos,:),'omitnan');
        end
        
        %% GLOBAL AVERAGE
        if job.globalavg
            AUX.data(:,1) = mean(data,2,'omitnan');
            AUX.data(:,2) = mean(data2,2,'omitnan');
            
            %% WRITE PARCOMP and BETAs
            tmpbeta(1:NC)=nan;
            tmpXm(tstart:tstop,1:NC)=nan;
            tmpXcorr(tstart:tstop,1:NC)=nan;
            for CH = 1:numel(CHhbo)
                tmpbeta(CHhbo(CH)) = dot(AUX.data(:,1),data(:,CH))/dot(AUX.data(:,1),AUX.data(:,1));
                tmpXm(:,CHhbo(CH)) = tmpbeta(CHhbo(CH)).* AUX.data(:,1);
                tmpXcorr(:,CHhbo(CH)) = data(:,CH) - tmpXm(:,CHhbo(CH));
            end
            for CH = 1:numel(CHhbr)
                tmpbeta(CHhbr(CH)) = dot(AUX.data(:,2),data2(:,CH))/dot(AUX.data(:,2),AUX.data(:,2));
                tmpXm(:,CHhbr(CH)) = tmpbeta(CHhbr(CH)).* AUX.data(:,2);
                tmpXcorr(:,CHhbr(CH)) = data2(:,CH) - tmpXm(:,CHhbr(CH));
            end
            
            tmpdata(tstart:tstop,1:NC)=nan;
            tmpdata(:,CHhbo)=data;
            tmpdata(:,CHhbr)=data2;
            
            
            PARCOMP(Prow).file= f;
            PARCOMP(Prow).filestr =  sprintf('Bl%02.0f', f);
            PARCOMP(Prow).label= ['GlobalAVG_' PARCOMP(Prow).filestr];
            PARCOMP(Prow).type = 'GLM';
            PARCOMP(Prow).beta = tmpbeta;
            PARCOMP(Prow).std= 0;
            PARCOMP(Prow).AUX.data= AUX.data;
            PARCOMP(Prow).AUX.label= {'GlobalAVG HbO' 'GlobalAVG HbR' };
            PARCOMP(Prow).data = tmpdata ;
            PARCOMP(Prow).Xm = tmpXm;
            PARCOMP(Prow).dataCORR = tmpXcorr;
            PARCOMP(Prow).indt = [tstart :tstop]; %indice de temps.
            PARCOMP(Prow).listgood =  [CHhbo;CHhbr];
            PARCOMP(Prow).module  = lst;
            PARCOMP(Prow).modulestr = NIRS.Dt.fir.pp(lst).pre;
            PARCOMP(Prow).ComponentToKeep = 1;
            PARCOMP(Prow).idreg = 1;
            PARCOMP(Prow).topo =  tmpbeta(PARCOMP(Prow).ComponentToKeep,:);
            Prow=Prow+1;
            clear tmp*
            
            %% WRITE DAT FILE
            AUX.infoBV.dataformat='BINARY';
            AUX.infoBV.dataorientation='VECTORIZED';
            AUX.infoBV.DataType='TIMEDOMAIN';
            AUX.infoBV.NumberOfChannels=size(AUX.data,2);
            AUX.infoBV.DataPoints=size(AUX.data,1);
            AUX.infoBV.SamplingInterval=(1/fs)*1000000;
            AUX.infoBV.BinaryFormat='IEEE_FLOAT_32';
            AUX.infoBV.name_ele{1}='GlobalAVG hbo';
            AUX.infoBV.name_ele{2}='GlobalAVG hbr';
            
            AUX.ind_dur_ch=[ NIRS.Dt.fir.aux5{f}(:,2) ones([size(NIRS.Dt.fir.aux5{f},1) 1]) zeros([size(NIRS.Dt.fir.aux5{f},1) 1])];
            for tn = 1:size(NIRS.Dt.fir.aux5{f},1)
                AUX.marker{tn,1}='trigger';
                AUX.marker{tn,2}=['S ' num2str(NIRS.Dt.fir.aux5{f}(tn,1))];
            end
            
            fileoutAUX=[job.multimodalPATH{filenb} 'AVG_b' num2str(f) '.dat'];
            fwrite_EEG(fileoutAUX,AUX,1,AUX.infoBV.DataPoints );
            disp(fileoutAUX)
            NIRS.Dt.AUX(AUXavg).pp.p{f,1}=fileoutAUX;
            NIRS.Dt.AUX(AUXavg).pp.sync_timesec{f,1}=0;
            clear AUX fileoutAUX
        end
        
        %% GLOBAL PCA
        if job.globalpca
            %u : temporal signature of the component
            %s : relative weigth of the component
            %v : spatial signature of the component
            dtemp = data'*data;
            [v,s,~]=svd(dtemp);
            u = data*v*inv(s);
            
            dtemp2 = data2'*data2; %2nd wavelength
            [v2,s2,~]=svd(dtemp2);
            u2 = data2*v2*inv(s2);
            
            for c=1:job.globalpca %component number
                compPCA=u(:,c)*s(c,c)*v(:,c)';
                compPCA2=u2(:,c)*s2(c,c)*v2(:,c)';
                AUX.data(:,c) = mean(compPCA,2,'omitnan');
                AUX.data(:,c+job.globalpca) = mean(compPCA2,2,'omitnan'); %2nd wavelength
                AUX.infoBV.name_ele{c}=['GlobalPCA comp' num2str(c) ' hbo'];
                AUX.infoBV.name_ele{c+job.globalpca}=['GlobalPCA comp' num2str(c) ' hbr']; %2nd wavelength
                
                %% WRITE PARCOMP and BETAs
                tmpbeta(1:NC)=nan;
                tmpXm(tstart:tstop,1:NC)=nan;
                tmpXcorr(tstart:tstop,1:NC)=nan;
                for CH = 1:numel(CHhbo)
                    tmpbeta(CHhbo(CH)) = dot(AUX.data(:,c),data(:,CH))/dot(AUX.data(:,c),AUX.data(:,c));
                    tmpXm(:,CHhbo(CH)) = tmpbeta(CHhbo(CH)).* AUX.data(:,c);
                    tmpXcorr(:,CHhbo(CH)) = data(:,CH) - tmpXm(:,CHhbo(CH));
                end
                for CH = 1:numel(CHhbr)
                    tmpbeta(CHhbr(CH)) = dot(AUX.data(:,c+job.globalpca),data2(:,CH))/dot(AUX.data(:,c+job.globalpca),AUX.data(:,c+job.globalpca));
                    tmpXm(:,CHhbr(CH)) = tmpbeta(CHhbr(CH)).* AUX.data(:,c+job.globalpca);
                    tmpXcorr(:,CHhbr(CH)) = data2(:,CH) - tmpXm(:,CHhbr(CH));
                end
                
                tmpdata(tstart:tstop,1:NC)=nan;
                tmpdata(:,CHhbo)=data;
                tmpdata(:,CHhbr)=data2;
                
                tmpPCA(tstart:tstop,1:NC)=nan;
                tmpPCA(:,CHhbo)=compPCA;
                tmpPCA(:,CHhbr)=compPCA2;
                
                PARCOMP(Prow).file= f;
                PARCOMP(Prow).filestr =  sprintf('Bl%02.0f', f);
                PARCOMP(Prow).label= ['GlobalPCA_C' num2str(c) '_' PARCOMP(Prow).filestr];
                PARCOMP(Prow).type = 'GLM';
                PARCOMP(Prow).beta = tmpbeta;
                PARCOMP(Prow).std= 0;
                PARCOMP(Prow).AUX.data= AUX.data(:,[c c+job.globalpca]);
                PARCOMP(Prow).AUX.Xm=tmpPCA;
                PARCOMP(Prow).AUX.dataCORR=tmpdata - tmpPCA;
                PARCOMP(Prow).AUX.label= {'GlobalPCA HbO' 'GlobalPCA HbR' };
                PARCOMP(Prow).data = tmpdata ;
                PARCOMP(Prow).Xm = tmpXm;
                PARCOMP(Prow).dataCORR = tmpXcorr;
                PARCOMP(Prow).indt = [tstart :tstop]; %indice de temps.
                PARCOMP(Prow).listgood =  [CHhbo;CHhbr];
                PARCOMP(Prow).module  = lst;
                PARCOMP(Prow).modulestr = NIRS.Dt.fir.pp(lst).pre;
                PARCOMP(Prow).ComponentToKeep = 1;
                PARCOMP(Prow).idreg = 1;
                PARCOMP(Prow).topo =  tmpbeta(PARCOMP(Prow).ComponentToKeep,:);
                Prow=Prow+1;
                clear tmp*

            end
            
            AUX.infoBV.dataformat='BINARY';
            AUX.infoBV.dataorientation='VECTORIZED';
            AUX.infoBV.DataType='TIMEDOMAIN';
            AUX.infoBV.NumberOfChannels=size(AUX.data,2);
            AUX.infoBV.DataPoints=size(AUX.data,1);
            AUX.infoBV.SamplingInterval=(1/fs)*1000000;
            AUX.infoBV.BinaryFormat='IEEE_FLOAT_32';
            
            AUX.ind_dur_ch=[ NIRS.Dt.fir.aux5{f}(:,2) ones([size(NIRS.Dt.fir.aux5{f},1) 1]) zeros([size(NIRS.Dt.fir.aux5{f},1) 1])];
            for tn = 1:size(NIRS.Dt.fir.aux5{f},1)
                AUX.marker{tn,1}='trigger';
                AUX.marker{tn,2}=['S ' num2str(NIRS.Dt.fir.aux5{f}(tn,1))];
            end
            
            fileoutAUX=[job.multimodalPATH{filenb} 'PCA' num2str(job.globalpca) 'comp_b' num2str(f) '.dat'];
            fwrite_EEG(fileoutAUX,AUX,1,AUX.infoBV.DataPoints );
            disp(fileoutAUX)
            NIRS.Dt.AUX(AUXpca).pp.p{f,1}=fileoutAUX;
            NIRS.Dt.AUX(AUXpca).pp.sync_timesec{f,1}=0;
            clear AUX fileoutAUX
        end
        
        
        %% GLOBAL spatial filtered PCA
        if job.spatialpca
            
           
                
            [filteredD,globalComp,individualComp]=spatialPCA(data,CHhbo,zone,0,0.8);
            [filteredD2,globalComp2,individualComp2]=spatialPCA(data2,CHhbo,zone,0,0.8); %as hbo and 
                %hbr have the same spatial position, easier to take the good hbo list instead of the hbr list
            
              AUX.data(:,1) = mean(globalComp,2,'omitnan');
                AUX.data(:,2) = mean(globalComp2,2,'omitnan'); %2nd wavelength
                
            %% WRITE PARCOMP and BETAs
                tmpbeta(1:NC)=nan;
                tmpXm(tstart:tstop,1:NC)=nan;
                tmpXcorr(tstart:tstop,1:NC)=nan;
                for CH = 1:numel(CHhbo)
                    tmpbeta(CHhbo(CH)) = dot(AUX.data(:,1),data(:,CH))/dot(AUX.data(:,1),AUX.data(:,1));
                    tmpXm(:,CHhbo(CH)) =globalComp(:,CH);
                    tmpXcorr(:,CHhbo(CH)) =filteredD(:,CH);
                    %tmpXm(:,CHhbo(CH)) = tmpbeta(CHhbo(CH)).* AUX.data(:,c);
                    %tmpXcorr(:,CHhbo(CH)) = data(:,CH) - tmpXm(:,CHhbo(CH));
                end
                for CH = 1:numel(CHhbr)
                    tmpbeta(CHhbr(CH)) = dot(AUX.data(:,2),data2(:,CH))/dot(AUX.data(:,2),AUX.data(:,2));
                    tmpXm(:,CHhbr(CH)) =globalComp2(:,CH);
                    tmpXcorr(:,CHhbr(CH)) =filteredD2(:,CH);
                    %tmpXm(:,CHhbr(CH)) = tmpbeta(CHhbr(CH)).* AUX.data(:,c+job.globalpca);
                    %tmpXcorr(:,CHhbr(CH)) = data2(:,CH) - tmpXm(:,CHhbr(CH));
                end
                
                tmpdata(tstart:tstop,1:NC)=nan;
                tmpdata(:,CHhbo)=data;
                tmpdata(:,CHhbr)=data2;
                
                PARCOMP(Prow).file= f;
                PARCOMP(Prow).filestr =  sprintf('Bl%02.0f', f);
                PARCOMP(Prow).label= ['SpatialPCA_' PARCOMP(Prow).filestr];
                PARCOMP(Prow).type = 'GLM';
                PARCOMP(Prow).beta = tmpbeta;
                PARCOMP(Prow).std= 0;
                PARCOMP(Prow).AUX.data= AUX.data;
                PARCOMP(Prow).AUX.IndividualComp={individualComp, individualComp2};
                PARCOMP(Prow).AUX.label= {'SpatialPCA HbO' 'SpatialPCA HbR' };
                PARCOMP(Prow).data = tmpdata ;
                PARCOMP(Prow).Xm = tmpXm;
                PARCOMP(Prow).dataCORR = tmpXcorr;
                PARCOMP(Prow).indt = [tstart :tstop]; %indice de temps.
                PARCOMP(Prow).listgood =  [CHhbo;CHhbr];
                PARCOMP(Prow).module  = lst;
                PARCOMP(Prow).modulestr = NIRS.Dt.fir.pp(lst).pre;
                PARCOMP(Prow).ComponentToKeep = 1;
                PARCOMP(Prow).idreg = 1;
                PARCOMP(Prow).topo =  tmpbeta(PARCOMP(Prow).ComponentToKeep,:);
                Prow=Prow+1;
                clear tmp*
            
        end
        
        disp(['block ' num2str(f) ' done'] );
        clear d badch chlistRegressor good* CH* dat* trigpos 
    end
    
    %save updated nirsmat
    fprintf('Update NIRSmat COMPLETED ...%s\n*\n**\n***\n',job.NIRSmat{filenb})
    save(job.NIRSmat{filenb},'NIRS');
    save(fullfile(nirsPATH,'SelectedFactors.mat'),'PARCOMP');
end
out.NIRSmat = job.NIRSmat;
end


function [filteredD,globalComp,individualComp]=spatialPCA(D,chlist,z,Visualize,kernel)
% SCRIPT based on the article from Noah et al. 2021 (https://doi.org/10.1117/1.NPh.8.1.015004)
% & Zhang et al. 2016 (https://doi.org/10.1117/1.NPh.3.1.015004)
%
if nargin <3
    Visualize=0;
    if nargin <4
        kernel=0.8;
    end
end

%% Calculate channel distance

nchan=length(chlist); %number of channels
Dist=zeros(nchan); 
rayon=zeros(nchan);
for a=1:nchan
    for b=a+1:nchan
        rayon(a,b)=sqrt( z.pos(chlist(a),1)^2 + z.pos(chlist(a),2)^2 + z.pos(chlist(a),3)^2);
        rayon(b,a)=rayon(a,b);
        
        tempdistance=sqrt((z.pos(chlist(a),1)-z.pos(chlist(b),1))^2 + (z.pos(chlist(a),2)-z.pos(chlist(b),2))^2+ (z.pos(chlist(a),3)-z.pos(chlist(b),3))^2);
        Dist(a,b)=2*rayon(a,b)*asin(tempdistance/(2*rayon(a,b))); %distance matrice between channels (around a sphere)
        Dist(b,a)=Dist(a,b);
    end
end
Dist=real(Dist);
rayon(rayon==0)=NaN;
mrayon=mean(rayon,'all','omitnan');  %the mean radius that allowed to measure the arc length distance between channels

%% KERNEL SIZE
% smoothing kernel set to 46 deg (or 0.8 rad) - Zhang et al 2017 (https://doi.org/10.1117/1.NPh.4.4.041409)
% as our distances between channels are all in centimeters, we need to
% convert the kernel (in degrees) into cm (arc length)
% kernel needs to be in radians (between 0 and 2pi)
if kernel < 2*pi %kernel in radiam
    ksigma=kernel*mrayon; %now in cm  
else  %kernel in degree
    ksigma=kernel*pi/180*mrayon; %now in cm  
end

%% VERIFY data configuration
sizeD = size(D);
if sizeD(2) == nchan
    %do nothing
elseif sizeD(1) == nchan
    D=permute(D, [2 1]);
else
    error(['The data matrice <D> doesn''t have the same number of channels as the zone file.\n'...
        'Please check your data matrice <D> (size: %d x %d points) and your zone file.\n'],sizeD(1),sizeD(2));
end

%% PCA DECOMPOSITION
%SpatialSig = spatial matrice. <spatial pattern of global components>
%             each column = a component // each row = each channel
%TemporalSig = temporal matrice. <temporal pattern of global components>
%             each column = a component  // each row = a time point
%ComponentWeigth = weigth of each component <diagonal matrice>
%             diagonal coordinates (square, eg. (1,1) (2,2)) = each component
squareD = D'*D;
[SpatialSig,ComponentWeigth,~]=svd(squareD);
TemporalSig = D*SpatialSig*inv(ComponentWeigth);

%% GAUSSIAN SMOOTHING OF THE SPATIAL MATRICE
SmoothSpatialSig=zeros(size(SpatialSig));

for ci=1:size(SpatialSig,1) %for each channel (row)
        for cj=1:size(SpatialSig,1)
            wij(ci,cj)=exp((-(Dist(ci,cj))^2)/(2*ksigma^2)); %calculate the multiplication factor based on the distance between channels ci and cj
            %lorsque ci==cj >> wij = 1
        end
        %wij(ci,:)=wij(ci,:)./sum(wij(ci,:)); %TEST1 
        %   NORMALISATION normalized the sum of multiplication factors to
        %   1. FOR EACH CHANNEL, the sum of its multiplication factors
        %   equals 1. For channels that have a lot of neighbors, it makes
        %   that the weigth of the current channel is smaller compared to
        %   the total weigth. For channels that have fewer neighbors,
        %   their own weigth is larger in proportion to the total.

end
wij=(wij./sum(wij,'all')).*size(SpatialSig,1); %NORMALISATION CHOSEN. 
% the sum of each column is approx equal to the sum of each column from the
% initial SpatialSig matrice.
% It takes the WHOLE MATRICE to do the normalization. therefore the sum of
% weigths for each channel is not necessarily 1 (the mean sum across channels = 1)

for vv=1:size(SpatialSig,2) %for each component (column)
        for ci=1:size(SpatialSig,1) %for each channel (row)
            for cj=1:size(SpatialSig,1)
            SmoothSpatialSig(ci,vv)= SmoothSpatialSig(ci,vv) + wij(ci,cj)*SpatialSig(cj,vv);
            if Visualize
                if vv==1 && ci==16 %help to visualize
                    fprintf('Distance: %0.2f -- kernel: %0.10f\n',Dist(ci,cj),wij(ci,cj));
                end
                
            end
        end
        end
end
globalComp=TemporalSig*ComponentWeigth*SmoothSpatialSig';
filteredD=D-globalComp;

for vv=1:size(SpatialSig,2) 
individualComp(:,:,vv)=TemporalSig(:,vv)*ComponentWeigth(vv,vv)*SmoothSpatialSig(:,vv)';
end
if Visualize
figure
subplot(3,1,1)
plot(D)
title('Original data')
subplot(3,1,2)
plot(globalComp)
title('Global component')
subplot(3,1,3)
plot(D-globalComp)
title('Derived neuronal component')
end
end