function out = nirs_run_GlobalPhysio(job)
% job:
%   - job.NIRSmat = {'path +  NIRS.mat'}
%   - job.physzone = {'path + filename + .zone'}
%   - job.multimodalPATH={'path that ends with global'}
%   - job.trig = [2 3 4] trigger or triggers related to the task. if
%       resting state, then enter 0.
%   - job.globalavg = [1 or 0] If we want to extract the global average.
%   - job.globalpca = [1 or 2 or 3 or 4 or 5 or 0] Normal PCA
%       Refers to the number of components that will be extracted. Up
%       to the first 5 most important components. If no PCA, then enter 0.
%
%   The current function doesnt take for account of the yellow
%   identification. If yellow identification accounts for block separation
%   (e.g. for resting-state data), you should first segment your data
%   in blocks or copy-paste this script and adapt it (can be based on
%   nirs_run_E_extractcomponent.m).
%
%
%   The current script is inspired from nirs_run_E_extractcomponent.m -
%   physiology section. Extract global average / global pca from the data
%   and create .dat file for each block. Update NIRS.mat information
%   accordingly regarding AUX files. Created 2021-02-23 by LCD.

for filenb=1:size(job.NIRSmat,1) %Loop over all subjects
    
    load(job.NIRSmat{filenb});
    
    lst = length(NIRS.Dt.fir.pp); %last step of nirs processing
    rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
    NC = NIRS.Cf.H.C.N; %number of channels
    fs = NIRS.Cf.dev.fs; %NIRS sampling rate
    wl=NIRS.Cf.H.C.wl';
    load(job.physzone{filenb}   ,'-mat');
    
    if ~contains(job.multimodalPATH{filenb}(end),'\')
        job.multimodalPATH{filenb}(end+1)='\';
    end
    if ~exist(job.multimodalPATH{filenb},'dir')
        mkdir(job.multimodalPATH{filenb});
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
    
    %     for iR = 1:numel(RegZone)    %%IN COMMENTS BECAUSE WE ONLY WANT TO
    %     EXTRACT THE GLOBAL REGRESSOR, NOT APPLY IT ON THE DATA YET!
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
        [cPATH,cFILE,cEXT] = fileparts(rDtp{f});
        tstart=1;
        tstop=size(d,1);
        badch = NIRS.Cf.H.C.ok(:,f);
        if badch==0 %if all channels are bad, switch to next block!
            continue
        end
        
        for ir = 1:numel(RegZone) %for each regressor
            %% DATA PREPARATION
            chlistRegressor = zone.plotLst{RegZone(ir)};
            
            %FIRST WAVELENGTH ~ HBO
            %Create a Channel regressor list with only the HbO channels!
            for tmpd=1:length(chlistRegressor)
                if wl(chlistRegressor(tmpd))==1 %check up if only HbO
                    goodhbo(tmpd)=tmpd;
                end
            end
            CHhbo =chlistRegressor(goodhbo);
            
            idbad = find(badch( CHhbo )==0); %remove exclude channel from regressor
            if ~isempty(idbad)
                CHhbo (idbad) = [];
                if isempty(CHhbo)
                    disp(['No good channel in the regressor zone please verify your zone ', zone.label{RegZone(ir)}])
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
                    disp('hello')
                    [t1,t2]=find(NIRS.Dt.fir.aux5{f}==job.trig(tt));
                    if t1
                        trigpos=NIRS.Dt.fir.aux5{f}(t1,t2+1);
                        break
                    end
                end
                data=data-mean(data(1:trigpos),'omitnan');
                data2=data2-mean(data2(1:trigpos),'omitnan');
            end
            
            %% GLOBAL AVERAGE
            if job.globalavg
                AUX.data(:,1) = mean(data,2,'omitnan');
                AUX.data(:,2) = mean(data2,2,'omitnan');
                
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
                tmpd = data'*data;
                [v,s,~]=svd(tmpd);
                u = data*v*inv(s);
                
                tmpd2 = data2'*data2; %2nd wavelength
                [v2,s2,~]=svd(tmpd2);
                u2 = data2*v2*inv(s2);
                
                
                for c=1:job.globalpca %component number
                AUX.data(:,c) = u(:,c);
                AUX.data(:,c+job.globalpca) = u2(:,c); %2nd wavelength
                AUX.infoBV.name_ele{c}=['GlobalPCA comp' num2str(c) ' hbo'];
                AUX.infoBV.name_ele{c+job.globalpca}=['GlobalPCA comp' num2str(c) ' hbr']; %2nd wavelength
                    
%                     REG.pca(1).wl='hbo';
%                     REG.pca(1).component(:,c)=u(:,c);
%                     REG.pca(1).spatialW(:,c)=v(:,c);
%                     REG.pca(1).Xm{c}=u(:,c)*s(c,c)*v(:,c)';
%                     
%                     REG.pca(2).wl='hbr';
%                     REG.pca(2).component(:,c)=u2(:,c); 
%                     REG.pca(2).spatialW(:,c)=v2(:,c);
%                     REG.pca(2).Xm{c}=u2(:,c)*s2(c,c)*v2(:,c)';
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
            
            
        end
        
        disp('block done');
        
    end
    
    %save updated nirsmat
    fprintf('Update NIRSmat COMPLETED ...%s\n*\n**\n***\n',job.NIRSmat{filenb})
    save(job.NIRSmat{filenb},'NIRS');
end
out.NIRSmat = job.NIRSmat;
end