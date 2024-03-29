%___________________________________________________________________
% Copyright (C) 2019 LION Lab, Centre de recherche CHU Sainte-Justine
% www.lionlab.umontreal.ca
%___________________________________________________________________
function out = nirs_run_filterAUX(job)
% Filter AUX (EEG) data based on parameters already used for the filtering
% of NIRS data. Create new .dat files and update the NIRS.mat info
% regarding the AUX files.
% INPUT:
%       job: structure with fieldname NIRSmat, containing at least one path+file of a NIRSmat file
%            eg.job.NIRSmat={'C:...\NIRS.mat'};
%          + job.outAUXfolder=string array (the new folder name for
%          filterAUX)
outAUXfolder=job.outAUXfolder;
AUXcolor={[0 0.4470 0.7410],[0.6350 0.0780 0.1840],[0.4660 0.6740 0.1880],[0.4940 0.1840 0.5560]};
for filenb=1:size(job.NIRSmat,1) %Loop over all subjects
    %Load NIRS.mat information
    NIRS = [];
    load(job.NIRSmat{1});
    
    
    liststeps={NIRS.Dt.fir.pp.pre};
    for xx = 1:length(liststeps)
        if contains(liststeps{xx},'Filter')
            fstep=xx;
        end
    end
    if ~fstep
        error('No filter applied to the NIRS data %s/ntherefore no filtering parameters', (job.NIRSmat{1}))
    end
    
    %% SET PARAMETERS
    if job.copynirs || ~isfield(job,'fparameters')
        [lowcut,applylowcut] = str2num(NIRS.Dt.fir.pp(fstep).job.lowcutfreq);
        [highcut ,applyhighcut] = str2num(NIRS.Dt.fir.pp(fstep).job.highcutfreq);
        paddingsym = NIRS.Dt.fir.pp(fstep).job.paddingsymfilter; %symetrie padding on the signal to avoid edge on the filtering
        interpolate = NIRS.Dt.fir.pp(fstep).job.interpolatebadfilter;  %interpolate bad intervals
        filt_ord = NIRS.Dt.fir.pp(fstep).job.filterorder;
        DelPreviousData  = NIRS.Dt.fir.pp(fstep).job.DelPreviousData;
    else
        DelPreviousData  = NIRS.Dt.fir.pp(fstep).job.DelPreviousData;
        if contains(num2str(job.fparameters.lowcut),'same') || contains(num2str(job.fparameters.applylowcut),'same')
            [lowcut,applylowcut] = str2num(NIRS.Dt.fir.pp(fstep).job.lowcutfreq);
        else
            lowcut=job.fparameters.lowcut;
            applylowcut=job.fparameters.applylowcut;
        end
        
        if contains(num2str(job.fparameters.highcut),'same') || contains(num2str(job.fparameters.applyhighcut),'same')
            [highcut ,applyhighcut] = str2num(NIRS.Dt.fir.pp(fstep).job.highcutfreq);
        else
            highcut=job.fparameters.highcut;
            applyhighcut=job.fparameters.applyhighcut;
        end
        
        if  contains(num2str(job.fparameters.paddingsym),'same')
            paddingsym = NIRS.Dt.fir.pp(fstep).job.paddingsymfilter; %symetrie padding on the signal to avoid edge on the filtering
        else
            paddingsym =job.fparameters.paddingsym;
        end
        
        if  contains(num2str(job.fparameters.interpolate),'same')
            interpolate = NIRS.Dt.fir.pp(fstep).job.interpolatebadfilter;  %interpolate bad intervals
        else
            interpolate = job.fparameters.interpolate;
        end
        
        if contains(num2str(job.fparameters.filt_ord),'same')
            filt_ord = NIRS.Dt.fir.pp(fstep).job.filterorder;
        else
            filt_ord = job.fparameters.filt_ord;
        end
    end
    
    %new parameters for interpolation 
    if isfield(job,'newparameters') %job.newparameters.interp_mode=  ;
        if isfield(job.newparameters,'interp_mode') 
            %1= for linear fit between start and end of bad segment
            %2= replace bad segments with the average of the block (good
            %portion without artifact)
            interp_mode=job.newparameters.interp_mode;
        else
            interp_mode=1; 
        end
    else
        interp_mode=1;
    end
        
    %use last step of preprocessing
    lst = length(NIRS.Dt.fir.pp);
    rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
    NC = NIRS.Cf.H.C.N; %number channels
    fs = NIRS.Cf.dev.fs; %sampling rate
    
    if ~isfield(NIRS.Dt,'AUX')
        error('No auxiliary attached to the NIRS file %s\nCheck in your NIRS.mat file... is there a <AUX> field in NIRS.Dt?\n', (job.NIRSmat{1}))
    end
    %% get AUX data + load for each block
    for iAUX = 1:numel(NIRS.Dt.AUX)
        newAUX=numel(NIRS.Dt.AUX)+1;
        if contains(NIRS.Dt.AUX(iAUX).label,'AUX') || contains(NIRS.Dt.AUX(iAUX).label,'EEG')
            fprintf('AUX found in NIRS.mat. Creating new AUX .dat file (downsampled, filtered, normalized)...\n')
            if ~isfield(NIRS.Dt.AUX(iAUX).pp(end),'sync_timesec') %check if there were synchronisation -take the last field
                disp('No segmentation have been made -- ensure that aux synchronisation are ok')
                disp('FilterAUX not completed')
            else
                
                %%%%%%%%%%plot
                figg=figure('units','normalized','outerposition',[0 0 1 1]);
                sizefig=ceil(length(NIRS.Dt.AUX(iAUX).pp(end).p)/4);
                figg=tiledlayout(7,4,'TileSpacing','Compact','Padding','Compact');
                nfig=1;
                %%%%%%%%%%plot
                
                for i=1:length(NIRS.Dt.AUX(iAUX).pp(end).p)
                    %%%%%%%%%%plot
                    if i==(7*4+1) || i==(7*4*2+1)
                        saveas(figg,[cPATH filesep outAUXfolder filesep 'AUXfigure_' cFILE '_' num2str(nfig) '.fig'])
                        saveas(figg,[cPATH filesep outAUXfolder filesep 'AUXfigure_' cFILE '_' num2str(nfig) '.png'])
                        close
                        clear figg
                        figg=figure('units','normalized','outerposition',[0 0 1 1]);
                        sizefig=ceil(length(NIRS.Dt.AUX(iAUX).pp(end).p)/4);
                        figg=tiledlayout(7,4,'TileSpacing','Compact','Padding','Compact');
                        nfig=nfig+1;
                    end
                    %%%%%%%%%%plot
                    nexttile;
                    %%%%%%%%%%plot
                    
                    nameAUX=NIRS.Dt.AUX(iAUX).pp(end).p{i};
                    moduleaux = numel(NIRS.Dt.AUX(iAUX).pp);
                    try
                        tstart=NIRS.Dt.AUX(iAUX).pp(end).sync_timesec{i};
                    catch
                        tstart = 0;
                    end
                    
                    timeSECaxe = (1/fs):(1/fs):(NIRS.Dt.fir.sizebloc{i}*1/fs);
                    tstop = tstart+timeSECaxe(end);
                    [cPATH,cFILE,cEXT]=fileparts(nameAUX); %current path file and extension
                    [data,infoBV,marker,auxind_dur_ch] = fopen_EEG(nameAUX, tstart, tstop);
                    
                    fsAUX =1/(infoBV.SamplingInterval/1000000); %Frequence echantillonage Hz
                    [~,q] = rat(fs/fsAUX,0.0001);
                    AUXupdate.ind_dur_ch=auxind_dur_ch((auxind_dur_ch(:,1)>tstart*fsAUX) & (auxind_dur_ch(:,1)<tstop*fsAUX),:);
                    AUXupdate.marker=marker((auxind_dur_ch(:,1)>tstart*fsAUX) & (auxind_dur_ch(:,1)<tstop*fsAUX),:);
                    
                    
                    %% SEGMENT FOR INTERPOLATION: DETERMINE BAD INTERVALS
                    %%BASED ON NIRS DATA. taken from the extractcomponent
                    %%function (physiology) and adapted
                    [nirsdir,nirsname,~] = fileparts(rDtp{i});
                    nirsvmrk_path = fullfile(nirsdir,[nirsname '.vmrk']);
                    [nirsind_dur_ch] = read_vmrk_find(nirsvmrk_path,'bad_step');
                    mrks = [];
                    ind = [];
                    noise =  logical(zeros([size(timeSECaxe,2) NC]));
                    if ~isempty(nirsind_dur_ch) &&  interpolate == 1
                        maxpoint  = nirsind_dur_ch(:,1)+nirsind_dur_ch(:,2);
                        badind = find(maxpoint>size(noise,1));
                        if ~isempty(badind)
                            disp(['Warning file ' nirsvmrk_path ' marker : ' num2str(badind') ' are out of range in the data file'])
                            nirsind_dur_ch(badind,2)=size(noise,2)- nirsind_dur_ch(badind,1);
                        end
                        for Idx = 1:size(noise,2)
                            mrks = find(nirsind_dur_ch(:,3)==Idx);
                            ind = nirsind_dur_ch(mrks,1);
                            indf = ind + nirsind_dur_ch(mrks,2) - 1;
                            if ~isempty(ind)
                                try
                                    for it = 1:numel(ind)
                                        noise(ind(it):indf(it),Idx) = 1;
                                    end
                                catch
                                    msgbox('Noise reading problem')
                                end
                            end
                        end
                        
                        % ici group channel with the same noise latency
                        
                        ind = find((sum(noise,2)./size(noise,2))>0.5);
                        inddiff = diff(ind);
                        
                        if isempty(ind) %No specific noisy event found in file
                            eventbadstartstop = [] ;
                        else
                            idsep = find(inddiff>1);
                            if isempty(idsep)
                                idstart  =[ind(1)];
                                idstop = [ind(end)];
                            else
                                idstart  =[ind(1);ind(idsep(1:end)+1)];
                                idstop =  [ind(idsep(1:end)-1);ind(end)];
                            end
                            
                            % add event start
                            
                            %%start bad event 2 data points before in AUX
                            %%(when possible!)
                            idstart=idstart-2;
                            idstart(idstart<1)=1;
                            eventbadstartstop = [idstart,idstop] ;
                        end
                    end
                    
                    for ich=1:numel(infoBV.name_ele)
                        
                        tmp = data(:,ich);
                        
                        
                        %% downsample to the same as fs NIRS
                        rstmp=downsample( tmp , q);
                        
                        %normalize (z score)
                        rstmp=(rstmp-mean(rstmp))/std(rstmp);
                        
                        %%%%%%%%%%plot
                        plot(rstmp,'color',AUXcolor{ich},'linewidth',.6);hold on;
                        %%%%%%%%%%plot
                        
                        %% INTERPOLATION OF BAD INTERVALS
                        if ~isempty(nirsind_dur_ch)&&  interpolate == 1
                            
                            if interp_mode==1 %linear fit
                                %eventbadstartstop = [idstart,idstop]
                                for sizebad = 1:size(eventbadstartstop,1)
                                    dur = eventbadstartstop(sizebad,2)-(eventbadstartstop(sizebad,1))+1;
                                    y1=rstmp((eventbadstartstop(sizebad,1)));
                                    y2=rstmp(eventbadstartstop(sizebad,2));
                                    interval = (eventbadstartstop(sizebad,1)):eventbadstartstop(sizebad,2);
                                    %y = ax + b
                                    a = (y2-y1)/dur;
                                    b = y1 - a*(eventbadstartstop(sizebad,1));
                                    interp = a.*interval + b;
                                    rstmp(interval) = interp;
                                end
                            elseif interp_mode==2 %average
                                
                                for sizebad = 1:size(eventbadstartstop,1)
                                    dur = eventbadstartstop(sizebad,2)-(eventbadstartstop(sizebad,1))+1;
                                    y1=rstmp((eventbadstartstop(sizebad,1)));
                                    y2=rstmp(eventbadstartstop(sizebad,2));
                                    interval = (eventbadstartstop(sizebad,1)):eventbadstartstop(sizebad,2);
                                    goodinterval=rstmp([1:eventbadstartstop(sizebad,1) eventbadstartstop(sizebad,2):end]);
                                    
                                    rstmp(interval) = mean(goodinterval);
                                end
                                
                            end
                        end
                        
                        %%%%%%%%%%plot
                        plot(rstmp,'color',AUXcolor{ich}+.1,'linewidth',.6);hold on;
                        %%%%%%%%%%plot
                        
                        %% FILTER
                        if paddingsym
                            
                            d = [fliplr(rstmp') rstmp' fliplr(rstmp')];
                            d=d';
                            tstartd = size(rstmp,1)+1;
                            tstopd = size(rstmp,1)*2;
                        else
                            d=rstmp;
                        end
                        
                        id = find(isnan(d));
                        if ~isempty(id)
                            d(id) = mean(d,'omitnan');
                        end
                        
                        %bandpass
                        %                         if applylowcut && applyhighcut %band pass (low and high)
                        %                             W2 = lowcut*2/fs;
                        %                             W1 = highcut*2/fs;
                        %                             [fb,fa]=butter(filt_ord,[W1, W2]);
                        %                             dfilt = filtfilt(fb,fa,d);
                        if applylowcut==true %&& applyhighcut==false %only low pass
                            Wn = lowcut*2/fs;
                            [fb,fa]=butter(filt_ord,Wn);
                            dfilt1 = filtfilt(fb,fa,d);
                        else
                            dfilt1=d;
                        end
                        if applyhighcut==true % && applylowcut==false %only high pass
                            Wn = highcut*2/fs;
                            [fb,fa]=butter(filt_ord,Wn,'high');
                            dfilt = filtfilt(fb,fa,dfilt1);
                        else
                            dfilt=dfilt1;
                        end
                        
                        if paddingsym
                            dfilt = dfilt(tstartd:tstopd);
                        end
                        
                        
                        %END OF FILTER SECTION
                        
                        %%%%%%%%%%plot
                        plot(dfilt,'color',AUXcolor{ich}+.2,'linewidth',1.5);hold on;
                        
                        xlim([0 numel(timeSECaxe)])
                        %%%%%%%%%%plot
                        tmpn=dfilt;
                        
                        % we cut aux to the data initial size
                        if numel(tmpn)<numel(timeSECaxe)
                            nplus = numel(timeSECaxe)-numel(tmpn);
                            try
                                tmpn = [tmpn; tmpn(end-nplus:end) ];
                            catch
                                try
                                    tmpn = [tmpn, tmpn(end-nplus:end) ];
                                catch
                                    
                                    msgbox(sprintf('Too short AUX\n%s b%d',cFILE,i))
                                end
                            end
                        elseif numel(tmpn)>numel(timeSECaxe)
                            tmpn = tmpn(1:numel(timeSECaxe));
                        end
                        
                        dataNEW(:,ich)=tmpn;
                    end
                    
                    %%%%%%%%%%plot
                    if sum(NIRS.Cf.H.C.ok(:,i),'all')<=10
                        text(.95,.9,['bad block #' num2str(i)],'HorizontalAlignment','right','units','normalized')
                    else
                        text(.95,.9,['block #' num2str(i)],'HorizontalAlignment','right','units','normalized')
                    end
                    yytemp=ylim;
                    if yytemp(1)<-5 && yytemp(2)>5
                        ylim([-5 5]);
                    elseif yytemp(1)<-5
                        ylim([-5 yytemp(2)]);
                    elseif yytemp(2)>5
                        ylim([yytemp(1) 5]);
                    end

                    %%%%%%%%%%plot
                    
                    %OVERWRITE INFOS in infoBV
                    infoBV.DataPoints=size(dataNEW,1);
                    infoBV.SamplingInterval=(1/fs)*1000000;
                    
                    AUXupdate.infoBV=infoBV;
                    AUXupdate.data=dataNEW;
                    
                    if ~exist([cPATH filesep outAUXfolder filesep],'dir')
                        mkdir([cPATH filesep outAUXfolder filesep])
                    end
                    fileoutAUX=[cPATH filesep outAUXfolder filesep cFILE 'b' num2str(i) cEXT];
                    fwrite_EEG(fileoutAUX,AUXupdate,1,AUXupdate.infoBV.DataPoints );
                    disp(fileoutAUX)
                    NIRS.Dt.AUX(newAUX).pp(moduleaux).p{i,1}=fileoutAUX;
                    NIRS.Dt.AUX(newAUX).pp(moduleaux).sync_timesec{i,1}=0;
                end
                NIRS.Dt.AUX(newAUX).label = ['fil' NIRS.Dt.AUX(iAUX).label ];
                
                %%%%%%%%%%plot savefig
                
                saveas(figg,[cPATH filesep outAUXfolder filesep 'AUXfigure_' cFILE '_' num2str(nfig) '.fig'])
                saveas(figg,[cPATH filesep outAUXfolder filesep 'AUXfigure_' cFILE '_' num2str(nfig) '.png'])
                close
                clear figg
                %%%%%%%%%%plot
                
            end
        end
    end
    
    fprintf('Update NIRSmat COMPLETED ...%s\n*\n**\n***\n',job.NIRSmat{1})
    save(job.NIRSmat{1},'NIRS');
    
end
out.NIRSmat = job.NIRSmat;
end
