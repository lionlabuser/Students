clear;close all
AUXpath='C:\Users\laura\AnalysesNIRS\MART0m\Multimodal\BB032';
AUXfile='BB032_0m_hist_AUX.dat';
idpart=split(AUXpath,filesep);
generalmultipath=fullfile(idpart{1:end-1}); %to get the general multimodal folder
idpart=idpart{end}; %to get participant ID

%open Following AUX correction file
allupdatefile=[generalmultipath filesep 'ALL_updateAUX_artefactcorrection.mat'];
if exist(allupdatefile,'file')
    load(allupdatefile,'out')
    outrow=length(out)+1;
else
    out=struct;
    outrow=1;
end

%import EEG file
[data,infoBV,marker,ind_dur_ch] = fopen_EEG([AUXpath filesep AUXfile]);

%find SAT and RESP
for id=1:length(infoBV.name_ele)
    if strcmpi(infoBV.name_ele{id}, 'SAT')
        dataid(1)=id;
    elseif strcmpi(infoBV.name_ele{id}, 'RESP')
        dataid(2)=id;
    end
end
chanlabels=infoBV.name_ele(dataid);
fsAUX=1/(infoBV.SamplingInterval/1000000);
fs=7.8125;
[~,q] = rat(fs/fsAUX,0.0001); %ratio sampling rate EEG (500Hz) vs nirs (7.8125)

%adjust the eeg data info to fit the new sampling rate!
ind_dur_ch(2:end,1)=ind_dur_ch(2:end,1)./q; %ajust trigs
infoBV.SamplingInterval=1/7.8125*1000000;
infoBV.DataPoints=infoBV.DataPoints/q;

for c=1:2
    newdata(:,c)=downsample(data(:,dataid(c)),q); %downsample to fit nirs
    newdata(:,c)=zscore(newdata(:,dataid(c))); %normalized
end

threshold=[.3 .7;... %moving STD threshold, respectively for SAT and RESP, in standardized units    .34 .67
    .4   .8];%...  %moving AVG threshold, respectively for SAT and RESP, in standardized units    .4   1.0
%    .2 .2 ]; %moving STD DIFFERENCE threshold, respectively for SAT and RESP, in standardized units    .34 .67
adjwindow=[1 3] ; %moving window length before and after X in seconds ( *2: if you want a 3 sec window = please write 1.5)
minlength=3; %minimum length of intervals (good or bad) in data points
pausetime=5;
question1=[0 0];
while any(question1==0)
    clear xma
    xma{1}=1;xma{2}=1;
    for x=1:size(newdata,1)
        
        tmpwindow=[x-ceil(fs*3.5) x+ceil(fs*3.5)];
        if tmpwindow(1) < 1
            tmpd=1-tmpwindow(1);
            tmpwindow=tmpwindow+tmpd;
        elseif tmpwindow(2) > size(newdata,1)
            tmpd=size(newdata,1)-tmpwindow(2);
            tmpwindow=tmpwindow+tmpd;
        end
        
        for c=1:2
            mstd(x,c)=std(newdata(tmpwindow(1):tmpwindow(2),c))  ;
            mavg(x,c)=mean(newdata(tmpwindow(1):tmpwindow(2),c))  ;
            % if x==1; mdstd(x,c)=0;
            % else; mdstd(x,c)=abs(diff([ mstd(x-1,c) mstd(x,c)])); end
            
            %identify which periods are considered an interval if at least one threshold is overpassed
            if mstd(x,c)>=threshold(1, c) || abs(mavg(x,c))>=threshold(2,c) %|| mdstd(x,c)>=threshold(3, c)
                MA(x,c)=true;
            else
                MA(x,c)=false;
            end
            
            if x>1 && diff([MA(x-1,c) MA(x,c)])
                xma{c}=[xma{c} x]; %identify data points == start and end of bad intervals
                if diff(xma{c}(end-1:end))<minlength %if the time points are too close in time (min length of interval, good or bad)
                    xma{c}(end-1:end)=[]; %merge with previous
                end
            end
            
        end
    end
    
    
    for c=1:2
        if ~question1(c)
            figaux=figure;
            figaux.Units='normalized';figaux.Position=[.01 .1 .97 .8]; figaux.Color=[1 1 1];
            figaux=tiledlayout(3,1,'TileSpacing','Compact','Padding','Compact'); %multi plot in a same figure
            
            nexttile;
            area(1:length(newdata(:,c)),MA(:,c)*10,0,'LineStyle','none','FaceColor',[.7 .7 .7]);hold on %box
            plot(newdata(:,c)); ylabel('AUX DATA')
            ylim([-5 5])
            nexttile;
            area(1:length(newdata(:,c)),MA(:,c)*10,0,'LineStyle','none','FaceColor',[.7 .7 .7]);hold on %box
            plot(mstd(:,c)); ylim([0 1.2]); yline(threshold(1,c)); ylabel('Moving STD')
            %nexttile;
            %area(1:length(newdata(:,c)),MA(:,c)*10,0,'LineStyle','none','FaceColor',[.7 .7 .7]);hold on %box
            %plot(mdstd(:,c)); ylim([0 1.2]); yline(threshold(3,c)); ylabel('Moving STD changes')
            nexttile;
            area(1:length(newdata(:,c)),MA(:,c)*10,0,'LineStyle','none','FaceColor',[.7 .7 .7]);hold on %box
            plot(mavg(:,c)); ylim([-1.5 1.5]); yline(threshold(2,c)); ylabel('Moving AVG');xlabel('Time (in data pts at 7.8Hz)')
            
            title(figaux,['Artifact Detection for : ' chanlabels{c}])
            minn=0:2000:length(newdata(:,c))-2000;
            maxx=minn+2000;
            for xx=1:length(minn)
                for a=1:size(figaux.Children,1)
                    nexttile(a); xlim([minn(xx) maxx(xx)]);
                end
                pause(pausetime)
            end
            question1(c)=input(['Are you satisfied with the identification of bad intervals for the current aux (' chanlabels{c} ')? [1=yes,0=no] ']);
            if ~question1(c)
                threshold(1,c)=input(['New threshold for moving STD? (previous = ' num2str(threshold(1,c)) '). Enter here: ']);
                %threshold(3,c)=input(['New threshold for moving STD changes? (previous = ' num2str(threshold(3,c)) '). Enter here: ']);
                threshold(2,c)=input(['New threshold for moving AVG? (previous = ' num2str(threshold(2,c)) '). Enter here: ']);
                adjwindow(c)=input(['New duration for moving window (in sec)? (previous = ' num2str(adjwindow(c)) '). Enter here: ']);
                pausetime=input(['The time window was changing every ' num2str(pausetime) ' sec. How many seconds do you want to see the data?. Enter here: ']);
                
            else
                savefig([AUXpath filesep 'ArtifactDetection_' chanlabels{c} '.fig'])
            end
            close
            clear figaux
        end
        
    end
end
fprintf('=========\nOut of first loop: Artifact identification completed.\n=========\n')
out(outrow).part=idpart;
out(outrow).originalfile=[AUXpath filesep AUXfile];
out(outrow).corrfile=[AUXpath filesep 'c' AUXfile];
out(outrow).channels=chanlabels;
out(outrow).mSTD_threshold=threshold(1,:);
%out(outrow).mdSTD_threshold=threshold(3,:);
out(outrow).mAVG_threshold=threshold(2,:);
out(outrow).timewindow=adjwindow;
out(outrow).minlength=minlength;

for c=1:2 %segmentation
    for m=1:length(xma{c})
        if m==length(xma{c})
            dsegment{m,c}= newdata(xma{c}(m):length(newdata(:,c)),c);
        else
            dsegment{m,c}= newdata(xma{c}(m):(xma{c}(m+1)-1),c);
            
        end
    end
end
pausetime=3;
int=[6 1]; %xaxis multiplier
pp=[.02 .01]; %cubic smoothing spline (pp=smoothing parameter);
%close to 0= least-square straight line fit ;
%close to 1= natural cubic spline interpolant;
question2=[0 0];
while any(question2==0)
    close
    clear figcorr
    for c=1:2 %spline interpolation
        for m=1:length(xma{c})
            if MA(xma{c}(m),c)
                idsegment{m,c}='bad';
                tmpsegment=csaps((1:length(dsegment{m,c}))*int(c),dsegment{m,c}',pp(c),(1:length(dsegment{m,c}))*int(c));
                %same function used by Scholkmann et al 2010 doi:
                % for more info: https://www.mathworks.com/help/curvefit/csaps.html
                newsegment{m,c}=dsegment{m,c}-(tmpsegment');
            else
                idsegment{m,c}='good';
                newsegment{m,c}=dsegment{m,c};
            end
        end
    end
    
    figcorr=figure;
    figcorr.Units='normalized';figcorr.Position=[.01 .1 .97 .8]; figcorr.Color=[1 1 1];
    figcorr=tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact'); %#ok<*NASGU> %multi plot in a same figure
    
    for c=1:2
        nexttile;
        corrdata(:,c)=vertcat(newsegment{1:length(xma{c}),c}); %#ok<*SAGROW>
        plot(newdata(:,c));hold on; plot(corrdata(:,c));
        ylabel(chanlabels{c})
    end
    title(figcorr,'Before and after data correction by spline interpolation')
    minn=0:1000:length(newdata(:,c))-999;
    maxx=minn+1000;
    for xx=1:length(minn)
        for a=1:size(figcorr.Children,1)
            nexttile(a); xlim([minn(xx) maxx(xx)]);
        end
        pause(pausetime)
    end
    
    for c=1:2
        
        question2(c)=input(['Are you satisfied with the identification of bad intervals for ' chanlabels{c} ' aux? [1=yes,0=no] Enter here: ']);
        if ~question2(c)
            int(c)=input(['New x-axis multiplier (previous = ' num2str(int(c)) '). Enter here: ']);
            pp(c)=input(['New smoothing parameter? (previous = ' num2str(pp(c)) '). Enter here: ']);
            pausetime=input(['The time window was changing every ' num2str(pausetime) ' sec. How many seconds do you want to see the data?. Enter here: ']);
        end
    end
end
savefig([AUXpath filesep 'ArtifactCorrection_' chanlabels{c} '.fig'])
close
clear figcorr

out(outrow).splineparameter=pp;
out(outrow).splineXfactor=int;

%SAVE NEW AUX (start with a 'c' for corrected, in the same folder as the previous!)
AUXupdate.infoBV=infoBV;
AUXupdate.data=corrdata;
AUXupdate.ind_dur_ch=ind_dur_ch;
AUXupdate.marker=marker;
fileoutAUX=out(outrow).corrfile;
fwrite_EEG(fileoutAUX,AUXupdate,1,AUXupdate.infoBV.DataPoints );
disp(fileoutAUX)

%save update file for following all participants' AUX artifact correction
save(allupdatefile,'out');


%% UPDATE INFO IN NIRS.MAT
nirsmat=['C:\Users\laura\AnalysesNIRS\MART0m\' idpart '\Segment\NormV2\NIRS.mat'];

NIRS = [];
load(nirsmat);


for iAUX = 1:numel(NIRS.Dt.AUX)
    newAUX=numel(NIRS.Dt.AUX)+1;
    if contains(NIRS.Dt.AUX(iAUX).label,'AUX') || contains(NIRS.Dt.AUX(iAUX).label,'EEG')
        fprintf('AUX found in NIRS.mat. Changing name of AUX to get the corrected AUX...\n')
        
        for i=1:length(NIRS.Dt.AUX(iAUX).pp(end).p)
            
            [cPATH,cFILE,cEXT]=fileparts(NIRS.Dt.AUX(iAUX).pp(end).p{i});
            NIRS.Dt.AUX(iAUX).pp(end).p{i}=fullfile(cPATH,['c' cFILE cEXT]);
        end
    end
end
save(nirsmat,'NIRS');
fprintf('Update+save completed: %s\n',nirsmat)

