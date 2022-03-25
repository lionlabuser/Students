%___________________________________________________________________
% Copyright (C) 2019 LION Lab, Centre de recherche CHU Sainte-Justine
% www.lionlab.umontreal.ca
%___________________________________________________________________
function out = nirs_run_AUX_artifactDetectionCorrection(job)

for filenb=1:size(job.NIRSmat,1)
    NIRS = [];
    load(job.NIRSmat{1});    
    %clear;close all
    idpart = job.ID;
    disp(['Computing AUX Artifact for ' idpart])
    AUXlabel = job.AUX;
    threshold = job.threshold;
    fsNIRS = NIRS.Cf.dev.fs; %sampling rate
    
    if ~isfield(NIRS.Dt,'AUX')
        error('No auxiliary attached to the NIRS file %s\nCheck in your NIRS.mat file... is there a <AUX> field in NIRS.Dt?\n', (job.NIRSmat{1}))
    end
    
    if max(contains({NIRS.Dt.AUX.label}, 'corr'))
        answer = questdlg('Corrected AUX already computed, what do you want to do?','Warning','Overwrite','Skip Participant','Skip Participant');
        switch answer
            case 'Overwrite'
                corroverwrite = 1;
                icorrAUX = find(contains({NIRS.Dt.AUX.label}, 'corr'));
                disp(['The old corrected AUX will be overwritten'])
            case 'Skip Participant'
                return
        end
    end
    
    if ~exist('corroverwrite','var')
            newAUX = numel(NIRS.Dt.AUX)+1;
        else
            newAUX = icorrAUX;
    end
    
    if numel(NIRS.Dt.AUX)>1
        listAUX = {NIRS.Dt.AUX.label};
        [indx, tf] = listdlg('PromptString',{'Multiple AUX version in the NIRSmat file.',...
            'Select on which one you want to apply the artifact detection.',''},'SelectionMode','single','ListString',listAUX);
        usedAUX = indx;
        if ~any(tf)
            fprintf('No AUX selected, participant will be skipped\n')
            return
        else
            AUXfile = NIRS.Dt.AUX(usedAUX).pp(1).p{1};
            [AUXpath,AUXname,AUXext] = fileparts(AUXfile);
            AUXname = [AUXname AUXext];
            cAUXpath = [AUXpath filesep 'corrAUX'];
            if ~isfolder(cAUXpath)
                mkdir(cAUXpath)
            end 
       end
    end 

    if ~isfield(NIRS.Dt.AUX(usedAUX).pp(end),'sync_timesec') %check if there were synchronisation -take the last field
        error('No data segmentation has been made -- ensure that aux synchronisation is ok')
    end    
    
    %open Following AUX correction file
    allupdatefile = [cAUXpath filesep 'corrAUX.mat'];
    if exist(allupdatefile,'file')
        load(allupdatefile,'out')
        outrow = length(out)+1;
    else
        out = struct;
        outrow = 1;
    end

    %import EEG file
    try
        [data,infoBV,marker,ind_dur_ch] = fopen_EEG([AUXpath filesep AUXname]);
    catch
        error('AUX not found in the location, consider using Folder Adjustment')
    end
    
    %find SAT, RESP & EKG
    for iAUX = 1:numel(job.AUX)
        idAUX = job.AUX{iAUX};
        for id = 1:length(infoBV.name_ele)
            if strcmpi(infoBV.name_ele{id}, idAUX)
                dataid(iAUX) = id;
                %             elseif strcmpi(infoBV.name_ele{id}, 'RESP')
                %                 dataid(2) = id;
                %             elseif strcmpi(infoBV.name_ele{id}, 'EKG')
                %                 dataid(3) = id;
            end
        end
    end
    %chanlabels = infoBV.name_ele(dataid);
    fsAUX = 1/(infoBV.SamplingInterval/1000000);
    [~,q] = rat(fsNIRS/fsAUX,0.0001); %ratio sampling rate EEG (500Hz) vs nirs (7.8125)
    
    %adjust the eeg data info to fit the new sampling rate!
    ind_dur_ch(2:end,1) = ind_dur_ch(2:end,1)./q; %ajust trigs
    infoBV.SamplingInterval = 1/fsNIRS*1000000;
    infoBV.DataPoints = infoBV.DataPoints/q;
    
    for ich = 1:numel(dataid)
        newdata(:,ich) = downsample(data(:,dataid(ich)),q); %downsample to fit nirs
        newdata(:,ich) = zscore(newdata(:,ich)); %normalized %mod_KR zscore(newdata(:,dataid(ich)))
    end
    
    timeaxe = linspace(0,size(newdata,1)/fsNIRS,size(newdata,1));
    
    %% Part 1 : Identification of artifact
    adjwindow = job.adjwindow ; %moving window length before and after X in seconds ( *2: if you want a 2 sec window = please write 1)
    minlength = round(fsNIRS*job.minlength); %minimum length of intervals (good or bad) in data points = here is 1.5 sec
    pausetime = job.pausetime;
    
    question1 = zeros(1,numel(AUXlabel)); %[0 0];
    xma = num2cell(ones(1,numel(AUXlabel))); %xma{1}=1; xma{2}=1; xma{3}=1;
    for ich = 1:numel(dataid) %for each AUX
       while any(question1(ich)==0)
            %clear xma
            xma{ich} = 1;
            
            %i = 1;
            for x = 1:size(newdata,1) %for each time point
                %for ich = 1:numel(dataid) %for each AUX
                tmpwindow = [x-ceil(fsNIRS*adjwindow(ich)) x+ceil(fsNIRS*adjwindow(ich))]; %take a window X sec before and after the data point
                if tmpwindow(1) < 1 % verifiy if the window extends below 0 sec or over the recording time
                    tmpd = 1 - tmpwindow(1);
                    tmpwindow = tmpwindow + tmpd;
                elseif tmpwindow(2) > size(newdata,1)
                    tmpd = size(newdata,1) - tmpwindow(2);
                    tmpwindow = tmpwindow + tmpd;
                end
                
                mstd(x,ich) = std(newdata(tmpwindow(1):tmpwindow(2),ich)); %compute std dev of the window
                mavg(x,ich) = mean(newdata(tmpwindow(1):tmpwindow(2),ich)); %compute mean of the window
                % if x==1; mdstd(x,c)=0;
                % else; mdstd(x,c)=abs(diff([ mstd(x-1,c) mstd(x,c)])); end
                
                %identify which periods are considered an interval if at least one threshold is overpassed
                if mstd(x,ich)>=threshold(1, ich) || abs(mavg(x,ich))>=threshold(2,ich) %|| mdstd(x,c)>=threshold(3, c)
                    MA(x,ich) = true;
                else
                    MA(x,ich) = false;
                end
                if x>1 && diff([MA(x-1,ich) MA(x,ich)]) %if there is a change in status (artifact or not)
                    xma{ich} = [xma{ich} x]; %identify data points == start and end of bad intervals
                    if diff(xma{ich}(end-1:end)) < minlength %if the time points are too close in time (min length of interval, good or bad)
                        %deleted{i,1} = [xma{ich}(end-1:end)];
                        %i = i + 1;
                        xma{ich}(end-1:end) = []; %merge with previous, delete the two last tp flagged
                        if isempty(xma{ich})
                            xma{ich} = 1;
                        end
                    end
                end
                %end
            end
            
           if ~question1(ich)
                figaux = figure;
                figaux.Units = 'normalized';
                figaux.Position = [.01 .1 .97 .8];
                figaux.Color = [1 1 1];
                figaux = tiledlayout(3,1,'TileSpacing','Compact','Padding','Compact'); %multi plot in a same figure
                
                nexttile;
                area(1:length(newdata(:,ich)),MA(:,ich)*10,0,'LineStyle','none','FaceColor',[.7 .7 .7]); hold on %box
                plot(newdata(:,ich)); ylabel('AUX DATA')
                ylim([-6 6])
                
                nexttile;
                area(1:length(newdata(:,ich)),MA(:,ich)*10,0,'LineStyle','none','FaceColor',[.7 .7 .7]); hold on %box
                plot(mstd(:,ich)); ylim([-3 3]); yline(threshold(1,ich)); ylabel('Moving STD')
                %nexttile;
                %area(1:length(newdata(:,c)),MA(:,c)*10,0,'LineStyle','none','FaceColor',[.7 .7 .7]);hold on %box
                %plot(mdstd(:,c)); ylim([0 1.2]); yline(threshold(3,c)); ylabel('Moving STD changes')
                
                nexttile;
                area(1:length(newdata(:,ich)),MA(:,ich)*10,0,'LineStyle','none','FaceColor',[.7 .7 .7]); hold on %box
                plot(mavg(:,ich)); ylim([-1.75 1.75]); yline(threshold(2,ich)); yline(-threshold(2,ich)); ylabel('Moving AVG');xlabel('Time (in data pts at 7.8Hz)')
                
                title(figaux,[idpart ' Artifact Detection for : ' AUXlabel{ich}])
                minn = 0:2000:length(newdata(:,ich));
                maxx = minn + 2000;
                for xx = 1:length(minn)
                    for a = 1:size(figaux.Children,1)
                        nexttile(a); xlim([minn(xx) maxx(xx)]); %Plot the artifact detection
                    end
                    pause(pausetime)
                end
                
                question1(ich) = input(['Are you satisfied with the identification of bad intervals for the current aux (' AUXlabel{ich} ')? [1=yes,0=no] ']);
                if ~question1(ich)
                    threshold(1,ich) = input(['New threshold for moving STD? (previous = ' num2str(threshold(1,ich)) '). Enter here: ']);
                    %threshold(3,c) = input(['New threshold for moving STD changes? (previous = ' num2str(threshold(3,c)) '). Enter here: ']);
                    threshold(2,ich) = input(['New threshold for moving AVG? (previous = ' num2str(threshold(2,ich)) '). Enter here: ']);
                    % adjwindow(c) = input(['New duration for moving window (in sec)? (previous = ' num2str(adjwindow(c)) '). Enter here: ']);
                    % pausetime = input(['The time window was changing every ' num2str(pausetime) ' sec. How many seconds do you want to see the data?. Enter here: ']);
                    
                else
                    for a = 1:size(figaux.Children,1)
                        nexttile(a); xlim([0 length(newdata(:,ich))]); %Plot the artifact detection
                    end
                    savefig([cAUXpath filesep idpart '_ArtifactDetection_' AUXlabel{ich} '.fig'])
               end
                close
                clear figaux
           end
       end
    end
    fprintf('Artifact identification completed.\n')
    out(outrow).part = idpart;
    out(outrow).originalfile = [AUXpath filesep AUXname];
    out(outrow).corrfile = [cAUXpath filesep 'corr' AUXname];
    out(outrow).channels = AUXlabel;
    out(outrow).mSTD_threshold = threshold(1,:);
    %out(outrow).mdSTD_threshold=threshold(3,:);
    out(outrow).mAVG_threshold = threshold(2,:);
    out(outrow).timewindow = adjwindow;
    out(outrow).minlength = minlength;
    
    %% Part 2: segmentation into good/bad
    for ich = 1:numel(dataid) %for each AUX
        for m = 1:length(xma{ich}) %for each interval
            if m == length(xma{ich})
                dsegment{m,ich} = newdata(xma{ich}(m):length(newdata(:,ich)),ich); %prendre du timepoint à la fin des données
            else
                dsegment{m,ich} = newdata(xma{ich}(m):(xma{ich}(m+1)-1),ich); %prendre du timepoint au prochain changement d'interval
                
            end
        end
    end
    
    %% Part 3: correction
    pausetime = job.pausetime; %time for wiewing before the automatic change of the time window
    int = job.int; %xaxis multiplier (because SAT has really high frequency artifacts,
    % I tried and the algorithm was better to correct the SAT artifacts when we multiply the Xaxis
    pp = job.pp; %cubic smoothing spline (pp=smoothing parameter);
    %close to 0= least-square straight line fit ;
    %close to 1= natural cubic spline interpolant;
   question3 = input('Do you need to do an offset adjustment after each bad interval? [1=yes,0=no] Enter here: '); %offset adjustment
   question4 = question3;  %offset adjustment
    question2 = zeros(1,numel(AUXlabel)); %[0 0];
    
    
   while any(question2 == 0) || ~(question3 == question4) %if insatisfied or change your mind about offset adjustment
        close
        clear figcorr
        for ich = 1:numel(dataid)%for each AUX
            for m = 1:length(xma{ich}) %for each interval
                if MA(xma{ich}(m),ich) %if marked as artifact
                    idsegment{m,ich} = 'bad'; %label as bad
                    tmpsegment = csaps(timeaxe(1:length(dsegment{m,ich}))*int(ich),dsegment{m,ich}',pp(ich),timeaxe(1:length(dsegment{m,ich}))*int(ich)); %spline interpolation
                    %same function used by Scholkmann et al 2010 doi:
                    % for more info: https://www.mathworks.com/help/curvefit/csaps.html
                    newsegment{m,ich} = dsegment{m,ich} - (tmpsegment'); %substract interpolated segment from original segment
                    if m>1 && question4 %offset adjustment
                        tmpdiff = mean(newsegment{m-1,ich}(end-minlength+1:end)) - mean(newsegment{m,ich}(1:minlength-1)); %difference between the beggining of the current segment and the end of the previous segment
                        newsegment{m,ich} = newsegment{m,ich} + tmpdiff; %add the difference to the current segment
                    end
                else
                    idsegment{m,ich} = 'good'; %label as good
                    newsegment{m,ich} = dsegment{m,ich}; %keep the original segment in the newdata
                    if m>1 && question4 %offset adjustment
                        tmpdiff = mean(newsegment{m-1,ich}(end-minlength+1:end)) - mean(newsegment{m,ich}(1:minlength-1)); %difference between the beggining of the current segment and the end of the previous segment
                        newsegment{m,ich} = newsegment{m,ich} + tmpdiff; %add the difference to the current segment
                    end
                end
            end
        end
        
        figcorr = figure;
        figcorr.Units = 'normalized';
        figcorr.Position=[.01 .1 .97 .8];
        figcorr.Color=[1 1 1];
        figcorr = tiledlayout(numel(dataid),1,'TileSpacing','Compact','Padding','Compact'); %#ok<*NASGU> %multi plot in a same figure
        
        for ich = 1:numel(dataid) %for each AUX
            nexttile;
            newcorrdata(:,ich) = vertcat(newsegment{1:length(xma{ich}),ich}); %#ok<*SAGROW> %concatenate all the new segments
            plot(newdata(:,ich)); ylim([-6 6]); hold on;
            plot(newcorrdata(:,ich)); ylim([-6 6]);
            legend({'Non-corrected signal' 'Corrected signal'})
            ylabel(AUXlabel{ich})
            
        end
        title(figcorr,[idpart ' Before and after data correction by spline interpolation'])
        minn = 0:1000:length(newdata(:,ich));
        maxx = minn+1000;
        for xx = 1:length(minn)
            for a = 1:numel(dataid)
                nexttile(a); xlim([minn(xx) maxx(xx)]); %plot artifact corrected and original AUX
            end
            pause(pausetime)
        end
        
        for ich = 1:numel(dataid) %for each AUX
            question2(ich) = input(['Are you satisfied with the correction of bad intervals for ' AUXlabel{ich} ' aux? [1=yes,0=no] Enter here: ']);
            if ~question2(ich) %if not satisfied
                int(ich) = input(['New x-axis multiplier (previous = ' num2str(int(ich)) '). Enter here: ']);
                %if ich == 2
                pp(ich) = input(['New smoothing parameter? (previous = ' num2str(pp(ich)) '). Enter here: ']);
                %end
                %  pausetime=input(['The time window was changing every ' num2str(pausetime) ' sec. How many seconds do you want to see the data?. Enter here: ']);
            end
        end
        if any(question2==0)
            question4 = input('Do you need to do an offset adjustment after each bad interval? [1=yes,0=no] Enter here: ');
        end
        if question2(ich) && ~question4
            for a = 1:numel(dataid)
                nexttile(a); xlim([0 length(newdata(:,ich))]); %plot artifact corrected and original AUX
            end
            savefig([cAUXpath filesep idpart '_ArtifactCorrection.fig']) %AUXlabels{ich}
       end
        close
        clear figcorr
    end
    
    fprintf('Artifact correction completed.\n')
    
    out(outrow).splineparameter = pp;
    out(outrow).splineXfactor = int;
    
    
    if size(data,2) ~= size(newdata,2) %if number of AUX corrected different than total number of AUX
        corrdata = zeros(size(newdata,1),size(data,2));
        for n = 1:size(data,2) %for each AUX
            if any(dataid(:) == n) %if AUX was corrected
                corrdata(:,n) = newcorrdata(:,find(dataid(:) == n));
            else %if AUX wasn't corrected
                tmpdata = downsample(data(:,n),q); %downsample to fit nirs
                tmpdata = zscore(tmpdata(:,n)); %normalized %mod_KR zscore(newdata(:,dataid(ich)))
                corrdata(:,n) = tmpdata(:,n);
            end
        end
    else
        corrdata = newcorrdata(:,[dataid]); %put again in initial order
    end
    
    %SAVE NEW AUX (starting with 'corr' for corrected, in the same folder as the previous!)
    AUXupdate.infoBV = infoBV;
    AUXupdate.data = corrdata;
    AUXupdate.ind_dur_ch = ind_dur_ch;
    AUXupdate.marker = marker;
    fileoutAUX = out(outrow).corrfile;
    fwrite_EEG(fileoutAUX,AUXupdate,1,AUXupdate.infoBV.DataPoints );
    fprintf('Save of corrAUX completed: %s\n',fileoutAUX)
    
    %save update file for following all participants' AUX artifact correction
    save(allupdatefile,'out');
    clearvars out
    %end
    
    %% UPDATE INFO IN NIRS.MAT
    % because some participants have 2 AUX files, it might be better to do this
    % step separately. This is why it is now only in <comments>. But the script
    % below works really well :) And you can use it!
    
    for i = 1:length(NIRS.Dt.AUX(usedAUX).pp(end).p) %pour chaque path listé dans la dernière entrée de path pour chaque version du AUX == bloc?
        moduleaux = numel(NIRS.Dt.AUX(usedAUX).pp);
        NIRS.Dt.AUX(newAUX).pp(moduleaux).p{i,1} = fileoutAUX;
        NIRS.Dt.AUX(newAUX).pp(moduleaux).sync_timesec{i,1} = NIRS.Dt.AUX(usedAUX).pp(end).sync_timesec{i}; %0
    end
    NIRS.Dt.AUX(newAUX).label = ['corr' NIRS.Dt.AUX(usedAUX).label];

    save(job.NIRSmat{filenb},'NIRS');
    fprintf('Update of NIRSmat completed: %s\n',job.NIRSmat{filenb})
end
%out.NIRSmat = job.NIRSmat{1};

%     function [output,p] = csaps(x,y,p,xx,w)
%         %CSAPS Cubic smoothing spline.
%         %
%         %   CSAPS(X,Y)  returns the ppform of a cubic smoothing spline for the
%         %   given data X,Y. The smoothing spline approximates, at the data site X(j),
%         %   the given data value Y(:,j), j=1:length(X). The data values may be
%         %   scalars, vectors, matrices, or even ND-arrays. Data points with the
%         %   same site are replaced by their (weighted) average, and this may affect the
%         %   smoothing spline.
%         %   For smoothing of gridded data, see below.
%         %
%         %   This smoothing spline  f  minimizes
%         %
%         %   P * sum_j W(j) |Y(:,j) - f(X(j))|^2  +  (1-P) * integral |D^2 f|^2 .
%         %
%         %   Here, the sum is over j=1:length(X);
%         %   X  and  Y  are the result of replacing any data points with the same site
%         %   by their weighted average, with its weight the sum of the corresponding
%         %   weights;
%         %   the integral is taken over the interval [min(X) .. max(X)];
%         %   |z|^2  is the sum of the squares of the entries of  z ;
%         %   D^2 f  is the second derivative of the function  f ;
%         %   W = ones(length(X),1)  is the default value for W; and
%         %   the default value for the smoothing parameter P is chosen,
%         %   in an ad hoc fashion and in dependence on X, as indicated in
%         %   the next paragraph. You can supply a specific value for P, by using
%         %   CSAPS(X,Y,P)  instead.
%         %
%         %   When P is 0, the smoothing spline is the least-squares straight line fit
%         %   to the data, while, at the other extreme, i.e., when P is 1, it is the
%         %   `natural' or variational cubic spline interpolant. The transition region
%         %   between these two extremes is usually only a rather small range of values
%         %   for P and its location strongly depends on the data sites. It is in this
%         %   small range that P is chosen when it is not supplied, or when an empty
%         %   P or a negative P is input.
%         %   If P > 1 , the corresponding solution of the above minimization problem
%         %   is returned, but this amounts to a roughening rather than a smoothing.
%         %
%         %   If the resulting smoothing spline, pp, is to be evaluated outside its basic
%         %   interval, it should be replaced by fnxtr(pp) to ensure that its second
%         %   derivative is zero outside that interval.
%         %
%         %   [OUT,P] = CSAPS(X,Y,...)  returns the value of P actually used, and this
%         %   is particularly useful when no P or an empty P was specified.
%         %
%         %   If you have difficulty choosing P but have some feeling for the size
%         %   of the noise in Y, consider using instead  spaps(X,Y,tol)  which, in
%         %   effect, chooses P in such a way that the roughness measure,
%         %                integral (D^2 f)^2 ,
%         %   is as small as possible subject to the condition that the error measure,
%         %                sum_i W(i)(Y(i) - f(X(i)))^2 ,
%         %   does not exceed the specified  tol . This usually means that the error
%         %   measure equals the specified  tol .
%         %
%         %   CSAPS(X,Y,P,XX)  returns the value(s) at XX of the cubic smoothing spline,
%         %   unless XX is empty, in which case the ppform of the cubic smoothing
%         %   spline is returned. This latter option is important when the user wants
%         %   the smoothing spline (rather than its values) corresponding to a specific
%         %   choice of error weights, as is discussed next.
%         %
%         %   CSAPS(X,Y,P,XX,W)  returns, depending on whether or not XX is empty, the
%         %   ppform, or the values at XX, of the cubic smoothing spline for the
%         %   specified weights W. Any negative weight is replaced by 0, and that
%         %   makes the resulting smoothing spline independent of the corresponding
%         %   data point. When data points with the same site are averaged, their
%         %   weights are summed.
%         %
%         %   See below for the case of GRIDDED data.
%         %
%         %   Example:
%         %      x = linspace(0,2*pi,21); y = sin(x)+(rand(1,21)-.5)*.1;
%         %      pp = csaps(x,y, .4, [], [ones(1,10), repmat(5,1,10), 0] );
%         %   returns a smooth fit to the data that is much closer to the data
%         %   in the right half, because of the much larger weight there, except for
%         %   the last data point, for which the weight is zero.
%         %
%         %   It is also possible to vary the smoothness requirement, by having P be a
%         %   sequence (of the same length as X) rather than a scalar.
%         %   In that case, the roughness measure is taken to be
%         %                   integral lambda(t)*(D^2 f(t))^2 dt ,
%         %   with the roughness weight  lambda  the piecewise constant function with
%         %   interior breaks X whose value on the interval (X(i-1),X(i)) is P(i)  for
%         %   i=2:length(x),  while P(1) continues to be taken as the smoothing
%         %   parameter, P.
%         %
%         %   Example:
%         %      pp1 = csaps(x,y, [.4,ones(1,10),repmat(.2,1,10)], [], ...
%         %                                  [ones(1,10), repmat(5,1,10), 0]);
%         %   uses the same data, smoothing parameter, and error weight as in the
%         %   earlier example, but chooses the roughness weight to be only .2 in the
%         %   right half of the interval and gives, correspondingly, a rougher but
%         %   better fit there, -- except for the last data point which is ignored.
%         %   A plot showing both examples for comparison could now be obtained by
%         %      fnplt(pp); hold on, fnplt(pp1,'r'), plot(x,y,'ok'), hold off
%         %      title('cubic smoothing spline, with right half treated differently:')
%         %      xlabel(['blue: larger error weights; ', ...
%         %              'red: also smaller roughness weights'])
%         %
%         %   CSAPS({X1,...,Xm},Y, ... )  provides a cubic smoothing spline to data
%         %   values Y on the m-dimensional rectangular grid specified by the  m
%         %   vectors X1, ..., Xm, and these may be of different lengths. Now,
%         %   Y has size [d,length(X1), ..., length(Xm)], with  d  the size of a
%         %   data value.
%         %   If Y is only of size [length(X1), ..., length(Xm)], i.e., the apparent
%         %   d is [], then  d  is taken to be [1], i.e., the function is scalar-valued.
%         %   As to the optional arguments,  P , XX , W , if present, they must be as
%         %   follows:
%         %
%         %   P must be a cell-array with  m  entries, or else an m-vector, except that
%         %   it may also be a scalar or empty, in which case it is converted to an
%         %   m-cell-array with all entries equal to the given P. The optional second
%         %   output argument will always be an m-cell-array.
%         %
%         %    XX  can either be a matrix with  m  rows, and then each of its columns
%         %   is taken as a point in m-space at which the smoothing spline is to be
%         %   evaluated; or else, XX must be a cell-array {XX1, ..., XXm} specifying
%         %   the m-dimensional grid at which to evaluate the smoothing spline.
%         %   With such an XX present, the values of the smoothing spline at the points
%         %   specified by XX are returned. If there is no XX or else XX is empty,
%         %   the ppform of the smoothing spline is returned instead.
%         %
%         %    W  must be a cell array of length  m , with each W{i} either a vector
%         %   of the same length as Xi, or else empty, and in that case the default
%         %   value, ones(1,length(Xi)), is used for W{i}.
%         %
%         %   Example:
%         %      x = {linspace(-2,3,51),linspace(-3,3,61)};
%         %      [xx,yy] = ndgrid(x{1},x{2}); y = peaks(xx,yy);
%         %      noisy = y+(rand(size(y))-.5);
%         %      [smooth,p] = csaps(x,noisy,[],x);
%         %      surf(x{1},x{2},smooth.')
%         %   adds uniform noise from the interval [-.5,.5] to the values of MATLAB's
%         %   PEAKS function on a 51-by-61 uniform grid, then obtains smoothed values
%         %   from CSAPS along with the smoothing parameters chosen by CSAPS, and plots
%         %   these smoothed values. -- Notice the use of NDGRID and the need to use the
%         %   transpose of the array  smooth  in the SURF command.
%         %   If the resulting surface does not strike you as smooth enough, try a
%         %   slightly smaller P than the one, .9998889, used, for each variable,
%         %   by CSAPS in this case:
%         %      smoother = csaps(x,noisy,.996,x);
%         %      figure, surf(x{1},x{2},smoother.')
%         %
%         %   See also SPAPS, CSAPSDEM, TPAPS.
%         
%         %   Copyright 1987-2010 The MathWorks, Inc.
%         
%         if nargin<3||isempty(p), p = -1; end
%         if nargin<4, xx = []; end
%         if nargin<5, w = []; end
%         
%         if iscell(x)     % we are to handle gridded data
%             
%             m = length(x);
%             sizey = size(y);
%             if length(sizey)<m
%                 error(message('SPLINES:CSAPS:toofewdims')), end
%             
%             if length(sizey)==m,  % grid values of a scalar-valued function
%                 if issparse(y), y = full(y); end
%                 sizey = [1 sizey];
%             end
%             
%             sizeval = sizey(1:end-m); sizey = [prod(sizeval), sizey(end-m+(1:m))];
%             y = reshape(y, sizey);
%             
%             if ~iscell(p)  % because of the possibility of weighted roughness measures
%                 % must have P be a cell array in the multivariate case.
%                 if length(p)~=m, p = repmat(p(1),1,m); end
%                 p = num2cell(p);
%             end
%             if isempty(w), w = cell(1,m); end
%             
%             v = y; sizev = sizey;
%             for i=m:-1:1   % carry out coordinatewise smoothing
%                 [cs,p{i}] = csaps1(x{i}, reshape(v,prod(sizev(1:m)),sizev(m+1)), ...
%                     p{i}, [], w{i});
%                 [breaks{i},v,l,k] = ppbrk(cs);
%                 sizev(m+1) = l*k; v = reshape(v,sizev);
%                 if m>1
%                     v = permute(v,[1,m+1,2:m]); sizev(2:m+1) = sizev([m+1,2:m]);
%                 end
%             end
%             % At this point, V contains the tensor-product pp coefficients;
%             % It remains to make up the formal description:
%             output = ppmak(breaks, v);
%             if length(sizeval)>1, output = fnchg(output,'dz',sizeval); end
%             if ~isempty(xx)
%                 output = fnval(output,xx);
%             end
%             
%         else             % we have univariate data
%             
%             [output,p] = csaps1(x,y,p,xx,w);
%             
%         end
%         
%         function [output,p] = csaps1(x,y,p,xx,w)
%             %CSAPS1 univariate cubic smoothing spline
%             
%             n=length(x); if isempty(w), w = ones(1,n); end
%             [xi,yi,sizeval,w,origint,p] = chckxywp(x,y,2,w,p);
%             n = size(xi,1); yd = size(yi,2); dd = ones(1,yd);
%             
%             dx = diff(xi); divdif = diff(yi)./dx(:,dd);
%             if n==2 % the smoothing spline is the straight line interpolant
%                 pp=ppmak(xi.',[divdif.' yi(1,:).'],yd); p = 1;
%             else % set up the linear system for solving for the 2nd derivatives at  xi .
%                 % this is taken from (XIV.6)ff of the `Practical Guide to Splines'
%                 % with the diagonal matrix D^2 there equal to diag(1/w) here.
%                 % Make use of sparsity of the system.
%                 
%                 dxol = dx;
%                 if length(p)>1
%                     lam = p(2:end).'; p = p(1);
%                     dxol = dx./lam;
%                 end
%                 
%                 R = spdiags([dxol(2:n-1), 2*(dxol(2:n-1)+dxol(1:n-2)), dxol(1:n-2)],...
%                     -1:1, n-2,n-2);
%                 odx=1./dx;
%                 Qt = spdiags([odx(1:n-2), -(odx(2:n-1)+odx(1:n-2)), odx(2:n-1)], ...
%                     0:2, n-2,n);
%                 % solve for the 2nd derivatives
%                 W = spdiags(1./w(:),0,n,n);
%                 Qtw = Qt*spdiags(1./sqrt(w(:)),0,n,n);
%                 if p<0 % we are to determine an appropriate P
%                     QtWQ = Qtw*Qtw.'; p = 1/(1+trace(R)/(6*trace(QtWQ)));
%                     % note that the resulting  p  behaves like
%                     %   1/(1 + w_unit*x_unit^3/lambda_unit)
%                     % as a function of the various units chosen
%                     u=((6*(1-p))*QtWQ+p*R)\diff(divdif);
%                 else
%                     u=((6*(1-p))*(Qtw*Qtw.')+p*R)\diff(divdif);
%                 end
%                 % ... and convert to pp form
%                 % Qt.'*u=diff([0;diff([0;u;0])./dx;0])
%                 yi = yi - ...
%                     (6*(1-p))*W*diff([zeros(1,yd)
%                     diff([zeros(1,yd);u;zeros(1,yd)])./dx(:,dd)
%                     zeros(1,yd)]);
%                 c3 = [zeros(1,yd);p*u;zeros(1,yd)];
%                 c2=diff(yi)./dx(:,dd)-dxol(:,dd).*(2*c3(1:n-1,:)+c3(2:n,:));
%                 if exist('lam','var')
%                     dxtl = dx.*lam;
%                     pp=ppmak(xi.',...
%                         reshape([(diff(c3)./dxtl(:,dd)).',3*(c3(1:n-1,:)./lam(:,dd)).', ...
%                         c2.',yi(1:n-1,:).'], (n-1)*yd,4),yd);
%                 else
%                     pp=ppmak(xi.',...
%                         reshape([(diff(c3)./dx(:,dd)).',3*c3(1:n-1,:).',c2.',yi(1:n-1,:).'],...
%                         (n-1)*yd,4),yd);
%                 end
%             end
%             
%             if ~isempty(origint), pp = fnchg(pp,'int',origint); end
%             if length(sizeval)>1, pp = fnchg(pp,'dz',sizeval); end
%             
%             if isempty(xx)
%                 output = pp;
%             else
%                 output = fnval(pp,xx);
%             end
%         end
%     end
end
