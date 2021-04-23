clear
[data,infoBV,marker,auxind_dur_ch] = fopen_EEG('C:\Users\laura\AnalysesNIRS\MART0m\Multimodal\BB043\BB043_0m_Mart_AUX.dat');

fs=7.8125;
[~,q] = rat(fs/500,0.0001);
for c=1:2
    newdata(:,c)=downsample(data(:,c),q);
    newdata(:,c)=zscore(newdata(:,c));
    
    % subplot(2,1,c)
    % histogram(newdata1(:,c));hold on;
    % xline([-3 -1 1 3])
    % xlim([-5 5])
end



threshold=[.35 .7;... %moving STD threshold, respectively for SAT and RESP, in standardized units    .34 .67
    .4   1.4]; %moving AVG threshold, respectively for SAT and RESP, in standardized units    .4   1.0
adjwindow=[1.5 3.5] ; %moving window length before and after X in seconds
%( *2: if you want a 3 sec window = please write 1.5)
minlength=3; %minimum length of intervals (good or bad) in data points

question1=[0 0];
while ~question1
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
            if mstd(x,c)>=threshold(1, c) || abs(mavg(x,c))>=threshold(2, c)
                MA(x,c)=true;
            else
                MA(x,c)=false;
                %mavg(x,c)=mean(newdata(tmpwindow(1):tmpwindow(2),c))  ;
            end
            if x>1 && diff([MA(x-1,c) MA(x,c)])
                xma{c}=[xma{c} x];
                if diff(xma{c}(end-1:end))<minlength
                    xma{c}(end-1:end)=[];
                end
            end
            
        end
    end
    
    
    for c=1:2
        if ~question1(c)
            figure
            subplot(3,1,1);
            area(1:length(newdata(:,c)),MA(:,c)*10,0,'LineStyle','none','FaceColor',[.7 .7 .7]);hold on %box
            plot(newdata(:,c));
            ylim([-5 5])
            subplot(3,1,2);
            area(1:length(newdata(:,c)),MA(:,c)*10,0,'LineStyle','none','FaceColor',[.7 .7 .7]);hold on %box
            plot(mstd(:,c)); ylim([0 1.2]); yline(threshold(1,c))
            subplot(3,1,3);
            area(1:length(newdata(:,c)),MA(:,c)*10,0,'LineStyle','none','FaceColor',[.7 .7 .7]);hold on %box
            plot(mavg(:,c)); ylim([-1.5 1.5]); yline(threshold(2,c))
            
            minn=0:2000:length(newdata(:,c))-2000;
            maxx=minn+2000;
            for xx=1:length(minn)
                for a=1:3
                    subplot(3,1,a); xlim([minn(xx) maxx(xx)]);
                end
                pause(5)
            end
            question1(c)=input('Are you satisfied with the identification of bad intervals for the current aux?  ');
            if ~question1(c)
                threshold(1,c)=input(['New threshold for moving STD? (previous = ' num2str(threshold(1,c)) ').Enter here: ']);
                threshold(2,c)=input(['New threshold for moving AVG? (previous = ' num2str(threshold(2,c)) ').Enter here: ']);
                adjwindow(c)=input(['New duration for moving window (in sec)? (previous = ' num2str(adjwindow(c)) ').Enter here: ']);
            end
            close
        end
        
    end
end
int=[12 1];
pp=[.02 .01];
for c=1:2 %segmentation and spline
    for m=1:length(xma{c})
        if m==length(xma{c})
            dsegment{m,c}= newdata(xma{c}(m):length(newdata(:,c)),c);
        else
            dsegment{m,c}= newdata(xma{c}(m):(xma{c}(m+1)-1),c);
            
        end
        
        if MA(xma{c}(m),c)
            idsegment{m,c}='bad';
            tmpsegment=csaps((1:length(dsegment{m,c}))*int(c),dsegment{m,c}',pp(c),(1:length(dsegment{m,c}))*int(c));
            newsegment{m,c}=dsegment{m,c}-(tmpsegment');
        else
            idsegment{m,c}='good';
            newsegment{m,c}=dsegment{m,c};
        end
    end
end

figure
for c=1:2
    subplot(2,1,c)
    corrdata(:,c)=vertcat(newsegment{1:length(xma{c}),c});
    plot(newdata(:,c));hold on; plot(corrdata(:,c));
end
minn=[0:1000:15000];
maxx=minn+1000;
for xx=1:length(minn)
    for a=1:2
        subplot(2,1,a); xlim([minn(xx) maxx(xx)]);
    end
    pause(3)
end



%
% figure;
% subplot(3,2,1);plot(newdata(:,1));
% subplot(3,2,2);plot(newdata(:,2));
% subplot(3,2,3);plot(mstd(:,1)); ylim([0 1])
% subplot(3,2,4);plot(mstd(:,2));ylim([0 1])
% subplot(3,2,5);plot(mavg(:,1)); ylim([-3 3])
% subplot(3,2,6);plot(mavg(:,2));ylim([-3 3])
%
% yline([0.3:.1:.5])
%
% figure
% for c=1:2
%    subplot(2,1,c)
% plot(newdata(:,c)); hold on;
% yline([-3 -1 1 3])
% yline([0 median(newdata1(:,c))], 'linewidth',2)
% %ylim([-4 4])
% end
% figure
% for c=1:2
%    subplot(2,1,c)
% plot(newdata1(:,c)); hold on;
% yline([-3 -1 1 3])
% yline([0 median(newdata1(:,c))], 'linewidth',2)
% ylim([-4 4])
% end
%
%
%
%     for c=1:2
%         tmpbins=std(newdata(tmpwindow(1):tmpwindow(2),c)     %prctile(newdata(tmpwindow(1):tmpwindow(2),c),[2 98])
%         if newdata(x,c) < tmpbins(1) || newdata(x,c) > tmpbins(2)
%             correcteddata(x,c)=median(newdata(tmpwindow(1):tmpwindow(2),c));
%         else
%             correcteddata(x,c)=newdata(x,c);
%         end
%         clear tmpbins
%     end
%     clear tmpwindow
% end
%
% figure
% subplot(2,1,1)
% plot(newdata(:,1));hold on; plot(correcteddata(:,1)); %xlim([0 1000])
% subplot(2,1,2)
% plot(newdata(:,2));hold on; plot(correcteddata(:,2)); %xlim([0 1000])
