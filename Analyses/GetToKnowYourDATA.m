%%%%%%%%%%%%%%GetToKnowYourDATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
disp('Computing GetToKnowYourDATA')

datapath = 'C:\data\CINC\CINC4m\test_nov21\';
load ([datapath 'workspace.mat'])
%load ([datapath 'workspacemat.mat'])

fileorderconnectogram = {'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Mixte.txt',...
                         'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Region.txt',...
                         'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Fonction.txt',...
                         'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Aire.txt'};

descrstatsmode = 1;
graphmode = 1;
channelmode = 1;
importROI = 0;
calculateROI = 0;

%%
%%Stats descriptives%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Computing Descriptive Statistics')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if descrstatsmode
    %Channels%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if importROI == 0 & channelmode == 1
        %données de FC séparée par groupe
        datachG1 = datach(idG1,:);
        datachG2 = datach(idG2,:);
        
        %normalité et valeurs extrêmes%
        meanG1ch = nanmean(datachG1);
        meanG2ch = nanmean(datachG1);
        stdG1ch = nanstd(datachG1);
        stdG2ch = nanstd(datachG1);
        sknG1ch = skewness(datachG1);
        sknG2ch = skewness(datachG1);
        krtG1ch = kurtosis(datachG1);
        krtG2ch = kurtosis(datachG1);

        tblgrmeanch = [array2table([meanG1ch],'VariableNames',labelch); array2table([meanG2ch],'VariableNames',labelch)];
        tblgrmeanch.Properties.RowNames = {'G1','G2'};
        tbldescrch = array2table([meanG1ch; meanG2ch; stdG1ch; stdG2ch; sknG1ch; sknG2ch; krtG1ch; krtG2ch],'VariableNames',labelch,'RowNames',{'meanG1','meanG2','stdG1','stdG2','sknG1','sknG2','krtG1','krtG2'});

        tf = tbldescrch{{'sknG1','sknG2','krtG1','krtG2'},:} >2 | tbldescrch{{'sknG1','sknG2','krtG1','krtG2'},:} <-2;
        n_abnormal = sum(any(tf));
        p_abnormal = n_abnormal/numel(tbldescrch(1,:))*100;
        fprintf('%.1f percent of the channels variables are not respecting the normal distribution\n',p_abnormal);

        nanzscore = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));
        zchALL = nanzscore(datach);
        tf = zchALL > 3.29 | zchALL <-3.29;
        n_extreme = sum(any(tf));
        p_extreme = n_extreme/numel(zchALL)*100;
        fprintf('%.1f percent of the channels variables have extreme values\n',p_extreme);

        %%NAN Values%%%%%%%%%%
        npart_nan = sum(isnan(datach),2); %nb de NAN channels par participant
        percentnan = sum(npart_nan)/numel(datach)*100; %pourcentage NAN channels pour tous participants
        fprintf('%.1f percent of the data is a missing value\n',percentnan)

        [sortpartnan,indexpartnan] = sort(npart_nan,'descend'); %ordonner les participants + au - NAN
        mostnanpart = sortpartnan > 50/100*size(datach,2); %50% ou plus de données manquantes
        pbnanpart = part(indexpartnan(mostnanpart)); %identifier les participants
        fprintf('%s has more than 50 percent missing values\n',pbnanpart)

        clear npart_nan percentnan sortpartnan indexpartnan mostnanpart pbnanpart


        x = 1;
        for p = 1:size(MATall,3) %pour chaque participant
            tfnanch = mean(isnan(MATall(:,:,p)),1) == 1; %trouver les canaux complètement NAN
            idxnanch = find(tfnanch); %trouver leur position dans la matrice
            nanch{x,1} = idxnanch; %liste des canaux manquants pour chaque participant
            x = x + 1;
        end

        nb = 0;
        for c = 1:size(MATall,1) %pour chaque canal
            for n = 1 : size(nanch,1) %pour chaque participant
                tf = sum(nanch{n,1} == c); %déterminer si le canal est manquant
                nb = nb + tf; % nombre de fois ou le canal est manquant
            end
            nanfreq{c,1} = nb; %liste du nb de fois ou chaque canal est manquant
            nb = 0;

            if calculateROI == 1
                R = 4; %Aroi
                for n = 1:size(roi{1,R},2) %pour chaque roi
                    tf = roi{1,R}{2,n} == c; %trouver si le canal fait partie de la roi
                    if any(tf) %si oui
                        nanfreq{c,2} = roi{1,R}(1,n); %ajouter le label
                    end
                end
            end
        end

        if calculateROI
            for n = 1:size(roi{1,R},2) % pour chaque roi
                tf = strcmp(roi{1,R}(1,n),[nanfreq{:,2}]); %trouver chaque canal qui provient de chaque région
                idx = find(tf); %trouver la position de chaque canal
                nantot(n,1) = sum([nanfreq{idx,1}]); %faire la somme des canaux NAN de chaque région
            end
            nantot = num2cell(nantot);
            nantot(:,2) = roi{1,R}(1,:);
        end

        [sortchnan,indexchnan] = sort(cell2mat(nanfreq(:,1)),'descend');
        mostnanch = sortchnan > 30/100*size(MATall,3); %30% ou plus de données manquantes
        pbnanch = Listch(indexchnan(mostnanch));
        if isempty(pbnanch) == 0
           fprintf('Channel %.0f has more than 30 percent missing values\n',pbnanch)
        else
        end

        figure
        X = categorical(1:size(MATall,1));
        bar(X,cell2mat(nanfreq(:,1)))
        yline(mean(cell2mat(nanfreq(:,1)),1))
        ylabel('Total number of missing')
        title('Total number of missing channels')
        savefig([savepath date '_NANCh']);
        exportgraphics(gcf,[savepath 'NANCh.png'])

        if calculateROI
            figure
            X = categorical(roi{1,R}(1,:));
            bar(X,cell2mat(nantot(:,1)))
            yline(mean(cell2mat(nantot(:,1)),1))
            ylabel('Total number of channels')
            title('Total number of channels with missing data in each region')
            savefig([savepath date '_NANChBYROI']);
            exportgraphics(gcf,[savepath 'NANChBYROI.png'])

            clear nantot
        end

        clear R X x p c n idx tf nanch tfnanch idxnanch nb sortchnan indexchnan mostnanch pbnanch nanfreq
    end

    %%%%%Roi%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ROImode
        if importROI == 1
            %normalité et valeurs extrêmes%
            meanG1roi = nanmean(dataroi(idG1,:));
            meanG2roi = nanmean(dataroi(idG2,:));
            stdG1roi = nanstd(dataroi(idG1,:));
            stdG2roi = nanstd(dataroi(idG2,:));
            sknG1roi = skewness(dataroi(idG1,:));
            sknG2roi = skewness(dataroi(idG2,:));
            krtG1roi = kurtosis(dataroi(idG1,:));
            krtG2roi = kurtosis(dataroi(idG2,:));

            tblgrmeanroi = [array2table([meanG1roi],'VariableNames',labelroi); array2table([meanG2roi],'VariableNames',labelroi)];
            tblgrmeanroi.Properties.RowNames = {'G1','G2'};
            tbldescrroi = array2table([meanG1roi; meanG2roi; stdG1roi; stdG2roi; sknG1roi; sknG2roi; krtG1roi; krtG2roi],'VariableNames',labelroi,'RowNames',{'meanG1','meanG2','stdG1','stdG2','sknG1','sknG2','krtG1','krtG2'});

            tf = tbldescrroi{{'sknG1','sknG2','krtG1','krtG2'},:} >2 | tbldescrroi{{'sknG1','sknG2','krtG1','krtG2'},:} <-2;
            n_abnormal = sum(any(tf));
            p_abnormal = n_abnormal/numel(tbldescrroi(1,:))*100;
            fprintf('%.1f percent of the ROIs variables are not respecting the normal distribution\n',p_abnormal);

            nanzscore = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));
            zdataroi = nanzscore(dataroi);
            tf = zdataroi > 3.29 | zdataroi <-3.29;
            n_extreme = sum(any(tf));
            p_extreme = n_extreme/numel(zdataroi)*100;
            fprintf('%.1f percent of the ROIs variables have extreme values\n',p_extreme);

            %%NAN Values%%%%%%%%%%
            npart_nan = sum(isnan(dataroi),2); %nb de NAN ROIs par participant
            percentnan = sum(npart_nan)/numel(dataroi)*100; %pourcentage NAN ROIs pour tous participants
            fprintf('%.1f percent of the data is a missing value\n',percentnan)

            [sortpartnan,indexpartnan] = sort(npart_nan,'descend'); %ordonner les participants + au - NAN
            mostnanpart = sortpartnan > 40/100*size(dataroi,2); %40% ou plus de données manquantes
            pbnanpart = part(indexpartnan(mostnanpart)); %identifier les participants
            fprintf('%s has more than 40 percent missing values\n',pbnanpart)

            clear npart_nan percentnan sortpartnan indexpartnan mostnanpart pbnanpart n_abnormal p_abnormal n_extreme p_extreme tf

            x = 1;
            for p = 1:size(MATall,3) %pour chaque participant
                tfnanroi = mean(isnan(MATall(:,:,p)),1) >0.9; %trouver les ROIs complètement NAN. AJUSTÉ CAR DIAG=0 MEME SI NAN
                idxnanroi = find(tfnanroi); %trouver leur position dans la matrice
                nanroi{x,1} = idxnanroi; %liste des ROIs manquants pour chaque participant
                x = x + 1;
            end

            nb = 0;
            for r = 1:size(MATall,1) %pour chaque ROI
                for n = 1 : size(nanroi,1) %pour chaque participant
                    tf = sum(nanroi{n,1} == r); %déterminer si le ROI est manquant
                    nb = nb + tf; % nombre de fois ou le ROI est manquant
                end

                nanfreq{r,1} = nb; %liste du nb de fois ou chaque ROI est manquant
                nb = 0;
            end

            pbnanroi = [];
            [sortroinan,indexroinan] = sort(cell2mat(nanfreq(:,1)),'descend');
            mostnanroi = sortroinan > 10/100*size(MATall,3); %10% ou plus de données manquantes
            pbnanroi = Listroi(indexroinan(mostnanroi));
            if isempty(pbnanroi) == 0
                fprintf('ROI %s has more than 10 percent missing values\n',pbnanroi)
            else
            end

            figure
            X = categorical(Listroi,Listroi);
            bar(X,cell2mat(nanfreq(:,1)))
            yline(mean(cell2mat(nanfreq(:,1)),1))
            ylabel('Total number of missing')
            title('Total number of missing ROIs')
            savefig([savepath date '_NANROI']);
            exportgraphics(gcf,[savepath 'NANROI.png'])

            clear x X p r n idx tf nanroi tfnanroi idxnanroi nb sortroinan indexroinan mostnanroi pbnanroi nanfreq
        end  

        if calculateROI ==1 & importROI == 0
            %normalité et valeurs extrêmes%
            ALLroi = [];
            for R = 1:4
                ALLroi = [ALLroi, dataroiALL{1,R}]; %%Fait pour tous les ROIS calculés ensemble%%%%%%
            end

            meanG1roi = nanmean(ALLroi(idG1,:));
            meanG2roi = nanmean(ALLroi(idG2,:));
            stdG1roi = nanstd(ALLroi(idG1,:));
            stdG2roi = nanstd(ALLroi(idG2,:));
            sknG1roi = skewness(ALLroi(idG1,:));
            sknG2roi = skewness(ALLroi(idG2,:));
            krtG1roi = kurtosis(ALLroi(idG1,:));
            krtG2roi = kurtosis(ALLroi(idG2,:));

            tblgrmeanroi = [array2table([meanG1roi],'VariableNames',[labelroiALL{1,1:4}]); array2table([meanG2roi],'VariableNames',[labelroiALL{1,1:4}])];
            tblgrmeanroi.Properties.RowNames = {'G1','G2'};
            tbldescrroi = array2table([meanG1roi; meanG2roi; stdG1roi; stdG2roi; sknG1roi; sknG2roi; krtG1roi; krtG2roi],'VariableNames',[labelroiALL{1,1:4}],'RowNames',{'meanG1','meanG2','stdG1','stdG2','sknG1','sknG2','krtG1','krtG2'});

            tf = tbldescrroi{{'sknG1','sknG2','krtG1','krtG2'},:} >2 | tbldescrroi{{'sknG1','sknG2','krtG1','krtG2'},:} <-2;
            n_abnormal = sum(any(tf));
            p_abnormal = n_abnormal/numel(tbldescrroi(1,:))*100;
            fprintf('%.1f percent of the ROIs variables are not respecting the normal distribution\n',p_abnormal);

            zALLroi = nanzscore(ALLroi);
            tf = zALLroi > 3.29 | zALLroi <-3.29;
            n_extreme = sum(any(tf));
            p_extreme = n_extreme/numel(zALLroi)*100;
            fprintf('%.1f percent of the ROIs variables have extreme values\n',p_extreme);

%             meanroiG1 = nanmean(dataroiALL{1,1}(idG1,:));
%             meanroiG2 = nanmean(dataroiALL{1,1}(idG2,:));
%             meanRroiG1 = nanmean(dataroiALL{1,2}(idG1,:));
%             meanRroiG2 = nanmean(dataroiALL{1,2}(idG2,:));
%             meanFroiG1 = nanmean(dataroiALL{1,3}(idG1,:));
%             meanFroiG2 = nanmean(dataroiALL{1,3}(idG2,:));
%             meanAroiG1 = nanmean(dataroiALL{1,4}(idG1,:));
%             meanAroiG2 = nanmean(dataroiALL{1,4}(idG2,:));

            clear X tf R ALLroi n_abnormal n_extreme p_abnormal p_extreme stdG1ch stdG2ch sknG1ch sknG2ch krtG1ch krtG2ch meanG1roi meanG2roi stdG1roi stdG2roi sknG1roi sknG2roi krtG1roi krtG2roi zchALL zroiALLmean zALLroi
        end
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%GRAPHIQUES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if graphmode
    disp('Computing Graphics')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%PAR CANAL%%%%%
    if channelmode
        %données de FC séparée par groupe
        datachG1 = datach(idG1,:);
        datachG2 = datach(idG2,:);

        %moyenne et std de la FC de chaque participant pour tous les canaux%
        meanchG1part = mean(datachG1,2,'omitnan');
        meanchG2part = mean(datachG2,2,'omitnan');
        MmeanchG1part = mean(meanchG1part);
        MmeanchG2part = mean(meanchG2part,'omitnan');
        stdchG1part = std(datachG1,0,2,'omitnan');
        stdchG2part = std(datachG2,0,2,'omitnan');
        MstdchG1part = mean(stdchG1part);
        MstdchG2part = mean(stdchG2part,'omitnan');

        %moyenne et std de la FC de chaque canal pour tous les participants%
        meanchG1ch = mean(datachG1,'omitnan');
        meanchG2ch = mean(datachG2,'omitnan');
        MmeanchG1ch = mean(meanchG1ch);
        MmeanchG2ch = mean(meanchG2ch);
        stdchG1ch = std(datachG1,0,1,'omitnan');
        stdchG2ch = std(datachG2,0,1,'omitnan');
        MstdchG1ch = mean(stdchG1ch);
        MstdchG2ch = mean(stdchG2ch);

        %moyenne et std de la FC pour tous les canaux et participants%
        MmeanchG1 = mean(datachG1,'all','omitnan');
        MmeanchG2 = mean(datachG2,'all','omitnan');
        MstdchG1 = std(datachG1,0,'all','omitnan');
        MstdchG2 = std(datachG2,0,'all','omitnan');

        figure
        hold on
        histg1 = histogram(meanchG1part);
        histg1.Normalization = 'pdf';
        histg1.BinWidth = 0.05;
        histg1.DisplayName = 'FC G1';
        xmin = floor((min(meanchG1part,[],'all'))*10)/10;
        xmax = ceil(max(meanchG1part,[],'all')*10)/10;
        x = xmin:0.05:xmax;
        pd = fitdist(meanchG1part,'Normal');
        y = pdf(pd,x);
        pg1 = plot(x,y,'LineWidth',1);
        pg1.Color = [0 0.4470 0.7410];
        histg2 = histogram(meanchG2part);
        histg2.Normalization = 'pdf';
        histg2.BinWidth = 0.05;
        histg2.DisplayName = 'FC G2';
        xmin = floor((min(meanchG2part,[],'all'))*10)/10;
        xmax = ceil(max(meanchG2part,[],'all')*10)/10;
        x = xmin:0.05:xmax;
        pd = fitdist(meanchG2part,'Normal');
        y = pdf(pd,x);
        pg2 = plot(x,y,'LineWidth',1);
        pg2.Color = [0.9290 0.6940 0.1250];
        dim = [.15 .6 .3 .3];
        str = sprintf('Mean G1 = %.3f\n Std G1 = %.3f\n Mean G2 = %.3f\n Std G2 = %.3f',MmeanchG1part,MstdchG1part,MmeanchG2part,MstdchG2part);
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
        legend;
        xlabel('Functional Connectivity Value');
        ylabel('Proportion of the participants (%)');
        title('Channels Histogram of the participants FC values')
        savefig([savepath date '_HistPartch']);
        exportgraphics(gcf,[savepath 'HistPartch.png'])    

        clear histg1 histg2 pg1 pg2 xmin xmax x pd y dim str xlabel ylabel

        figure
        hold on
        histg1 = histogram(meanchG1ch);
        histg1.Normalization = 'pdf';
        histg1.BinWidth = 0.05;
        histg1.DisplayName = 'FC G1';
        xmin = floor((min(meanchG1ch,[],'all'))*10)/10;
        xmax = ceil(max(meanchG1ch,[],'all')*10)/10;
        x = xmin:0.05:xmax;
        pd = fitdist(meanchG1ch','Normal');
        y = pdf(pd,x);
        pg1 = plot(x,y,'LineWidth',1);
        pg1.Color = [0 0.4470 0.7410];
        histg2 = histogram(meanchG2ch);
        histg2.Normalization = 'pdf';
        histg2.BinWidth = 0.05;
        histg2.DisplayName = 'FC G2';
        xmin = floor((min(meanchG2ch,[],'all'))*10)/10;
        xmax = ceil(max(meanchG2ch,[],'all')*10)/10;
        x = xmin:0.05:xmax;
        pd = fitdist(meanchG2ch','Normal');
        y = pdf(pd,x);
        pg2 = plot(x,y,'LineWidth',1);
        pg2.Color = [0.9290 0.6940 0.1250];
        dim = [.15 .6 .3 .3];
        str = sprintf('Mean G1 = %.3f\n Std G1 = %.3f\n Mean G2 = %.3f\n Std G2 = %.3f',MmeanchG1ch,MstdchG1ch,MmeanchG2ch,MstdchG2ch);
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
        legend;
        xlabel('Functional Connectivity Value');
        ylabel('Proportion of the Channels (%)');
        title('Channels Histogram of the Channels FC values')
        savefig([savepath date '_HistCHch']);
        exportgraphics(gcf,[savepath 'HistCHch.png'])

        clear histg1 histg2 pg1 pg2 xmin xmax x pd y dim str xlabel ylabel

        figure
        hold on
        histg1 = histogram(datachG1);
        histg1.Normalization = 'pdf';
        histg1.BinWidth = 0.05;
        histg1.DisplayName = 'FC G1';
        limits = histg1.BinLimits;
        x = limits(1,1):histg1.BinWidth:limits(1,2);
        mu = MmeanchG1ch;
        sigma = MstdchG1ch;
        f = exp(-(x-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
        pg1 = plot(x,f,'LineWidth',1);
        pg1.Color = [0 0.4470 0.7410];
        histg2 = histogram(datachG2);
        histg2.Normalization = 'pdf';
        histg2.BinWidth = 0.05;
        histg2.DisplayName = 'FC G2';
        limits = histg2.BinLimits;
        x = limits(1,1):histg2.BinWidth:limits(1,2);
        mu = MmeanchG2ch;
        sigma = MstdchG2ch;
        f = exp(-(x-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
        pg2 = plot(x,f,'LineWidth',1);
        pg2.Color = [0.9290 0.6940 0.1250];
        legend;
        xlabel('Functional Connectivity Value');
        ylabel('Proportion of the FC (%)');
        title('Channels Histogram of the FC values')
        dim = [.15 .6 .3 .3];
        str = sprintf('Mean G1 = %.3f\n Std G1 = %.3f\n Mean G2 = %.3f\n Std G2 = %.3f',MmeanchG1ch,MstdchG1ch,MmeanchG2ch,MstdchG2ch);
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
        savefig([savepath date '_HistFCch']);
        exportgraphics(gcf,[savepath 'HistFCch.png'])

        clear histg1 histg2 pg1 pg2 x mu sigma f dim str limits MmeanchG1 MmeanchG1part MmeanchG1ch MmeanchG2 MmeanchG2part MmeanchG2ch MstdchG1 MstdchG1part MstdchG1ch MstdchG2 MstdchG2part MstdchG2ch

    end

    %%%%%%%%%ROI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if importROI
        %données de FC séparée par groupe
        dataroiG1 = dataroi(idG1,:);
        dataroiG2 = dataroi(idG2,:);

        %moyenne et std de la FC de chaque participant pour toutes les régions%
        meanroiG1part = mean(dataroiG1,2,'omitnan');
        meanroiG2part = mean(dataroiG2,2,'omitnan');
        MmeanroiG1part = mean(meanroiG1part);
        MmeanroiG2part = mean(meanroiG2part,'omitnan');
        stdroiG1part = std(dataroiG1,0,2,'omitnan');
        stdroiG2part = std(dataroiG2,0,2,'omitnan');
        MstdroiG1part = mean(stdroiG1part);
        MstdroiG2part = mean(stdroiG2part,'omitnan');

        %moyenne et std de la FC de chaque région pour tous les participants%
        meanroiG1roi = mean(dataroiG1,'omitnan');
        meanroiG2roi = mean(dataroiG2,'omitnan');
        MmeanroiG1roi = mean(meanroiG1roi);
        MmeanroiG2roi = mean(meanroiG2roi);
        stdroiG1roi = std(dataroiG1,0,1,'omitnan');
        stdroiG2roi = std(dataroiG2,0,1,'omitnan');
        MstdroiG1roi = mean(stdroiG1roi);
        MstdroiG2roi = mean(stdroiG2roi);

        %moyenne et std de la FC pour toutes les régions et participants%
        MmeanroiG1 = mean(dataroiG1,'all','omitnan');
        MmeanroiG2 = mean(dataroiG2,'all','omitnan');
        MstdroiG1 = std(dataroiG1,0,'all','omitnan');
        MstdroiG2 = std(dataroiG2,0,'all','omitnan');

        %%HISTOGRAM%%%%%%%%%%%

        figure
        hold on
        histg1 = histogram(meanroiG1part);
        histg1.Normalization = 'pdf';
        histg1.BinWidth = 0.05;
        histg1.DisplayName = 'FC G1';
        xmin = floor((min(meanroiG1part,[],'all'))*10)/10;
        xmax = ceil(max(meanroiG1part,[],'all')*10)/10;
        x = xmin:0.05:xmax;
        pd = fitdist(meanroiG1part,'Normal');
        y = pdf(pd,x);
        pg1 = plot(x,y,'LineWidth',1);
        pg1.Color = [0 0.4470 0.7410];
        histg2 = histogram(meanroiG2part);
        histg2.Normalization = 'pdf';
        histg2.BinWidth = 0.05;
        histg2.DisplayName = 'FC G2';
        xmin = floor((min(meanroiG2part,[],'all'))*10)/10;
        xmax = ceil(max(meanroiG2part,[],'all')*10)/10;
        x = xmin:0.05:xmax;
        pd = fitdist(meanroiG2part,'Normal');
        y = pdf(pd,x);
        pg2 = plot(x,y,'LineWidth',1);
        pg2.Color = [0.9290 0.6940 0.1250];
        dim = [.15 .6 .3 .3];
        str = sprintf('Mean G1 = %.3f\n Std G1 = %.3f\n Mean G2 = %.3f\n Std G2 = %.3f',MmeanroiG1part,MstdroiG1part,MmeanroiG2part,MstdroiG2part);
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
        legend;
        xlabel('Functional Connectivity Value');
        ylabel('Proportion of the participants (%)');
        title('Roi Histogram of the participants FC values')
        savefig([savepath date '_HistROIPart']);
        exportgraphics(gcf,[savepath 'HistROIPart.png'])

        clear histg1 histg2 pg1 pg2 xmin xmax x pd y dim str xlabel ylabel

        figure
        hold on
        histg1 = histogram(meanroiG1roi);
        histg1.Normalization = 'pdf';
        histg1.BinWidth = 0.05;
        histg1.DisplayName = 'FC G1';
        xmin = floor((min(meanroiG1roi,[],'all'))*10)/10;
        xmax = ceil(max(meanroiG1roi,[],'all')*10)/10;
        x = xmin:0.05:xmax;
        pd = fitdist(meanroiG1roi','Normal');
        y = pdf(pd,x);
        pg1 = plot(x,y,'LineWidth',1);
        pg1.Color = [0 0.4470 0.7410];
        histg2 = histogram(meanroiG2roi);
        histg2.Normalization = 'pdf';
        histg2.BinWidth = 0.05;
        histg2.DisplayName = 'FC G2';
        xmin = floor((min(meanroiG2roi,[],'all'))*10)/10;
        xmax = ceil(max(meanroiG2roi,[],'all')*10)/10;
        x = xmin:0.05:xmax;
        pd = fitdist(meanroiG2roi','Normal');
        y = pdf(pd,x);
        pg2 = plot(x,y,'LineWidth',1);
        pg2.Color = [0.9290 0.6940 0.1250];
        dim = [.15 .6 .3 .3];
        str = sprintf('Mean G1 = %.3f\n Std G1 = %.3f\n Mean G2 = %.3f\n Std G2 = %.3f',MmeanroiG1roi,MstdroiG1roi,MmeanroiG2roi,MstdroiG2roi);
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
        legend;
        xlabel('Functional Connectivity Value');
        ylabel('Proportion of the ROIs (%)');
        title('Roi Histogram of the ROIs FC values')
        savefig([savepath date '_HistROIRoi']);
        exportgraphics(gcf,[savepath 'HistROIRoi.png'])

        clear histg1 histg2 pg1 pg2 xmin xmax x pd y dim str xlabel ylabel

        figure
        hold on
        histg1 = histogram(dataroiG1);
        histg1.Normalization = 'pdf';
        histg1.BinWidth = 0.05;
        histg1.DisplayName = 'FC G1';
        limits = histg1.BinLimits;
        x = limits(1,1):0.05:limits(1,2);
        mu = MmeanroiG1;
        sigma = MstdroiG1;
        f = exp(-(x-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
        pg1 = plot(x,f,'LineWidth',1);
        pg1.Color = [0 0.4470 0.7410];
        histg2 = histogram(dataroiG2);
        histg2.Normalization = 'pdf';
        histg2.BinWidth = 0.05;
        histg2.DisplayName = 'FC G2';
        limits = histg2.BinLimits;
        x = limits(1,1):0.05:limits(1,2);
        mu = MmeanroiG2;
        sigma = MstdroiG2;
        f = exp(-(x-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
        pg2 = plot(x,f,'LineWidth',1);
        pg2.Color = [0.9290 0.6940 0.1250];
        legend;
        xlabel('Functional Connectivity Value');
        ylabel('Proportion of the FC (%)');
        title('Channels Histogram of the FC values')
        dim = [.15 .6 .3 .3];
        str = sprintf('Mean G1 = %.3f\n Std G1 = %.3f\n Mean G2 = %.3f\n Std G2 = %.3f',MmeanroiG1,MstdroiG1,MmeanroiG2,MstdroiG2);
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
        savefig([savepath date '_HistROIFC']);
        exportgraphics(gcf,[savepath 'HistROIFC.png'])

        clear histg1 histg2 pg1 pg2 x mu sigma f dim str limits

        %%%HISTFIT%%%%%%%%%%%

        figure
        hold on
        histg1 = histfit(meanroiG1part);
        alpha(histg1,.5)
        histg1(2).Color = [0 0.4470 0.7410];
        histg2 = histfit(meanroiG2part);
        alpha(histg2,.5)
        histg2(2).Color = [0.9290 0.6940 0.1250];
        dim = [.15 .6 .3 .3];
        str = sprintf('Mean G1 = %.3f\n Std G1 = %.3f\n Mean G2 = %.3f\n Std G2 = %.3f',MmeanroiG1part,MstdroiG1part,MmeanroiG2part,MstdroiG2part);
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
        legend('FC G1','gaussian G1','FC G2','gaussian G2');
        xlabel('Functional Connectivity Value');
        ylabel('Proportion of the participants (%)');
        savefig([savepath date '_HistFitROIPart']);
        exportgraphics(gcf,[savepath 'HistFitROIPart.png'])

        clear histg1 histg2 dim str xlabel ylabel

        figure
        hold on
        histg1 = histfit(meanroiG1roi);
        alpha(histg1,.5)
        histg1(2).Color = [0 0.4470 0.7410];
        histg2 = histfit(meanroiG2roi);
        alpha(histg2,.5)
        histg2(2).Color = [0.9290 0.6940 0.1250];
        dim = [.15 .6 .3 .3];
        str = sprintf('Mean G1 = %.3f\n Std G1 = %.3f\n Mean G2 = %.3f\n Std G2 = %.3f',MmeanroiG1roi,MstdroiG1roi,MmeanroiG2roi,MstdroiG2roi);
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
        legend('FC G1','gaussian G1','FC G2','gaussian G2');
        xlabel('Functional Connectivity Value');
        ylabel('Proportion of the ROIs (%)');
        title('Distribution of the FC values between groups')
        savefig([savepath date '_HistFitROIRoi']);
        exportgraphics(gcf,[savepath 'HistFitROIRoi.png'])

        clear histg1 histg2 dim str xlabel ylabel

        %%HISTFIT ROI BY ROI ONLY IF FEW ROIS%%%
    %     for r = 1:size(dataroi,2)
    %         figure
    %         hold on
    %         histg1 = histfit(dataroi(idG1,r));
    %         alpha(histg1,.5)
    %         histg1(2).Color = [0 0.4470 0.7410];
    %         histg2 = histfit(dataroi(idG2,r));
    %         alpha(histg2,.5)
    %         histg2(2).Color = [0.9290 0.6940 0.1250];
    %         dim = [.15 .6 .3 .3];
    %         meanG1 = meanroiG1roi(:,r);
    %         stdG1 = stdroiG1roi(:,r);
    %         meanG2 = meanroiG2roi(:,r);
    %         stdG2 = stdroiG2roi(:,r);
    %         str = sprintf('Mean G1 = %.3f\n Std G1 = %.3f\n Mean G2 = %.3f\n Std G2 = %.3f',meanG1,stdG1,meanG2,stdG2);
    %         annotation('textbox',dim,'String',str,'FitBoxToText','on');
    %         legend('FC G1','gaussian G1','FC G2','gaussian G2');
    %         xlabel('Functional Connectivity Value');
    %         ylabel('Proportion of the participants (%)');
    %     end

          clear r histg1 histg2 dim str meang1 meang2 stdg1 stdg2 xlabel ylabel


    end

    if calculateROI
        for R = 1:numel(roi)
            ROIname = lgndroi{1,R};
            %données de FC séparée par groupe
            dataroiG1{:,R} = dataroiALL{1,R}(idG1,:);
            dataroiG2{:,R} = dataroiALL{1,R}(idG2,:);

            %moyenne et std de la FC de chaque participant pour tous les canaux%
            meanroiG1p{:,R} = nanmean(dataroiG1{:,R},2);
            meanroiG2p{:,R} = nanmean(dataroiG2{:,R},2);
            stdroiG1p{:,R} = std(dataroiG1{:,R},0,2,'omitnan');
            stdroiG2p{:,R} = std(dataroiG2{:,R},0,2,'omitnan');
            MmeanroiG1p{:,R} = mean(meanroiG1p{:,R});
            MmeanroiG2p{:,R} = mean(meanroiG2p{:,R},'omitnan');
            MstdroiG1p{:,R} = mean(stdroiG1p{:,R});
            MstdroiG2p{:,R} = mean(stdroiG2p{:,R},'omitnan');

            %moyenne de la FC de chaque région pour tous les participants%
            meanroiG1r{:,R} = nanmean(dataroiG1{:,R});
            meanroiG2r{:,R} = nanmean(dataroiG2{:,R});
            stdroiG1r{:,R} = std(dataroiG1{:,R},0,1,'omitnan');
            stdroiG2r{:,R} = std(dataroiG2{:,R},0,1,'omitnan');
            MmeanroiG1r{:,R} = mean(meanroiG1r{:,R});
            MmeanroiG2r{:,R} = mean(meanroiG2r{:,R},'omitnan');
            MstdroiG1r{:,R} = mean(stdroiG1r{:,R});
            MstdroiG2r{:,R} = mean(stdroiG2r{:,R},'omitnan');

            %moyenne de la FC pour tous les canaux et participants%
            MmeanroiG1{:,R} = mean(dataroiG1{:,R},'all','omitnan');
            MmeanroiG2{:,R} = mean(dataroiG2{:,R},'all','omitnan');
            MstdroiG1{:,R} = std(dataroiG1{:,R},0,'all','omitnan');
            MstdroiG2{:,R} = std(dataroiG2{:,R},0,'all','omitnan');
            
            figure
            hold on
            histg1 = histogram(meanroiG1p{:,R});
            histg1.Normalization = 'pdf';
            histg1.BinWidth = 0.05;
            histg1.DisplayName = 'FC G1';
            xmin = floor((min(meanroiG1p{:,R},[],'all'))*10)/10;
            xmax = ceil(max(meanroiG1p{:,R},[],'all')*10)/10;
            x = xmin:0.05:xmax;
            pd = fitdist(meanroiG1p{:,R},'Normal');
            y = pdf(pd,x);
            pg1 = plot(x,y,'LineWidth',1);
            pg1.Color = [0 0.4470 0.7410];
            histg2 = histogram(meanroiG2p{:,R});
            histg2.Normalization = 'pdf';
            histg2.BinWidth = 0.05;
            histg2.DisplayName = 'FC G2';
            xmin = floor((min(meanroiG2p{:,R},[],'all'))*10)/10;
            xmax = ceil(max(meanroiG2p{:,R},[],'all')*10)/10;
            x = xmin:0.05:xmax;
            pd = fitdist(meanroiG2p{:,R},'Normal');
            y = pdf(pd,x);
            pg2 = plot(x,y,'LineWidth',1);
            pg2.Color = [0.9290 0.6940 0.1250];
            dim = [.15 .6 .3 .3];
            str = sprintf('Mean G1 = %.3f\n Std G1 = %.3f\n Mean G2 = %.3f\n Std G2 = %.3f',MmeanroiG1p{:,R},MstdroiG1p{:,R},MmeanroiG2p{:,R},MstdroiG2p{:,R});
            annotation('textbox',dim,'String',str,'FitBoxToText','on');
            legend;
            xlabel('Functional Connectivity Value');
            ylabel('Proportion of the participants (%)');
            str = sprintf('%s Histogram of the participants FC values',ROIname);
            title(str);
            %savefig([savepath date '_HistROIPart']);
            %exportgraphics(gcf,[savepath 'HistROIPart.png'])
            
            clear histg1 histg2 pg1 pg2 xmin xmax x pd y dim str xlabel ylabel
            
            figure
            hold on
            histg1 = histogram(meanroiG1r{:,R});
            histg1.Normalization = 'pdf';
            histg1.BinWidth = 0.05;
            histg1.DisplayName = 'FC G1';
            xmin = floor((min(meanroiG1r{:,R},[],'all'))*10)/10;
            xmax = ceil(max(meanroiG1r{:,R},[],'all')*10)/10;
            x = xmin:0.05:xmax;
            pd = fitdist(meanroiG1r{:,R}','Normal');
            y = pdf(pd,x);
            pg1 = plot(x,y,'LineWidth',1);
            pg1.Color = [0 0.4470 0.7410];
            histg2 = histogram(meanroiG2r{:,R});
            histg2.Normalization = 'pdf';
            histg2.BinWidth = 0.05;
            histg2.DisplayName = 'FC G2';
            xmin = floor((min(meanroiG2r{:,R},[],'all'))*10)/10;
            xmax = ceil(max(meanroiG2r{:,R},[],'all')*10)/10;
            x = xmin:0.05:xmax;
            pd = fitdist(meanroiG2r{:,R}','Normal');
            y = pdf(pd,x);
            pg2 = plot(x,y,'LineWidth',1);
            pg2.Color = [0.9290 0.6940 0.1250];
            dim = [.15 .6 .3 .3];
            str = sprintf('Mean G1 = %.3f\n Std G1 = %.3f\n Mean G2 = %.3f\n Std G2 = %.3f',MmeanroiG1r{:,R},MstdroiG1r{:,R},MmeanroiG2r{:,R},MstdroiG2r{:,R});
            annotation('textbox',dim,'String',str,'FitBoxToText','on');
            legend;
            xlabel('Functional Connectivity Value');
            ylabel('Proportion of the ROIs (%)');
            str = sprintf('%s Histogram of the ROIs FC values',ROIname);
            title(str);
            %savefig([savepath date '_HistROIRoi']);
            %exportgraphics(gcf,[savepath 'HistROIRoi.png'])
            
            clear histg1 histg2 pg1 pg2 xmin xmax x pd y dim str xlabel ylabel MmeanroiG1 MmeanroiG1p MmeanroiG1r MmeanroiG2 MmeanroiG2p MmeanroiG2r MstdroiG1 MstdroiG1p MstdroiG1r MstdroiG2 MstdroiG2p MstdroiG2r
            
            %     figure
            %     A = meanG1Aroi;
            %     B = meanG2Aroi;
            %     C = [A' B'];
            %     x = 1;
            %     for r = 1:size(roi{1,4},2)
            %         for rr = (r+1):size(roi{1,4},2)
            %             namesAroi{1,x} =  [roi{1,4}{1,r} '-' roi{1,4}{1,rr}];
            %             x = x + 1;
            %         end
            %     end
            %     X = categorical(namesAroi);
            %     p1 = bar(X,C,'stacked');
            %     p1(1).FaceColor = 'r';
            %     p1(2).FaceColor = 'b';
            %     ylabel('Pearson correlation');
            %     %savefig([savepath date '_SigCh']);
            %     %exportgraphics(gcf,[savepath 'SigCh.png'])
            %
            %     clear A B C x r rr X p1 ylabel
        end
    end
    %%%%% Différence de FC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if importROI | channelmode
    
        MATmeanG1G2 = MATmeanG1 - MATmeanG2;
        MATmeanG2G1 = MATmeanG2 - MATmeanG1;

        figure
        c = max([max(MATmeanG1G2,[],'all'), abs(min(MATmeanG1G2,[],'all'))]);
        cmin = -c + -c/10;
        cmax = c + c/10;
        clims = [cmin cmax];
        hG1G2 = imagesc(MATmeanG1G2,clims);
        colorbar
        colormap(jet)
        set(gca,'xtick', 1:size(MATmeanG1G2,1));
        set(gca,'xticklabel', 1:46);
        %xtickangle(90);
        set(gca,'ytick', 1:size(MATmeanG1G2,1));
        set(gca,'yticklabel', 1:46);
        ax = gca;
        ax.FontSize = 12;
        if importROI
            title('ROIs Mean G1-G2');
        elseif channelmode
            title('Channels Mean G1-G2')
        end
        savefig([savepath date '_MatG1G2']);
        exportgraphics(gcf,[savepath 'MatG1G2.png'])
        
        clear c cmin cmax clims ax hG1G2
        
        meandiff = mean(MATmeanG1G2,'all','omitnan');
        stddiff = std(MATmeanG1G2,0,'all','omitnan');
        tr = stddiff*2; %threshold de différence = 2 écarts type
        tf = (MATmeanG1G2>tr|MATmeanG1G2<-tr);
        MATG1G2 = MATmeanG1G2.*tf;

        figure
        c = max([max(MATG1G2,[],'all'), abs(min(MATG1G2,[],'all'))]);
        cmin = -c + -c/10;
        cmax = c + c/10;
        clims = [cmin cmax];
        hG1G2int = imagesc(MATG1G2,clims);
        colorbar
        colormap(jet)
        set(gca,'xtick', 1:size(MATG1G2,1));
        set(gca,'xticklabel', 1:46);
        %xtickangle(90);
        set(gca,'ytick', 1:size(MATG1G2,1));
        set(gca,'yticklabel', 1:46);
        ax = gca;
        ax.FontSize = 12;
        if importROI
            title('ROIs Mean G1-G2 > 2SD');
        elseif channelmode
            title('Channels Mean G1-G2 > 2SD')
        end
        savefig([savepath date '_MatSDG1G2']);
        exportgraphics(gcf,[savepath 'MatSDG1G2.png'])
        
        clear meandiff stddiff tr tf c cmin cmax clims ax hG1G2int

        datachG1G2 = MATvTBL.MAT2TBL(MATG1G2); %transformer la matrice en tableau
        tf = datachG1G2 ~= 0 & datachG1G2 > 0;
        idxpos = find(tf);
        tf = datachG1G2 ~= 0 & datachG1G2 < 0;
        idxneg = find(tf);

        figure
        hold on
        if importROI
            X = categorical(labelroi(idxpos));
        elseif channelmode
            X = categorical(labelch(idxpos));
        end
        bar(X,datachG1G2(1,idxpos), 'r')
        if importROI
            X = categorical(labelroi(idxneg));
        elseif channelmode
            X = categorical(labelch(idxneg));
        end
        bar(X,datachG1G2(1,idxneg), 'b')
        ylabel('Pearson correlation Difference G1-G2');
        if importROI
            xlabel('Pair of ROIs');
        elseif channelmode
            xlabel('Pair of Channels')
        end
        if importROI
            title('ROIs Largest FC difference between G1 and G2');
        elseif channelmode
            title('Channels Largest FC difference between G1 and G2')
        end
        savefig([savepath date '_TblSDG1G2']);
        exportgraphics(gcf,[savepath 'TblSDG1G2.png'])
        
        if find(MATG1G2)
            id = 1;
            List = strvcat(DATA{id}.ZoneList); %liste des paires SD
            ML = DATA{id}.zone.ml; %Loader S/D/ROI/Gr
            plotLst = DATA{id}.zone.plotLst;
            label =  DATA{id}.zone.label;
            plotconnectogram(fileorderconnectogram{1,1},MATG1G2,List,label,plotLst,ML)
            str = sprintf('%s Largest FC difference between G1 and G2');
            title(str)
            savefig([savepath date '_ConnectSDG1G2']);
            %exportgraphics(gcf,[savepath 'ConnectSDG1G2.png'])
        else
        end
        
        clear MATG1G2 tf idxpos idxneg X id List plotLst label str
        
    end
    
    
    if calculateROI
        for R = 1:numel(roi)
            ROIname = sprintf('%s',lgndroi{1,R});
            
            dataroiG1G2{:,R} = meanroiG1r{:,R} - meanroiG2r{:,R};
            dataroiG2G1{:,R} = meanroiG2r{:,R} - meanroiG1r{:,R};           
            %convertir les données en matrices
            MATroimeanG1G2{:,R} = MATvTBL.TBL2MAT(dataroiG1G2{:,R});
            
            figure
            c = max([max(MATroimeanG1G2{:,R},[],'all'), abs(min(MATroimeanG1G2{:,R},[],'all'))]);
            cmin = -c + -c/10;
            cmax = c + c/10;
            clims = [cmin cmax];
            hG1G2 = imagesc(MATroimeanG1G2{:,R},clims);
            colorbar
            colormap(jet)
            set(gca,'xtick', 1:size(MATroimeanG1G2{:,R},1));
            set(gca,'xticklabel', roi{1,R}(1,:));
            xtickangle(90);
            set(gca,'ytick', 1:size(MATroimeanG1G2{:,R},1));
            set(gca,'yticklabel', roi{1,R}(1,:));
            ax = gca;
            ax.FontSize = 12;
            str = sprintf('%s Mean G1-G2',ROIname);
            title(str)
            str = sprintf('_MatG1G2 %s',ROIname);
            savefig([savepath date str]);
            exportgraphics(gcf,[savepath str '.png'])
            clear c cmin cmax clims ax str hG1G2
            
            meandiff = mean(MATroimeanG1G2{:,R},'all','omitnan');
            stddiff = std(MATroimeanG1G2{:,R},0,'all','omitnan');
            tr = stddiff*2; %threshold de différence = 2 écarts type
            tf = (MATroimeanG1G2{:,R}>tr|MATroimeanG1G2{:,R}<-tr);
            MATG1G2 = MATroimeanG1G2{:,R}.*tf;
            
            figure
            c = max([max(MATG1G2,[],'all'), abs(min(MATG1G2,[],'all'))]);
            cmin = -c + -c/10;
            cmax = c + c/10;
            clims = [cmin cmax];
            hG1G2int = imagesc(MATG1G2,clims);
            colorbar
            colormap(jet)
            set(gca,'xtick', 1:size(MATG1G2,1));
            set(gca,'xticklabel', roi{1,R}(1,:));
            xtickangle(90);
            set(gca,'ytick', 1:size(MATG1G2,1));
            set(gca,'yticklabel', roi{1,R}(1,:));
            ax = gca;
            ax.FontSize = 12;
            str = sprintf('%s Mean G1-G2 > 2SD',ROIname);
            title(str)
            str = sprintf('_MatSDG1G2 %s',ROIname);
            savefig([savepath date str]);
            exportgraphics(gcf,[savepath str '.png'])
        
            clear meandiff stddiff tr tf c cmin cmax clims ax str hG1G2int
            
            DATG1G2 = MATvTBL.MAT2TBL(MATG1G2); %transformer la matrice en tableau
            tf = DATG1G2 ~= 0 & DATG1G2 > 0;
            idxpos = find(tf);
            tf = DATG1G2 ~= 0 & DATG1G2 < 0;
            idxneg = find(tf);
            
            figure
            hold on
            X = categorical(labelroiALL{1,R}(idxpos));
            bar(X,DATG1G2(1,idxpos), 'r')
            X = categorical(labelroiALL{1,R}(idxneg));
            bar(X,DATG1G2(1,idxneg), 'b')
            ylabel('Pearson correlation Difference G1-G2');
            xlabel('Pair of ROIs')
            str = sprintf('%s Largest FC difference between G1 and G2',ROIname);
            title(str);
            str = sprintf('_TblSDG1G2 %s',ROIname);
            savefig([savepath date str]);
            exportgraphics(gcf,[savepath str '.png'])
        
            if find(MATG1G2)
            plotLst = roi{:,R}(2,:);
            label =  roi{:,R}(1,:);
            plotconnectogramroi(fileorderconnectogram{:,R},MATG1G2,label,plotLst)
            str = sprintf('_ConnectSDG1G2 %s',ROIname);
            savefig([savepath date str]);
            %exportgraphics(gcf,[savepath str '.png'])
            else
            end 
            
            clear MATG1G2 DATG1G2 tf idxpos idxneg X str plotLst label  R ROIname
        end
    end
end

save([savepath 'workspace.mat'])

X = ['Results saved in ', savepath];
disp(X)
clear X
toc