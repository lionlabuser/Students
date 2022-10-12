%%%%%%%%%%%%%%GetToKnowYourDATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

savepathi = 'C:\data\Malnutrition\Resting\NIRS\Analyses\CORRmatrice0,01_0,08\PCAPW\CorrPairC0,1 ExcY\nofisher\';
load ([savepathi 'workspace.mat'])
%load ([datapath 'workspacemat.mat'])

disp(['Computing GetToKnowYourDATA on ' savepathi])

fileorderconnectogram = {'C:\data\Malnutrition\Resting\NIRS\Analyses\Connectogram_Mixte.txt',...
                         'C:\data\Malnutrition\Resting\NIRS\Analyses\Connectogram_Region.txt',...
                         'C:\data\Malnutrition\Resting\NIRS\Analyses\Connectogram_Fonction.txt',...
                         'C:\data\Malnutrition\Resting\NIRS\Analyses\Connectogram_Aire.txt'};

descrstatsmode = 1;
graphmode = 1;
channelmode = 1;
importROI = 0;
calculateROI = 1;

%%
%%Stats descriptives%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Computing Descriptive Statistics')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if descrstatsmode
    %Channels%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if importROI == 0 & channelmode == 1

        %normalité et valeurs extrêmes%
        stdchG1 = nanstd(datachG1);
        stdchG2 = nanstd(datachG2);
        sknchG1 = skewness(datachG1);
        sknchG2 = skewness(datachG2);
        krtchG1 = kurtosis(datachG1);
        krtchG2 = kurtosis(datachG2);

        tblgrmeanch = [array2table([datameanchG1],'VariableNames',labelch); array2table([datameanchG2],'VariableNames',labelch)];
        tblgrmeanch.Properties.RowNames = {'G1','G2'};
        tbldescrch = array2table([datameanchG1; datameanchG2; stdchG1; stdchG2; sknchG1; sknchG2; krtchG1; krtchG2],'VariableNames',labelch,'RowNames',{'meanG1','meanG2','stdG1','stdG2','sknG1','sknG2','krtG1','krtG2'});

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
        pbnanch =  ListCh(indexchnan(mostnanch));
        if isempty(pbnanch) == 0
           fprintf('Channel %.0f has more than 30 percent missing values\n',pbnanch)
        else
        end

        fig = figure;
        X = categorical(1:size(MATall,1));
        bar(X,cell2mat(nanfreq(:,1)))
        yline(mean(cell2mat(nanfreq(:,1)),1))
        ylabel('Total number of missing')
        title('Total number of missing channels')
        fig.WindowState = 'maximized';
        savefig([savepath date '_NANCh']);
        exportgraphics(gcf,[savepath 'NANCh.png'])
        close

        if calculateROI
            fig = figure;
            X = categorical(roi{1,R}(1,:));
            bar(X,cell2mat(nantot(:,1)))
            yline(mean(cell2mat(nantot(:,1)),1))
            ylabel('Total number of channels')
            title('Total number of channels with missing data in each region')
            pbaspect([1 1 1])
            fig.WindowState = 'maximized';
            savefig([savepath date '_NANChBYROI']);
            exportgraphics(gcf,[savepath 'NANChBYROI.png'])
            close

            clear nantot
        end

        clear R X x p c n idx tf nanch tfnanch idxnanch nb sortchnan indexchnan mostnanch pbnanch nanfreq
    end

    %%%%%Roi%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ROImode
        if importROI == 1
            %normalité et valeurs extrêmes%
            stdroiG1 = nanstd(dataroi(idG1,:));
            stdroiG2 = nanstd(dataroi(idG2,:));
            sknroiG1 = skewness(dataroi(idG1,:));
            sknroiG2 = skewness(dataroi(idG2,:));
            krtroiG1 = kurtosis(dataroi(idG1,:));
            krtroiG2 = kurtosis(dataroi(idG2,:));

            tblgrmeanroi = [array2table([datameanroiG1],'VariableNames',labelroi); array2table([datameanroiG2],'VariableNames',labelroi)];
            tblgrmeanroi.Properties.RowNames = {'G1','G2'};
            tbldescrroi = array2table([datameanroiG1; datameanroiG2; stdroiG1; stdroiG2; sknroiG1; sknroiG2; krtroiG1; krtroiG2],'VariableNames',labelroi,'RowNames',{'meanG1','meanG2','stdG1','stdG2','sknG1','sknG2','krtG1','krtG2'});

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

            fig = figure;
            X = categorical(Listroi,Listroi);
            bar(X,cell2mat(nanfreq(:,1)))
            yline(mean(cell2mat(nanfreq(:,1)),1))
            ylabel('Total number of missing')
            title('Total number of missing ROIs')
            fig.WindowState = 'maximized';
            savefig([savepath date '_NANROI']);
            exportgraphics(gcf,[savepath 'NANROI.png'])
            close

            clear x X p r n idx tf nanroi tfnanroi idxnanroi nb sortroinan indexroinan mostnanroi pbnanroi nanfreq
        end  

        if calculateROI ==1 & importROI == 0
            %normalité et valeurs extrêmes%
            ALLroi = [];
            for R = 1:numel(roi)
                ALLroi = [ALLroi, dataroiALL{1,R}]; %%Fait pour tous les ROIS calculés ensemble%%%%%%
            end

            meanALLroiG1 = nanmean(ALLroi(idG1,:));
            meanALLroiG2 = nanmean(ALLroi(idG2,:));
            stdALLroiG1 = nanstd(ALLroi(idG1,:));
            stdALLroiG2 = nanstd(ALLroi(idG2,:));
            sknALLroiG1 = skewness(ALLroi(idG1,:));
            sknALLroiG2 = skewness(ALLroi(idG2,:));
            krtALLroiG1 = kurtosis(ALLroi(idG1,:));
            krtALLroiG2 = kurtosis(ALLroi(idG2,:));

            tblgrmeanroi = [array2table([meanALLroiG1],'VariableNames',[labelroiALL{1,1:numel(roi)}]); array2table([meanALLroiG2],'VariableNames',[labelroiALL{1,1:numel(roi)}])];
            tblgrmeanroi.Properties.RowNames = {'G1','G2'};
            tbldescrroi = array2table([meanALLroiG1; meanALLroiG2; stdALLroiG1; stdALLroiG2; sknALLroiG1; sknALLroiG2; krtALLroiG1; krtALLroiG2],'VariableNames',[labelroiALL{1,1:numel(roi)}],'RowNames',{'meanG1','meanG2','stdG1','stdG2','sknG1','sknG2','krtG1','krtG2'});

            tf = tbldescrroi{{'sknG1','sknG2','krtG1','krtG2'},:} >2 | tbldescrroi{{'sknG1','sknG2','krtG1','krtG2'},:} <-2;
            n_abnormal = sum(any(tf));
            p_abnormal = n_abnormal/numel(tbldescrroi(1,:))*100;
            fprintf('%.1f percent of the ROIs variables are not respecting the normal distribution\n',p_abnormal);

            zALLroi = nanzscore(ALLroi);
            tf = zALLroi > 3.29 | zALLroi <-3.29;
            n_extreme = sum(any(tf));
            p_extreme = n_extreme/numel(zALLroi)*100;
            fprintf('%.1f percent of the ROIs variables have extreme values\n',p_extreme);

            clear X tf R ALLroi n_abnormal n_extreme p_abnormal p_extreme...
                stdchG1 stdchG2 sknchG1 sknchG2 krtchG1 krtchG2...
                meanroiG1 meanroiG2 stdroiG1 stdroiG2 sknroiG1 sknroiG2 krtroiG1 krtroiG2...
                meanALLroiG1 meanALLroiG2 stdALLroiG1 stdALLroiG2 sknALLroiG1 sknALLroiG2 krtALLroiG1 krtALLroiG2...
                zchALL zroiALLmean zALLroi
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
        
        %moyenne et std de la FC de chaque participant pour tous les canaux%
        meanchG1p = mean(datachG1,2,'omitnan');
        meanchG2p = mean(datachG2,2,'omitnan');
        MmeanchG1p = mean(meanchG1p);
        MmeanchG2p = mean(meanchG2p,'omitnan');
        stdchG1p = std(datachG1,0,2,'omitnan');
        stdchG2p = std(datachG2,0,2,'omitnan');
        MstdchG1p = mean(stdchG1p);
        MstdchG2p = mean(stdchG2p,'omitnan');

        %moyenne et std de la FC de chaque canal pour tous les participants%
        meanchG1c = mean(datachG1,'omitnan');
        meanchG2c = mean(datachG2,'omitnan');
        MmeanchG1c = mean(meanchG1c);
        MmeanchG2c = mean(meanchG2c);
        stdchG1c = std(datachG1,0,1,'omitnan');
        stdchG2c = std(datachG2,0,1,'omitnan');
        MstdchG1c = mean(stdchG1c);
        MstdchG2c = mean(stdchG2c);

        %moyenne et std de la FC pour tous les canaux et participants%
        MmeanchG1 = mean(datachG1,'all','omitnan');
        MmeanchG2 = mean(datachG2,'all','omitnan');
        MstdchG1 = std(datachG1,0,'all','omitnan');
        MstdchG2 = std(datachG2,0,'all','omitnan');

        fig = figure;
        hold on
        histg1 = histogram(meanchG1p);
        histg1.Normalization = 'pdf';
        histg1.BinWidth = 0.05;
        histg1.DisplayName = 'FC G1';
        xmin = floor((min(meanchG1p,[],'all'))*10)/10;
        xmax = ceil(max(meanchG1p,[],'all')*10)/10;
        x = xmin:0.05:xmax;
        pd = fitdist(meanchG1p,'Normal');
        y = pdf(pd,x);
        pg1 = plot(x,y,'LineWidth',1);
        pg1.Color = [0 0.4470 0.7410];
        histg2 = histogram(meanchG2p);
        histg2.Normalization = 'pdf';
        histg2.BinWidth = 0.05;
        histg2.DisplayName = 'FC G2';
        xmin = floor((min(meanchG2p,[],'all'))*10)/10;
        xmax = ceil(max(meanchG2p,[],'all')*10)/10;
        x = xmin:0.05:xmax;
        pd = fitdist(meanchG2p,'Normal');
        y = pdf(pd,x);
        pg2 = plot(x,y,'LineWidth',1);
        pg2.Color = [0.9290 0.6940 0.1250];
        dim = [.15 .6 .3 .3];
        str = sprintf('Mean G1 = %.3f\n Std G1 = %.3f\n Mean G2 = %.3f\n Std G2 = %.3f',MmeanchG1p,MstdchG1p,MmeanchG2p,MstdchG2p);
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
        legend;
        xlabel('Functional Connectivity Value');
        ylabel('Proportion of the participants (%)');
        title('Channels Histogram of the participants FC values')
        fig.WindowState = 'maximized';
        savefig([savepath date '_HistPartch']);
        exportgraphics(gcf,[savepath 'HistPartch.png'])
        close

        clear histg1 histg2 pg1 pg2 xmin xmax x pd y dim str xlabel ylabel

        fig = figure;
        hold on
        histg1 = histogram(meanchG1c);
        histg1.Normalization = 'pdf';
        histg1.BinWidth = 0.05;
        histg1.DisplayName = 'FC G1';
        xmin = floor((min(meanchG1c,[],'all'))*10)/10;
        xmax = ceil(max(meanchG1c,[],'all')*10)/10;
        x = xmin:0.05:xmax;
        pd = fitdist(meanchG1c','Normal');
        y = pdf(pd,x);
        pg1 = plot(x,y,'LineWidth',1);
        pg1.Color = [0 0.4470 0.7410];
        histg2 = histogram(meanchG2c);
        histg2.Normalization = 'pdf';
        histg2.BinWidth = 0.05;
        histg2.DisplayName = 'FC G2';
        xmin = floor((min(meanchG2c,[],'all'))*10)/10;
        xmax = ceil(max(meanchG2c,[],'all')*10)/10;
        x = xmin:0.05:xmax;
        pd = fitdist(meanchG2c','Normal');
        y = pdf(pd,x);
        pg2 = plot(x,y,'LineWidth',1);
        pg2.Color = [0.9290 0.6940 0.1250];
        dim = [.15 .6 .3 .3];
        str = sprintf('Mean G1 = %.3f\n Std G1 = %.3f\n Mean G2 = %.3f\n Std G2 = %.3f',MmeanchG1c,MstdchG1c,MmeanchG2c,MstdchG2c);
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
        legend;
        xlabel('Functional Connectivity Value');
        ylabel('Proportion of the Channels (%)');
        title('Channels Histogram of the Channels FC values')
        fig.WindowState = 'maximized';
        savefig([savepath date '_HistCHch']);
        exportgraphics(gcf,[savepath 'HistCHch.png'])
        close

        clear histg1 histg2 pg1 pg2 xmin xmax x pd y dim str xlabel ylabel

        fig = figure;
        hold on
        histg1 = histogram(datachG1);
        histg1.Normalization = 'pdf';
        histg1.BinWidth = 0.05;
        histg1.DisplayName = 'FC G1';
        limits = histg1.BinLimits;
        x = limits(1,1):histg1.BinWidth:limits(1,2);
        mu = MmeanchG1c;
        sigma = MstdchG1c;
        f = exp(-(x-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
        pg1 = plot(x,f,'LineWidth',1);
        pg1.Color = [0 0.4470 0.7410];
        histg2 = histogram(datachG2);
        histg2.Normalization = 'pdf';
        histg2.BinWidth = 0.05;
        histg2.DisplayName = 'FC G2';
        limits = histg2.BinLimits;
        x = limits(1,1):histg2.BinWidth:limits(1,2);
        mu = MmeanchG2c;
        sigma = MstdchG2c;
        f = exp(-(x-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
        pg2 = plot(x,f,'LineWidth',1);
        pg2.Color = [0.9290 0.6940 0.1250];
        legend;
        xlabel('Functional Connectivity Value');
        ylabel('Proportion of the FC (%)');
        title('Channels Histogram of the FC values')
        dim = [.15 .6 .3 .3];
        str = sprintf('Mean G1 = %.3f\n Std G1 = %.3f\n Mean G2 = %.3f\n Std G2 = %.3f',MmeanchG1c,MstdchG1c,MmeanchG2c,MstdchG2c);
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
        fig.WindowState = 'maximized';
        savefig([savepath date '_HistFCch']);
        exportgraphics(gcf,[savepath 'HistFCch.png'])
        close

        clear histg1 histg2 pg1 pg2 x mu sigma f dim str limits MmeanchG1 MmeanchG1p MmeanchG1c MmeanchG2 MmeanchG2p MmeanchG2c MstdchG1 MstdchG1p MstdchG1c MstdchG2 MstdchG2p MstdchG2c

    end

    %%%%%%%%%ROI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if importROI
        %moyenne et std de la FC de chaque participant pour toutes les régions%
        meanroiG1p = mean(dataroiG1,2,'omitnan');
        meanroiG2p = mean(dataroiG2,2,'omitnan');
        MmeanroiG1p = mean(meanroiG1p);
        MmeanroiG2p = mean(meanroiG2p,'omitnan');
        stdroiG1p = std(dataroiG1,0,2,'omitnan');
        stdroiG2p = std(dataroiG2,0,2,'omitnan');
        MstdroiG1p = mean(stdroiG1p);
        MstdroiG2p = mean(stdroiG2p,'omitnan');

        %moyenne et std de la FC de chaque région pour tous les participants%
        meanroiG1r = mean(dataroiG1,'omitnan');
        meanroiG2r = mean(dataroiG2,'omitnan');
        MmeanroiG1r = mean(meanroiG1r);
        MmeanroiG2r = mean(meanroiG2r);
        stdroiG1r = std(dataroiG1,0,1,'omitnan');
        stdroiG2r = std(dataroiG2,0,1,'omitnan');
        MstdroiG1r = mean(stdroiG1r);
        MstdroiG2r = mean(stdroiG2r);

        %moyenne et std de la FC pour toutes les régions et participants%
        MmeanroiG1 = mean(dataroiG1,'all','omitnan');
        MmeanroiG2 = mean(dataroiG2,'all','omitnan');
        MstdroiG1 = std(dataroiG1,0,'all','omitnan');
        MstdroiG2 = std(dataroiG2,0,'all','omitnan');

        %%HISTOGRAM%%%%%%%%%%%

        fig = figure;
        hold on
        histg1 = histogram(meanroiG1p);
        histg1.Normalization = 'pdf';
        histg1.BinWidth = 0.05;
        histg1.DisplayName = 'FC G1';
        xmin = floor((min(meanroiG1p,[],'all'))*10)/10;
        xmax = ceil(max(meanroiG1p,[],'all')*10)/10;
        x = xmin:0.05:xmax;
        pd = fitdist(meanroiG1p,'Normal');
        y = pdf(pd,x);
        pg1 = plot(x,y,'LineWidth',1);
        pg1.Color = [0 0.4470 0.7410];
        histg2 = histogram(meanroiG2p);
        histg2.Normalization = 'pdf';
        histg2.BinWidth = 0.05;
        histg2.DisplayName = 'FC G2';
        xmin = floor((min(meanroiG2p,[],'all'))*10)/10;
        xmax = ceil(max(meanroiG2p,[],'all')*10)/10;
        x = xmin:0.05:xmax;
        pd = fitdist(meanroiG2p,'Normal');
        y = pdf(pd,x);
        pg2 = plot(x,y,'LineWidth',1);
        pg2.Color = [0.9290 0.6940 0.1250];
        dim = [.15 .6 .3 .3];
        str = sprintf('Mean G1 = %.3f\n Std G1 = %.3f\n Mean G2 = %.3f\n Std G2 = %.3f',MmeanroiG1p,MstdroiG1p,MmeanroiG2p,MstdroiG2p);
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
        legend;
        xlabel('Functional Connectivity Value');
        ylabel('Proportion of the participants (%)');
        title('Roi Histogram of the participants FC values')
        fig.WindowState = 'maximized';
        savefig([savepath date '_HistROIPart']);
        exportgraphics(gcf,[savepath 'HistROIPart.png'])
        close

        clear histg1 histg2 pg1 pg2 xmin xmax x pd y dim str xlabel ylabel

        fig = figure;
        hold on
        histg1 = histogram(meanroiG1r);
        histg1.Normalization = 'pdf';
        histg1.BinWidth = 0.05;
        histg1.DisplayName = 'FC G1';
        xmin = floor((min(meanroiG1r,[],'all'))*10)/10;
        xmax = ceil(max(meanroiG1r,[],'all')*10)/10;
        x = xmin:0.05:xmax;
        pd = fitdist(meanroiG1r','Normal');
        y = pdf(pd,x);
        pg1 = plot(x,y,'LineWidth',1);
        pg1.Color = [0 0.4470 0.7410];
        histg2 = histogram(meanroiG2r);
        histg2.Normalization = 'pdf';
        histg2.BinWidth = 0.05;
        histg2.DisplayName = 'FC G2';
        xmin = floor((min(meanroiG2r,[],'all'))*10)/10;
        xmax = ceil(max(meanroiG2r,[],'all')*10)/10;
        x = xmin:0.05:xmax;
        pd = fitdist(meanroiG2r','Normal');
        y = pdf(pd,x);
        pg2 = plot(x,y,'LineWidth',1);
        pg2.Color = [0.9290 0.6940 0.1250];
        dim = [.15 .6 .3 .3];
        str = sprintf('Mean G1 = %.3f\n Std G1 = %.3f\n Mean G2 = %.3f\n Std G2 = %.3f',MmeanroiG1r,MstdroiG1r,MmeanroiG2r,MstdroiG2r);
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
        legend;
        xlabel('Functional Connectivity Value');
        ylabel('Proportion of the ROIs (%)');
        title('Roi Histogram of the ROIs FC values')
        fig.WindowState = 'maximized';
        savefig([savepath date '_HistROIRoi']);
        exportgraphics(gcf,[savepath 'HistROIRoi.png'])
        close

        clear histg1 histg2 pg1 pg2 xmin xmax x pd y dim str xlabel ylabel

        fig = figure;
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
        fig.WindowState = 'maximized';
        savefig([savepath date '_HistROIFC']);
        exportgraphics(gcf,[savepath 'HistROIFC.png'])
        close

        clear histg1 histg2 pg1 pg2 x mu sigma f dim str limits

        %%%HISTFIT%%%%%%%%%%%

        fig = figure;
        hold on
        histg1 = histfit(meanroiG1p);
        alpha(histg1,.5)
        histg1(2).Color = [0 0.4470 0.7410];
        histg2 = histfit(meanroiG2p);
        alpha(histg2,.5)
        histg2(2).Color = [0.9290 0.6940 0.1250];
        dim = [.15 .6 .3 .3];
        str = sprintf('Mean G1 = %.3f\n Std G1 = %.3f\n Mean G2 = %.3f\n Std G2 = %.3f',MmeanroiG1p,MstdroiG1p,MmeanroiG2p,MstdroiG2p);
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
        legend('FC G1','gaussian G1','FC G2','gaussian G2');
        xlabel('Functional Connectivity Value');
        ylabel('Proportion of the participants (%)');
        fig.WindowState = 'maximized';
        savefig([savepath date '_HistFitROIPart']);
        exportgraphics(gcf,[savepath 'HistFitROIPart.png'])
        close

        clear histg1 histg2 dim str xlabel ylabel

        fig = figure;
        hold on
        histg1 = histfit(meanroiG1r);
        alpha(histg1,.5)
        histg1(2).Color = [0 0.4470 0.7410];
        histg2 = histfit(meanroiG2r);
        alpha(histg2,.5)
        histg2(2).Color = [0.9290 0.6940 0.1250];
        dim = [.15 .6 .3 .3];
        str = sprintf('Mean G1 = %.3f\n Std G1 = %.3f\n Mean G2 = %.3f\n Std G2 = %.3f',MmeanroiG1r,MstdroiG1r,MmeanroiG2r,MstdroiG2r);
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
        legend('FC G1','gaussian G1','FC G2','gaussian G2');
        xlabel('Functional Connectivity Value');
        ylabel('Proportion of the ROIs (%)');
        title('Distribution of the FC values between groups')
        fig.WindowState = 'maximized';
        savefig([savepath date '_HistFitROIRoi']);
        exportgraphics(gcf,[savepath 'HistFitROIRoi.png'])
        close

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
            
            %moyenne et std de la FC de chaque participant pour tous les canaux%
            meanroiG1p{:,R} = nanmean(dataroiG1{:,R},2); %moyenne de la FC de chaque participant pour tous les canaux%
            meanroiG2p{:,R} = nanmean(dataroiG2{:,R},2);
            stdroiG1p{:,R} = std(dataroiG1{:,R},0,2,'omitnan');
            stdroiG2p{:,R} = std(dataroiG2{:,R},0,2,'omitnan');
            MmeanroiG1p{:,R} = mean(meanroiG1p{:,R});
            MmeanroiG2p{:,R} = mean(meanroiG2p{:,R},'omitnan');
            MstdroiG1p{:,R} = mean(stdroiG1p{:,R});
            MstdroiG2p{:,R} = mean(stdroiG2p{:,R},'omitnan');

            %moyenne et std de la FC de chaque région pour tous les participants%
            meanroiG1r{:,R} = nanmean(dataroiG1{:,R});%moyenne de la FC de chaque région pour tous les participants%
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
            
            fig = figure;
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
            fig.WindowState = 'maximized';
            savefig([savepath date '_HistROIPart']);
            exportgraphics(gcf,[savepath 'HistROIPart.png'])
            close
            
            clear histg1 histg2 pg1 pg2 xmin xmax x pd y dim str xlabel ylabel
            
            fig = figure;
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
            fig.WindowState = 'maximized';
            savefig([savepath date '_HistROIRoi']);
            exportgraphics(gcf,[savepath 'HistROIRoi.png'])
            close
            
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
        fig = figure;
        c = max([max(MATallG1G2,[],'all'), abs(min(MATallG1G2,[],'all'))]);
        cmin = -c + -c/10;
        cmax = c + c/10;
        clims = [cmin cmax];
        hG1G2 = imagesc(MATallG1G2,clims);
        colorbar
        colormap(jet)
        set(gca,'xtick', 1:size(MATallG1G2,1));
        set(gca,'xticklabel', ListCh);
        %xtickangle(90);
        set(gca,'ytick', 1:size(MATallG1G2,1));
        set(gca,'yticklabel', ListCh);
        ax = gca;
        ax.FontSize = 12;
        if importROI
            title('ROIs Mean G1-G2');
        elseif channelmode
            title('Channels Mean G1-G2')
        end
        pbaspect([1 1 1])
        fig.WindowState = 'maximized';
        savefig([savepath date '_MatG1G2']);
        exportgraphics(gcf,[savepath 'MatG1G2.png'])
        close
        
        clear c cmin cmax clims ax hG1G2
        
        meandiff = mean(MATallG1G2,'all','omitnan');
        stddiff = std(MATallG1G2,0,'all','omitnan');
        tr = stddiff*2; %threshold de différence = 2 écarts type
        tf = (MATallG1G2>tr|MATallG1G2<-tr);
        MATSDG1G2 = MATallG1G2.*tf;

        fig = figure;
        c = max([max(MATSDG1G2,[],'all'), abs(min(MATSDG1G2,[],'all'))]);
        cmin = -c + -c/10;
        cmax = c + c/10;
        clims = [cmin cmax];
        hG1G2int = imagesc(MATSDG1G2,clims);
        colorbar
        colormap(jet)
        set(gca,'xtick', 1:size(MATSDG1G2,1));
        set(gca,'xticklabel', ListCh);
        %xtickangle(90);
        set(gca,'ytick', 1:size(MATSDG1G2,1));
        set(gca,'yticklabel', ListCh);
        ax = gca;
        ax.FontSize = 12;
        if importROI
            title('ROIs Mean G1-G2 > 2SD');
        elseif channelmode
            title('Channels Mean G1-G2 > 2SD')
        end
        pbaspect([1 1 1])
        fig.WindowState = 'maximized';
        savefig([savepath date '_MatSDG1G2']);
        exportgraphics(gcf,[savepath 'MatSDG1G2.png'])
        close
        
        clear meandiff stddiff tr tf c cmin cmax clims ax hG1G2int

        DATSDG1G2 = MATvTBL.MAT2TBL(MATSDG1G2); %transformer la matrice en tableau
        tf = DATSDG1G2 ~= 0 & DATSDG1G2 > 0;
        idxpos = find(tf);
        tf = DATSDG1G2 ~= 0 & DATSDG1G2 < 0;
        idxneg = find(tf);

        fig = figure;
        hold on
        if importROI
            X = categorical(labelroi(idxpos));
        elseif channelmode
            X = categorical(labelch(idxpos));
        end
        bar(X,DATSDG1G2(1,idxpos), 'r')
        if importROI
            X = categorical(labelroi(idxneg));
        elseif channelmode
            X = categorical(labelch(idxneg));
        end
        bar(X,DATSDG1G2(1,idxneg), 'b')
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
        fig.WindowState = 'maximized';
        savefig([savepath date '_TblSDG1G2']);
        exportgraphics(gcf,[savepath 'TblSDG1G2.png'])
        close
        
        if find(MATSDG1G2)
            id = 1;
            List = ZoneList; %strvcat(DATA{id}.ZoneList); %liste des paires SD
            ML = DATA{id}.zone.ml; %Loader S/D/ROI/HBO-HBR
            plotLst = DATA{id}.zone.plotLst;
            label =  DATA{id}.zone.label;
            label = strrep(label, '_', ' ');
            plotconnectogram(fileorderconnectogram{1,1},MATSDG1G2(ListChinv,ListChinv),List,label,plotLst,ML)
            str = sprintf('%s Largest FC difference between G1 and G2');
            title(str)
            fig = gcf;
            pbaspect([1 1 1])
            fig.WindowState = 'maximized';
            savefig([savepath date '_ConnectSDG1G2']);
            exportgraphics(gcf,[savepath 'ConnectSDG1G2.png'])
            close
        else
        end
        
        clear MATSDG1G2 tf idxpos idxneg X id List ML plotLst label str
        
    end
    
    
    if calculateROI
        for R = 1:numel(roi)
            ROIname = sprintf('%s',lgndroi{1,R});
            
            fig = figure;
            c = max([max(MATroiG1G2{:,R},[],'all'), abs(min(MATroiG1G2{:,R},[],'all'))]);
            cmin = -c + -c/10;
            cmax = c + c/10;
            clims = [cmin cmax];
            hG1G2 = imagesc(MATroiG1G2{:,R},clims);
            colorbar
            colormap(jet)
            set(gca,'xtick', 1:size(MATroiG1G2{:,R},1));
            set(gca,'xticklabel', roi{1,R}(1,:));
            xtickangle(90);
            set(gca,'ytick', 1:size(MATroiG1G2{:,R},1));
            set(gca,'yticklabel', roi{1,R}(1,:));
            ax = gca;
            ax.FontSize = 12;
            str = sprintf('%s Mean G1-G2',ROIname);
            title(str)
            str = sprintf('_MatG1G2 %s',ROIname);
            pbaspect([1 1 1])
            fig.WindowState = 'maximized';
            savefig([savepath date str]);
            exportgraphics(gcf,[savepath str(2:end) '.png'])
            close
            clear c cmin cmax clims ax str hG1G2
            
            meandiff = mean(MATroiG1G2{:,R},'all','omitnan');
            stddiff = std(MATroiG1G2{:,R},0,'all','omitnan');
            tr = stddiff*2; %threshold de différence = 2 écarts type
            tf = (MATroiG1G2{:,R}>tr|MATroiG1G2{:,R}<-tr);
            MATroiSDG1G2 = MATroiG1G2{:,R}.*tf;
            
            fig = figure;
            c = max([max(MATroiSDG1G2,[],'all'), abs(min(MATroiSDG1G2,[],'all'))]);
            cmin = -c + -c/10;
            cmax = c + c/10;
            clims = [cmin cmax];
            hG1G2int = imagesc(MATroiSDG1G2,clims);
            colorbar
            colormap(jet)
            set(gca,'xtick', 1:size(MATroiSDG1G2,1));
            set(gca,'xticklabel', roi{1,R}(1,:));
            xtickangle(90);
            set(gca,'ytick', 1:size(MATroiSDG1G2,1));
            set(gca,'yticklabel', roi{1,R}(1,:));
            ax = gca;
            ax.FontSize = 12;
            str = sprintf('%s Mean G1-G2 > 2SD',ROIname);
            title(str)
            str = sprintf('_MatSDG1G2 %s',ROIname);
            pbaspect([1 1 1])
            fig.WindowState = 'maximized';
            savefig([savepath date str]);
            exportgraphics(gcf,[savepath str(2:end) '.png'])
            close
        
            clear meandiff stddiff tr tf c cmin cmax clims ax str hG1G2int
            
            DATroiSDG1G2 = MATvTBL.MAT2TBL(MATroiSDG1G2); %transformer la matrice en tableau
            tf = DATroiSDG1G2 ~= 0 & DATroiSDG1G2 > 0;
            idxpos = find(tf);
            tf = DATroiSDG1G2 ~= 0 & DATroiSDG1G2 < 0;
            idxneg = find(tf);
            
            fig = figure;
            hold on
            X = categorical(labelroiALL{1,R}(idxpos));
            bar(X,DATroiSDG1G2(1,idxpos), 'r')
            X = categorical(labelroiALL{1,R}(idxneg));
            bar(X,DATroiSDG1G2(1,idxneg), 'b')
            ylabel('Pearson correlation Difference G1-G2');
            xlabel('Pair of ROIs')
            str = sprintf('%s Largest FC difference between G1 and G2',ROIname);
            title(str);
            str = sprintf('_TblSDG1G2 %s',ROIname);
            fig.WindowState = 'maximized';
            savefig([savepath date str]);
            exportgraphics(gcf,[savepath str(2:end) '.png'])
            close
        
            if find(MATroiSDG1G2)
            plotLst = roi{:,R}(2,:);
            label =  roi{:,R}(1,:);
            plotconnectogramroi(fileorderconnectogram{:,R},MATroiSDG1G2,label,plotLst)
            str = sprintf('_ConnectSDG1G2 %s',ROIname);
            fig = gcf;
            pbaspect([1 1 1])
            fig.WindowState = 'maximized';
            savefig([savepath date str]);
            exportgraphics(gcf,[savepath str(2:end) '.png'])
            close
            else
            end 
            
            clear MATSDG1G2 DATSDG1G2 MATroiSDG1G2 DATroiSDG1G2 tf idxpos idxneg X str plotLst label  R ROIname
        end
    end
end

%save([savepath 'workspace.mat'])

X = ['Results saved in ', savepath];
disp(X)
clear all
toc