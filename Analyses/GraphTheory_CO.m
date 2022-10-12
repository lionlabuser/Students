%%%%%%%%%%%%%%%%%%%%%%%%%% Graph Theory%%%%%%%%%%%%%%%%%%%%
%tic
savepathi = 'C:\data\Malnutrition\Resting\NIRS\Analyses\CORR0,01_0,08\PCAPW\CorrPairC0,25 ExcY\nofisher\';

load ([savepathi 'workspace.mat'])
load ([savepathi 'workspacemat.mat'])


%% PARAMETERS TO MODIFY%%%%%%%%%%%
calculatemode = 1;
computegraphmode = 1;
graphmode = 1; %only if calculatemode = 1
analysemode = 1;
normalizemode = 1;
savemode = 1;
nointermode = 0;
fdrmode = 1;

weightedmode = 1;
binarizedmode = 1;
absolutethresh = 1;
proportionalthresh = 1;

nrand = 100;
nonegmode = 1
absmode = 0
rescalemode = 0

%BCT threshold
%Weigthed = garder les valeurs de corr, Binary = remplacer tout par 0 ou 1.%
itra = 0.01:0.01:0.85;% Seuil variable BCT soustraction
itrp = 0.01:0.01:0.85;
p = 0.05;
fixetr = 0.2; %pour les figures
thrrange = 1:85; %for analysis
factors = {gr,ses};
labelfactors = {'Group','SES'};

%%%%%%%%% END MODIFY %%%%%%%%%%%%%%

savepath = fullfile(savepath, 'GraphTheory');
if nonegmode
    savepath = [savepath '_Noneg'];
elseif absmode
    savepath = [savepath '_Abs'];
end
if rescalemode
    savepath = [savepath 'Resc'];
end

savepath = fullfile(savepath, filesep);
if ~isfolder(savepath)
    mkdir(savepath)
end

graphpath = fullfile(savepath,'Graphs', filesep);
if ~isfolder(graphpath)
    mkdir(graphpath)
end

%% Calcul Graph Theory
if calculatemode == 1
    if computegraphmode == 1
        disp(['Computing GraphTheory on ' savepathi])
    
        A = zeros(size(MATall));
        Amean = mean(MATall,3,'omitnan');
        for isubject = 1:size(MATall,3)
            Ai = squeeze(MATall(:,:,isubject)); %prendre la matrice de chaque participant%
            NaNch{:,isubject} = find(sum(Ai,'omitnan') == 0);
            %Ai(isnan(Ai)) = 0; %% Mettre les NAN à 0
            Ai(isnan(Ai)) = Amean(isnan(Ai)); %Remplacer les NaN par la moyenne
            if any(isnan(Ai),'all')
                fprintf('Warning NaNs left in the matrix for participant %s',num2str(isubject))
                Ai(isnan(Ai)) = 0; %% Mettre les NAN à 0
            end
            if nonegmode
                Ai(Ai<0) = 0;
            elseif absmode
                Ai = abs(Ai);
            end
            %         if rescalemode
            %            Ai = weight_conversion(Ai, 'normalize'); %Ai = rescale(Ai);
            %         end
            A(:,:,isubject) = Ai;
            Amaxmin = min(max(A,[],[1 2]));
            AMAX = max(A,[],'all');
            Aminmax = max(min(A,[],[1 2]));
            AMIN = min(A,[],'all');
        end
        if rescalemode
            A = weight_conversion(A, 'normalize'); %Ai = rescale(Ai);
        end
    end

    %% Absolute threshold, weighted
    if absolutethresh ==1 && weightedmode ==1 %(WTA)
        if computegraphmode
            disp('Computing absolute threshold & weighted metrics, running')
            LE_wta = zeros(size(MATall,1),size(MATall,3), numel(itra));
            LERRall_wta = zeros(size(MATall,1),nrand,size(MATall,3), numel(itra),'single');
            GE_wta = zeros(size(MATall,3), numel(itra));
            GERRall_wta = zeros(nrand,size(MATall,3),numel(itra),'single');
            LL_wta = zeros(size(MATall,3), numel(itra));
            LLRRall_wta = zeros(nrand,size(MATall,3),  numel(itra),'single');
            CC_wta = zeros(size(MATall,1),size(MATall,3), numel(itra));
            CCRRall_wta = zeros(size(MATall,1), nrand, size(MATall,3), numel(itra),'single');
            K_wta = zeros(size(MATall,1),size(MATall,3), numel(itra));
            KRRall_wta = zeros(size(MATall,1), nrand, size(MATall,3), numel(itra),'single');
            RRi_wta = zeros(size(MATall,1),size(MATall,2),nrand,'single');
            if savemode == 1
                A_wta = zeros(size(MATall,1),size(MATall,2),numel(itra),size(MATall,3),'single');
                RR_wta = zeros(size(MATall,1),size(MATall,2),nrand,numel(itra),size(MATall,3),'single');
                D_wta = zeros(size(MATall,1),size(MATall,2),size(MATall,3),numel(itra),'single');
                DRRall_wta = zeros(size(MATall,1),size(MATall,2),nrand,size(MATall,3),numel(itra),'single');
            end
    
            for isubject = 1:size(MATall,3)
                fprintf('\t Computing subject %.0f %s\n',isubject,part(isubject));
                fprintf(1,'\t \t Computing threshold:     \n');
                for idtr = 1:numel(itra)
                    itrv = itra(idtr);
                    str = num2str(itrv);
                    if length(str) < 4
                        if length(str) == 2
                            str = append(str,'  ');
                        elseif length(str) == 3
                            str = append(str,' ');
                        end
                    end
                    fprintf(1,'\b\b\b\b%s',str);
    
                    Ai = squeeze(A(:,:,isubject)); %prendre la matrice de chaque participant%
                    Ai_wta = threshold_absolute(Ai, itrv); %Thresholder la matrice (weigthed)%
                    if savemode == 1
                        A_wta(:,:,idtr,isubject) = Ai_wta;
                        %A_wta(isnan(A_wta)) = 0;  %% Mettre les NAN à 0
                    end
    
                    for indrand = 1:nrand %% For each subject, 100 random networks that preserved the degree distribution of the original network were estimated%
                        if ~any(Ai_wta(:)) %if the matrix does not contain any connexion (all 0)
                            fprintf(1,'\n \t \t Matrix %s is all zeros. Rdm matrix set as NaN.\n', num2str(indrand))
                            RRi_wta(:,:,:) = NaN;
                            if itrv~=itra(end)
                                fprintf(1,'\t Computing threshold:     ')
                            end
                            if savemode == 1
                                RR_wta(:,:,:,idtr,isubject) = nan(size(RRi_wta));
                                %RR_wta (:,:, indrand) = NaN;
                            end
                            break
                        end
                        RRi_wta(:,:,indrand) = randmio_und(Ai_wta,1); %% random matrix
                        if isnan(RRi_wta(:,:,indrand)) %(RR_wta (:,:, indrand)) %if the random matrix did not convert and was set as NaN
                            fprintf(1,'\n \t \t Rdm matrix # %s wont converge. Set as NaN.\n', num2str(indrand))
                            RRi_wta(:,:,:) = NaN;
                            if itrv~=itra(end)
                                fprintf(1,'\t Computing threshold:     ')
                            end
                            if savemode == 1
                                RR_wta(:,:,:,idtr,isubject) = nan(size(RRi_wta));
                                %RR_wta (:,:, indrand) = NaN;
                            end
                            break
                        end
                        if savemode == 1
                            RR_wta(:,:,indrand,idtr,isubject) = RRi_wta(:,:,indrand);
                        end
                    end
    
                    RRi_wta(isnan(RRi_wta)) = 0; %Mettre les NAN à 0
    
                    K_wta(:,isubject, idtr) = degrees_und(Ai_wta);
                    LE_wta(:,isubject, idtr) = efficiency_wei(Ai_wta,2); %calculer le LE pour chaque canal, de chaque matrice binarisée de chaque participant%
                    GE_wta(isubject, idtr) = efficiency_wei(Ai_wta,0); %calculer le GE pour chaque matrice binarisée de chaque participant%
                    D = distance_wei(weight_conversion(Ai_wta,'lengths')); %calculer la matrice de distance%
                    [lambda,~,~,~,~] = charpath(D,0,0);
                    LL_wta(isubject, idtr) = lambda; %calculer le LL pour chaque matrice de distance binarisée de chaque participant%
                    if savemode == 1
                        D_wta(:,:,isubject,idtr) = D;
                    end
                    CC_wta(:,isubject, idtr) = clustering_coef_wu(Ai_wta); %calculer le CC pour chaque canal, de chaque matrice binarisée de chaque participant%
    
                    for indrand = 1:nrand
                        %RRi_wta(:,:,indrand) = squeeze(RR_wta(:,:,indrand,idtr,isubject));
                        KRRall_wta(:,indrand,isubject,idtr) = degrees_und(RRi_wta(:,:,indrand));
                        LERRall_wta(:,indrand,isubject,idtr) = efficiency_wei(RRi_wta(:,:,indrand),2); %calculer le LE pour chaque canal, de chaque matrice de chaque participant%
                        GERRall_wta(indrand,isubject,idtr) = efficiency_wei(RRi_wta(:,:,indrand)); %calculer le GE pour chaque matrice de chaque participant%
                        DRR = distance_wei(weight_conversion(RRi_wta(:,:,indrand),'lengths'));
                        [lambdaRR,~,~,~,~] = charpath(DRR,0,0);
                        LLRRall_wta(indrand, isubject,idtr) = lambdaRR; %calculer le LL pour chaque matrice de distance de chaque participant%
                        if savemode == 1
                            DRRall_wta(:,:,indrand,isubject,idtr) = DRR;
                        end
                        CCRRall_wta(:,indrand,isubject,idtr) = clustering_coef_wu(RRi_wta(:,:,indrand)); %calculer le CC pour chaque canal, de chaque matrice de chaque participant%
                    end
                end
                fprintf(1,'\n')
            end
            clear RRi_wta
    
            if savemode == 1
                if ~isfile([savepath date '_RR.mat'])
                    save([savepath date '_RR.mat'],'RR_wta','-v7.3')
                else
                    save([savepath date '_RR.mat'],'RR_wta','-append')
                end
                if ~isfile([savepath date '_A.mat'])
                    save([savepath date '_A.mat'],'A_wta','-v7.3')
                else
                    save([savepath date '_A.mat'],'A_wta','-append')
                end
                if ~isfile([savepath date '_D.mat'])
                    save([savepath date '_D.mat'],'D_wta','DRRall_wta','-v7.3')
                else
                    save([savepath date '_D.mat'],'D_wta','DRRall_wta','-append')
                end
                clear RR_wta A_wta D_wta DRRall_wta
            end
            
            disp('Computing absolute threshold & weighted metrics, done')
        end

        %% Metriques%

        %%Local Efficiency Input Distance matrix, une valeur par canal, par threshold, par participant%
        fprintf('\t Local Efficiency_wta, running\n')

        LEmean_wta = squeeze(mean(LE_wta,1)); %faire la moyenne du LE de tous les canaux, pour chaque matrice de chaque participant%

        LERRmean_wta = squeeze(mean(LERRall_wta,2)); %compute the mean of the 100 random metrics
        LERRstd_wta = squeeze(std(LERRall_wta,0,2)); %compute the std of the 100 random metrics
        LEmeanRRall_wta = squeeze(mean(LERRmean_wta,1)); %compute the mean of the channels

        %LENormall_wta = LE_wta./LERRmean_wta; %divide the real metric by the mean random metric of each participant
        LENormall_wta = LE_wta./mean(LERRmean_wta,2,'omitnan'); %divide the real metric of each participant by the mean random metric
        LENormall_wta(isinf(LENormall_wta)) = NaN;
        LEmeanNormall_wta = squeeze(mean(LENormall_wta,1,'omitnan'));
        LENormallgood_wta = LENormall_wta;
        LENormallgood_wta(isnan(LENormallgood_wta)) = 0;
        LEmeanNormallgood_wta = LEmeanNormall_wta;
        LEmeanNormallgood_wta(isnan(LEmeanNormallgood_wta)) = 0;
        
        ZLE_wta = (LE_wta-LERRmean_wta)./LERRstd_wta;
        ZLE_wta(isinf(ZLE_wta)) = NaN;
        ZLEmean_wta = squeeze(mean(ZLE_wta,1,'omitnan'));
        ZLEgood_wta = ZLE_wta;
        ZLEgood_wta(isnan(ZLEgood_wta)) = 0;
        ZLEmeangood_wta = ZLEmean_wta;
        ZLEmeangood_wta(isnan(ZLEmeangood_wta)) = 0;

        fprintf('\t Local Efficiency_wta, done\n')

        %%Global efficiency, une valeur par threshold par participant
        fprintf('\t Global Efficiency_wta, running\n')

        GERRmean_wta = squeeze(mean(GERRall_wta,1));
        GERRstd_wta = squeeze(std(GERRall_wta,0,1));

        %GENormall_wta = GE_wta./GERRmean_wta;
        GENormall_wta = GE_wta./mean(GERRmean_wta,2,'omitnan');
        GENormall_wta(isinf(GENormall_wta)) = NaN;
        GENormallgood_wta = GENormall_wta;
        GENormallgood_wta(isnan(GENormallgood_wta)) = 0;

        ZGE_wta = (GE_wta-GERRmean_wta)./GERRstd_wta;
        ZGE_wta(isinf(ZGE_wta)) = NaN;
        ZGEgood_wta = ZGE_wta;
        ZGEgood_wta(isnan(ZGEgood_wta)) = 0;

        fprintf('\t Global Efficiency_wta, done\n')

        %%Characteristic pathlenght charpath(D)(BU, BD, WU, WD networks), une valeur par threshold par participant.
        fprintf('\t Characteristic path length_wta, running\n')
        
        LLgood_wta = LL_wta;
        LLgood_wta(isnan(LLgood_wta)) = 0;
        LLRRmean_wta = squeeze(mean(LLRRall_wta,1,'omitnan'));
        LLRRstd_wta = squeeze(std(LLRRall_wta,0,1,'omitnan'));

        LLNormall_wta = LL_wta./mean(LLRRmean_wta,2,'omitnan');
        LLNormall_wta(isinf(LLNormall_wta)) = NaN;
        LLNormallgood_wta = LLNormall_wta;
        LLNormallgood_wta(isnan(LLNormallgood_wta)) = 0;

        ZLL_wta = (LL_wta-LLRRmean_wta)./LLRRstd_wta;
        ZLL_wta(isinf(ZLL_wta)) = NaN;
        ZLLgood_wta = ZLL_wta;
        ZLLgood_wta(isnan(ZLLgood_wta)) = 0;

        fprintf('\t Characteristic path length_wta, done\n')

        %%Clustering Coefficient clustering_coef_bu(A), une valeur par canal, par threshold par participant
        fprintf('\t Clustering coefficient_wta, running\n')

        CCmean_wta = squeeze(mean(CC_wta,1)); %faire la moyenne du CC de tous les canaux, pour chaque matrice binarisée de chaque participant%
        CCRRmean_wta = squeeze(mean(CCRRall_wta,2));
        CCRRstd_wta = squeeze(std(CCRRall_wta,0,2));
        CCmeanRRall_wta = squeeze(mean(CCRRmean_wta,1));

        CCNormall_wta = CC_wta./mean(CCRRmean_wta,2,'omitnan');
        CCNormall_wta(isinf(CCNormall_wta)) = NaN;
        CCmeanNormall_wta = squeeze(mean(CCNormall_wta,1,'omitnan'));
        CCNormallgood_wta = CCNormall_wta;
        CCNormallgood_wta(isnan(CCNormallgood_wta)) = 0;
        CCmeanNormallgood_wta = CCmeanNormall_wta;
        CCmeanNormallgood_wta(isnan(CCmeanNormallgood_wta)) = 0;
        
        ZCC_wta = (CC_wta-CCRRmean_wta)./CCRRstd_wta;
        ZCC_wta(isinf(ZCC_wta)) = NaN;
        ZCCmean_wta = squeeze(mean(ZCC_wta,1,'omitnan'));
        ZCCgood_wta = ZCC_wta;
        ZCCgood_wta(isnan(ZCCgood_wta)) = 0;
        ZCCmeangood_wta = ZCCmean_wta;
        ZCCmeangood_wta(isnan(ZCCmeangood_wta)) = 0;

        fprintf('\t Clustering coefficient_wta, done\n')

        %%Small worldness index%
        fprintf('\t Small Worldness index, running\n')

        SW_wta = CCmean_wta./LL_wta;
        SW_wta(isinf(SW_wta)) = NaN;
        SWgood_wta = SW_wta;
        SWgood_wta(isnan(SWgood_wta)) = 0;

        SWNormall_wta = CCmeanNormall_wta./LLNormall_wta;
        SWNormall_wta(isinf(SWNormall_wta)) = NaN;
        SWNormallgood_wta = SWNormall_wta;
        SWNormallgood_wta(isnan(SWNormallgood_wta)) = 0;

        ZSW_wta = ZCCmean_wta./ZLL_wta;
        ZSW_wta(isinf(ZSW_wta)) = NaN;
        ZSWgood_wta = ZSW_wta;
        ZSWgood_wta(isnan(ZSWgood_wta)) = 0;

        SWmean_wta = mean(SWNormall_wta,1,'omitnan');
        SWmeanone_wta = find(SWmean_wta >= 1,1,'first');

        fprintf('\t Small Worldness index, done\n')

        %%Degree
        fprintf('\t Degree_wta, running\n')

        Kmean_wta = squeeze(mean(K_wta,1)); %compute the mean K of all channels
        KmeanG1_wta = squeeze(mean(K_wta(:,idG1,:),1));
        KmeanG2_wta = squeeze(mean(K_wta(:,idG2,:),1));

        if graphmode == 1
        fig = figure;
        hold on
        histg1 = histogram(KmeanG1_wta);
        histg1.Normalization = 'pdf';
        %histg1.BinWidth = 0.05;
        histg1.DisplayName = 'Mean degree G1';
        xmin = floor((min(KmeanG1_wta,[],'all'))*10)/10;
        xmax = ceil(max(KmeanG1_wta,[],'all')*10)/10;
        x = xmin:0.05:xmax;
        pd = fitdist(KmeanG1_wta(:),'Normal');
        y = pdf(pd,x);
        pg1 = plot(x,y,'LineWidth',1);
        pg1.Color = [0 0.4470 0.7410];
        histg2 = histogram(KmeanG2_wta);
        histg2.Normalization = 'pdf';
        %histg2.BinWidth = 0.05;
        histg2.DisplayName = 'Mean degree G2';
        xmin = floor((min(KmeanG2_wta,[],'all'))*10)/10;
        xmax = ceil(max(KmeanG2_wta,[],'all')*10)/10;
        x = xmin:0.05:xmax;
        pd = fitdist(KmeanG2_wta(:),'Normal');
        y = pdf(pd,x);
        pg2 = plot(x,y,'LineWidth',1);
        pg2.Color = [0.9290 0.6940 0.1250];
        %dim = [.15 .6 .3 .3];
        %str = sprintf('Mean G1 = %.3f\n Std G1 = %.3f\n Mean G2 = %.3f\n Std G2 = %.3f',MmeanchG1p,MstdchG1p,MmeanchG2p,MstdchG2p);
        %annotation('textbox',dim,'String',str,'FitBoxToText','on');
        legend;
        xlabel('Degree Value');
        ylabel('Proportion of the participants (%)');
        title('Histogram of the participants degree values for wta')
        fig.WindowState = 'maximized';
        savefig([graphpath date '_HistK_wta']);
        exportgraphics(gcf,[graphpath 'HistK_wta.png'])
        close
        end

        fprintf('\t Degree_wta, done\n')

        %%Criteria
        for isubject = 1:size(K_wta,2)
            goodch = ~ismember(1:size(K_wta,1), NaNch{1,isubject});
            Kmin_wta(isubject,:) = min(squeeze(K_wta(goodch,isubject,:)));
            CCmin_wta(isubject,:) = min(squeeze(CC_wta(goodch,isubject,:)));
            try
            Kone_wta(isubject,:) = find(Kmin_wta(isubject,:) > 0,1,'last');
            catch
                disp(['Not finding any K>0 for subject ' num2str(isubject)])
               Kone_wta(isubject,1) = nan;
            end
            try
            CCone_wta(isubject,:) = find(CCmin_wta(isubject,:) > 0,1,'last'); 
            catch
               disp(['Not finding any CC>0 for subject ' num2str(isubject)])
               CCone_wta(isubject,1) = nan;
            end
        end

        Kmeanone_wta = min(Kone_wta);
        CCmeanone_wta = min(CCone_wta);
        Kmeanmin_wta = min(Kmean_wta);
        CCmeanmin_wta = min(CCmean_wta);
        KminM_wta = mean(Kmin_wta);
        CCminM_wta = mean(CCmin_wta);
        CCmeanNormallM_wta = mean(CCmeanNormall_wta,'omitnan');
        LLNormallM_wta = mean(LLNormall_wta,'omitnan');

        if graphmode == 1
            figure
            hold on
            plot(KminM_wta,'r','displayname','Degree');
            plot(CCminM_wta+1,'b','displayname','ClustCoeff');
            plot(SWmean_wta,'g','displayname','SmallWorld');
            ylim([0.5 2])
            yline(1)
            xlabel('Threshold','fontsize',12)
            ylabel('Value','fontsize',12)
            legend
            box on
            %fig.WindowState = 'maximized';
            savefig([graphpath date '_Plotmin_wta']);
            exportgraphics(gcf,[graphpath 'Plotmin_wta.png'])
            close
        end
        
                clear b D d1 d1R d2 d2R DRR f id idtr indrand isubject itrv lambda lambdaRR meand1 meand2 tr x y ...
            Ai Ai_wta RRi_wta RRmeani_wta ...
            fig histg1 xmin xmax x y pd pg1  histg2 pg2 dim

        %%save%
        if ~isfile([savepath date '_workspace.mat'])
            save([savepath date '_workspace.mat'],'-v7.3')
        else
            save([savepath date '_workspace.mat'],'-append')
        end
        fprintf('Data saved in %s on %s\n', savepath, date);

         %% AUC
        for ich = 1:size(LE_wta,1)
            AUCLE_wta(:,ich) = trapz(squeeze(LE_wta(ich,:,thrrange)),2);
        end
        AUCLEmean_wta = trapz(LEmean_wta(:,thrrange),2);
        for ich = 1:size(LENormallgood_wta,1)
            AUCLENormall_wta(:,ich) = trapz(squeeze(LENormallgood_wta(ich,:,thrrange)),2);
        end
        AUCLEmeanNormall_wta = trapz(LEmeanNormallgood_wta(:,thrrange),2);
        for ich = 1:size(ZLEgood_wta,1)
            AUCZLE_wta(:,ich) = trapz(squeeze(ZLEgood_wta(ich,:,thrrange)),2);
        end
        AUCZLEmean_wta = trapz(ZLEmeangood_wta(:,thrrange),2);
        AUCGE_wta = trapz(GE_wta(:,thrrange),2);
        AUCGENormall_wta = trapz(GENormallgood_wta(:,thrrange),2);
        AUCZGE_wta = trapz(ZGEgood_wta(:,thrrange),2);
        AUCLL_wta = trapz(LLgood_wta(:,thrrange),2);
        AUCLLNormall_wta = trapz(LLNormallgood_wta(:,thrrange),2);
        AUCZLL_wta = trapz(ZLLgood_wta(:,thrrange),2);
        for ich = 1:size(CC_wta,1)
            AUCCC_wta(:, ich) = trapz(squeeze(CC_wta(ich,:,thrrange)),2);
        end
        AUCCCmean_wta = trapz(CCmean_wta(:,thrrange),2);
        for ich = 1:size(CCNormallgood_wta,1)
            AUCCCNormall_wta(:,ich) = trapz(squeeze(CCNormallgood_wta(ich,:,thrrange)),2);
        end
        AUCCCmeanNormall_wta = trapz(CCmeanNormallgood_wta(:,thrrange),2);
        for ich = 1:size(ZCCgood_wta,1)
            AUCZCC_wta(:,ich) = trapz(squeeze(ZCCgood_wta(ich,:,thrrange)),2);
        end
        AUCZCCmean_wta = trapz(ZCCmeangood_wta(:,thrrange),2);
        AUCSW_wta = trapz(SWgood_wta(:,thrrange),2);
        AUCSWNormall_wta = trapz(SWNormallgood_wta(:,thrrange),2);
        AUCZSW_wta = trapz(ZSWgood_wta(:,thrrange),2);
        AUCKmean_wta = trapz(Kmean_wta(:,thrrange),2);

clear ich

        %%save%
        if ~isfile([savepath date '_workspace.mat'])
            save([savepath date '_workspace.mat'],'-v7.3')
        else
            save([savepath date '_workspace.mat'],'-append')
        end
        fprintf('Data saved in %s on %s\n', savepath, date);
    end

    %% Absolute threshold, binarized
    if absolutethresh == 1 && binarizedmode ==1 %BTA
        if computegraphmode
            disp('Computing absolute threshold & binarized metrics, running')
    
            LE_bta = zeros(size(MATall,1),size(MATall,3), numel(itra));
            LERRall_bta = zeros(size(MATall,1), nrand, size(MATall,3), numel(itra),'single');
            GE_bta = zeros(size(MATall,3), numel(itra));
            GERRall_bta = zeros(nrand, size(MATall,3), numel(itra),'single');
            LL_bta = zeros(size(MATall,3), numel(itra));
            LLRRall_bta = zeros(nrand, size(MATall,3),  numel(itra),'single');
            CC_bta = zeros(size(MATall,1),size(MATall,3), numel(itra));
            CCRRall_bta = zeros(size(MATall,1), nrand, size(MATall,3), numel(itra),'single');
            K_bta = zeros(size(MATall,1),size(MATall,3), numel(itra));
            KRRall_bta = zeros(size(MATall,1), nrand, size(MATall,3), numel(itra),'single');
            RRi_bta = zeros(size(MATall,1),size(MATall,2),nrand,'single');
            if savemode == 1
                A_bta = zeros(size(MATall,1),size(MATall,2),numel(itra),size(MATall,3),'single');
                RR_bta = zeros(size(MATall,1),size(MATall,2),nrand,numel(itra),size(MATall,3),'single'); %RR_bta = zeros(size(MATall,1),size(MATall,2),100);
                D_bta = zeros(size(MATall,1),size(MATall,2),size(MATall,3),numel(itra),'single');
                DRRall_bta = zeros(size(MATall,1),size(MATall,2),nrand,size(MATall,3),numel(itra),'single');
            end
    
            for isubject = 1:size(MATall,3)
                fprintf('\t Computing subject %.0f %s\n',isubject,part(isubject));
                fprintf(1,'\t \t Computing threshold:     \n');
                for idtr = 1:numel(itra)
                    itrv = itra(idtr);
                    str = num2str(itrv);
                    if length(str) < 4
                        if length(str) == 2
                            str = append(str,'  ');
                        elseif length(str) == 3
                            str = append(str,' ');
                        end
                    end
                    fprintf(1,'\b\b\b\b%s',str);
                    
                    Ai = squeeze(A(:,:,isubject)); %A %prendre la matrice de chaque participant%
                    Ai_wta = threshold_absolute(Ai, itrv); %A_wta %Thresholder la matrice (weigthed)% %Ai_wta = A_wta(:,:,idtr,isubject);
                    Ai_bta = weight_conversion(Ai_wta,'binarize'); %A_bta %Binarizer la matrice (binarized
                    if savemode == 1
                        A_bta(:,:,idtr,isubject) = Ai_bta;
                        %A_bta(isnan(A_bta)) = 0;  %% Mettre les NAN à 0
                    end
    
                    for indrand = 1:nrand %% For each subject, 100 random networks that preserved the degree distribution of the original network were estimated%
                        if ~any(Ai_bta(:))
                            fprintf(1,'\n \t \t Matrix %s is all zeros. Rdm matrix set as NaN.\n', num2str(indrand))
                            RRi_bta(:,:,:) = NaN;
                            if itrv~=itra(end)
                                fprintf(1,'\t Computing threshold:     ')
                            end
                            if savemode == 1
                                RR_bta(:,:,:,idtr,isubject) = nan(size(RRi_bta)); %RR_bta
                            end
                            break
                        end
                        RRi_bta(:,:,indrand) = randmio_und(Ai_bta,1); %RR_bta(:,:,indran) %% random matrix
                        if isnan(RRi_bta(:,:,indrand))
                            fprintf(1,'\n \t \t Rdm matrix # %s wont converge. Set as NaN.\n', num2str(indrand))
                            RRi_bta(:,:,:) = NaN;
                            if itrv~=itra(end)
                                fprintf(1,'\t Computing threshold:     ')
                            end
                            if savemode == 1
                                RR_bta(:,:,:,idtr,isubject) = nan(size(RRi_bta)); %RR_bta
                            end
                            break
                        end
                        if savemode == 1
                            RR_bta(:,:,indrand,idtr,isubject) = RRi_bta(:,:,indrand); %RR_bta(:,:,indran) %% random matrix
                        end
                    end
    
                    RRi_bta(isnan(RRi_bta)) = 0; %Mettre les NAN à 0
    
                    K_bta(:,isubject, idtr) = degrees_und(Ai_bta);
                    LE_bta(:,isubject, idtr) = efficiency_bin(Ai_bta,1); %calculer le LE pour chaque canal, de chaque matrice binarisée de chaque participant%
                    GE_bta(isubject, idtr) = efficiency_bin(Ai_bta); %calculer le GE pour chaque matrice de chaque participant%
                    D = distance_bin(Ai_bta); %calculer la matrice de distance%
                    [lambda,~,~,~,~] = charpath(D,0,0);
                    if savemode == 1
                        D_bta(:,:,isubject, idtr) = D;
                    end
                    LL_bta(isubject, idtr) = lambda; %calculer le LL pour chaque matrice de distance binarisée de chaque participant%
                    CC_bta(:,isubject, idtr) = clustering_coef_bu(Ai_bta); %calculer le CC pour chaque canal, de chaque matrice de chaque participant%
    
                    for indrand = 1:nrand
                        %RRi_bta(:,:,indrand) = squeeze(RR_bta(:,:,indrand,idtr,isubject));
                        KRRall_bta(:,indrand,isubject,idtr) = degrees_und(RRi_bta(:,:,indrand));
                        LERRall_bta(:,indrand,isubject,idtr) = efficiency_bin(RRi_bta(:,:,indrand),1); %calculer le LE pour chaque canal, de chaque matrice de chaque participant%
                        GERRall_bta(indrand,isubject,idtr) = efficiency_bin(RRi_bta(:,:,indrand)); %calculer le GE pour chaque matrice de chaque participant%
                        DRR = distance_bin(RRi_bta(:,:,indrand));
                        [lambdaRR,~,~,~,~] = charpath(DRR,0,0);
                        if savemode == 1
                            DRRall_bta(:,:,indrand, idtr, isubject) = DRR;
                        end
                        LLRRall_bta(indrand, isubject,idtr) = lambdaRR; %calculer le LL pour chaque matrice de distance de chaque participant%
                        CCRRall_bta(:,indrand,isubject,idtr) = clustering_coef_bu(RRi_bta(:,:,indrand)); %calculer le CC pour chaque canal, de chaque matrice de chaque participant%
                    end
                end
                fprintf(1,'\n')
            end
            clear RRi_bta
    
            if savemode == 1
                if ~isfile([savepath date '_RR.mat'])
                    save([savepath date '_RR.mat'],'RR_bta','-v7.3')
                else
                    save([savepath date '_RR.mat'],'RR_bta','-append')
                end
                if ~isfile([savepath date '_A.mat'])
                    save([savepath date '_A.mat'],'A_bta','-v7.3')
                else
                    save([savepath date '_A.mat'],'A_bta','-append')
                end
                if ~isfile([savepath date '_D.mat'])
                    save([savepath date '_D.mat'],'D_bta','DRRall_bta','-v7.3')
                else
                    save([savepath date '_D.mat'],'D_bta','DRRall_bta','-append')
                end
                clear RR_bta A_bta D_bta DRRall_bta
            end
    
            disp('Computing absolute threshold & binarized metrics, done')
        end

        %% Metriques
        % Local Efficiency Input Distance matrix, une valeur par canal, par threshold, par participant%
        fprintf(' \t Local Efficiency_bta, running\n');

        %faire la moyenne du LE de tous les canaux, pour chaque matrice de chaque participant%
        LEmean_bta = squeeze(mean(LE_bta,1));
        LERRmean_bta = squeeze(mean(LERRall_bta,2));
        LERRstd_bta = squeeze(std(LERRall_bta,0,2));
        LEmeanRRall_bta = squeeze(mean(LERRmean_bta,1));

        LENormall_bta = LE_bta./mean(LERRmean_bta,2,'omitnan');
        LENormall_bta(isinf(LENormall_bta)) = NaN;
        LEmeanNormall_bta = squeeze(mean(LENormall_bta,1,'omitnan'));
        LENormallgood_bta = LENormall_bta;
        LENormallgood_bta(isnan(LENormallgood_bta)) = 0;
        LEmeanNormallgood_bta = LEmeanNormall_bta;
        LEmeanNormallgood_bta(isnan(LEmeanNormallgood_bta)) = 0;
        
        ZLE_bta = (LE_bta-LERRmean_bta)./LERRstd_bta;
        ZLE_bta(isinf(ZLE_bta)) = NaN;
        ZLEmean_bta = squeeze(mean(ZLE_bta,1,'omitnan'));
        ZLEgood_bta = ZLE_bta;
        ZLEgood_bta(isnan(ZLEgood_bta)) = 0;
        ZLEmeangood_bta = ZLEmean_bta;
        ZLEmeangood_bta(isnan(ZLEmeangood_bta)) = 0;

        fprintf('\t Local Efficiency_bta, done\n');

        %%Global efficiency, une valeur par threshold par participant
        fprintf('\t Global Efficiency_bta, running\n');

        GERRmean_bta = squeeze(mean(GERRall_bta,1));
        GERRstd_bta = squeeze(std(GERRall_bta,0,1));

        GENormall_bta = GE_bta./mean(GERRmean_bta,2,'omitnan');
        GENormall_bta(isinf(GENormall_bta)) = NaN;
        GENormallgood_bta = GENormall_bta;
        GENormallgood_bta(isnan(GENormallgood_bta)) = 0;

        ZGE_bta = (GE_bta-GERRmean_bta)./GERRstd_bta;
        ZGE_bta(isinf(ZGE_bta)) = NaN;
        ZGEgood_bta = ZGE_bta;
        ZGEgood_bta(isnan(ZGEgood_bta)) = 0;

        fprintf('\t Global Efficiency_bta, done\n');

        %%Characteristic pathlenght charpath(D)(BU, BD, WU, WD networks), une valeur par threshold par participant.
        fprintf('\t Characteristic path length_bta, running\n');

        LLgood_bta = LL_bta;
        LLgood_bta(isnan(LLgood_bta)) = 0;
        LLRRmean_bta = squeeze(mean(LLRRall_bta,1,'omitnan'));
        LLRRstd_bta = squeeze(std(LLRRall_bta,0,1,'omitnan'));

        LLNormall_bta = LL_bta./mean(LLRRmean_bta,2,'omitnan');
        LLNormall_bta(isinf(LLNormall_bta)) = NaN;
        LLNormallgood_bta = LLNormall_bta;
        LLNormallgood_bta(isnan(LLNormallgood_bta)) = 0;

        ZLL_bta = (LL_bta-LLRRmean_bta)./LLRRstd_bta;
        ZLL_bta(isinf(ZLL_bta)) = NaN;
        ZLLgood_bta = ZLL_bta;
        ZLLgood_bta(isnan(ZLLgood_bta)) = 0;

        fprintf('\t Characteristic path length_bta, done\n');

        %%Clustering Coefficient clustering_coef_bu(A), une valeur par canal, par threshold par participant
        fprintf('\t Clustering coefficient_bta, running\n');

        CCmean_bta = squeeze(mean(CC_bta,1)); %faire la moyenne du CC de tous les canaux, pour chaque matrice binarisée de chaque participant%
        CCRRmean_bta = squeeze(mean(CCRRall_bta,2));
        CCRRstd_bta = squeeze(std(CCRRall_bta,0,2));
        CCmeanRRall_bta = squeeze(mean(CCRRmean_bta,1));
        
        CCNormall_bta = CC_bta./mean(CCRRmean_bta,2,'omitnan');
        CCNormall_bta(isinf(CCNormall_bta)) = NaN;
        CCmeanNormall_bta = squeeze(mean(CCNormall_bta,1,'omitnan'));
        CCNormallgood_bta = CCNormall_bta;
        CCNormallgood_bta(isnan(CCNormallgood_bta)) = 0;
        CCmeanNormallgood_bta = CCmeanNormall_bta;
        CCmeanNormallgood_bta(isnan(CCmeanNormallgood_bta)) = 0;
        
        ZCC_bta = (CC_bta-CCRRmean_bta)./CCRRstd_bta;
        ZCC_bta(isinf(ZCC_bta)) = NaN;
        ZCCmean_bta = squeeze(mean(ZCC_bta,1,'omitnan'));
        ZCCgood_bta = ZCC_bta;
        ZCCgood_bta(isnan(ZCCgood_bta)) = 0;
        ZCCmeangood_bta = ZCCmean_bta;
        ZCCmeangood_bta(isnan(ZCCmeangood_bta)) = 0;

        fprintf('\t Clustering coefficient_bta, done\n');

        %%Small worldness index%
        fprintf('\t Small Worldness index_bta, running\n');

        SW_bta = CCmean_bta./LL_bta;
        SW_bta(isinf(SW_bta)) = NaN;
        SWgood_bta = SW_bta;
        SWgood_bta(isnan(SWgood_bta )) = 0;

        SWNormall_bta = CCmeanNormall_bta./LLNormall_bta;
        SWNormall_bta(isinf(SWNormall_bta)) = NaN;
        SWNormallgood_bta = SWNormall_bta;
        SWNormallgood_bta(isnan(SWNormallgood_bta)) = 0;

        ZSW_bta = ZCCmean_bta./ZLL_bta;
        ZSW_bta(isinf(ZSW_bta)) = NaN;
        ZSWgood_bta = ZSW_bta;
        ZSWgood_bta(isnan(ZSWgood_bta)) = 0;

        SWmean_bta = mean(SWNormall_bta,1,'omitnan');
        SWmeanone_bta = find(SWmean_bta >= 1,1,'first');

        fprintf('\t Small Worldness index_bta, done\n');

        %%Degree
        fprintf('\t Degree_bta, running\n')

        Kmean_bta = squeeze(mean(K_bta,1));
        KmeanG1_bta = squeeze(mean(K_bta(:,idG1,:),1));
        KmeanG2_bta = squeeze(mean(K_bta(:,idG2,:),1));

        if graphmode == 1
        fig = figure;
        hold on
        histg1 = histogram(KmeanG1_bta);
        histg1.Normalization = 'pdf';
        %histg1.BinWidth = 0.05;
        histg1.DisplayName = 'Mean degree G1';
        xmin = floor((min(KmeanG1_bta,[],'all'))*10)/10;
        xmax = ceil(max(KmeanG1_bta,[],'all')*10)/10;
        x = xmin:0.05:xmax;
        pd = fitdist(KmeanG1_bta(:),'Normal');
        y = pdf(pd,x);
        pg1 = plot(x,y,'LineWidth',1);
        pg1.Color = [0 0.4470 0.7410];
        histg2 = histogram(KmeanG2_bta);
        histg2.Normalization = 'pdf';
        %histg2.BinWidth = 0.05;
        histg2.DisplayName = 'Mean degree G2';
        xmin = floor((min(KmeanG2_bta,[],'all'))*10)/10;
        xmax = ceil(max(KmeanG2_bta,[],'all')*10)/10;
        x = xmin:0.05:xmax;
        pd = fitdist(KmeanG2_bta(:),'Normal');
        y = pdf(pd,x);
        pg2 = plot(x,y,'LineWidth',1);
        pg2.Color = [0.9290 0.6940 0.1250];
        %dim = [.15 .6 .3 .3];
        %str = sprintf('Mean G1 = %.3f\n Std G1 = %.3f\n Mean G2 = %.3f\n Std G2 = %.3f',MmeanchG1p,MstdchG1p,MmeanchG2p,MstdchG2p);
        %annotation('textbox',dim,'String',str,'FitBoxToText','on');
        legend;
        xlabel('Degree Value');
        ylabel('Proportion of the participants (%)');
        title('Histogram of the participants degree values for bta')
        fig.WindowState = 'maximized';
        savefig([graphpath date '_HistK_bta']);
        exportgraphics(gcf,[graphpath 'HistK_bta.png'])
        close
        end

        fprintf('\t Degree_bta, done\n')

        %%Criteria
        for isubject = 1:size(K_bta,2)
            goodch = ~ismember(1:size(K_bta,1), NaNch{1,isubject});
            Kmin_bta(isubject,:) = min(squeeze(K_bta(goodch,isubject,:)));
            CCmin_bta(isubject,:) = min(squeeze(CC_bta(goodch,isubject,:)));
            try
            Kone_bta(isubject,:) = find(Kmin_bta(isubject,:) > 0,1,'last');
            catch
               disp(['Not finding any K>0 for subject ' num2str(isubject)])
               Kone_bta(isubject,1) = nan;
            end
            try
            CCone_bta(isubject,:) = find(CCmin_bta(isubject,:) > 0,1,'last');
            catch
                disp(['Not finding any CC>0 for subject ' num2str(isubject)])
               CCone_bta(isubject,1) = nan;
            end
        end
        %             Kmin_bta(:,idtr) = min(K_bta(:,:,idtr),[],[1 2]);
        %             CCmin_bta(:,idtr) = min(CC_bta(:,:,idtr),[],[1 2]);

        Kmeanone_bta = min(Kone_bta);
        CCmeanone_bta = min(CCone_bta);
        Kmeanmin_bta = min(Kmean_bta);
        CCmeanmin_bta = min(CCmean_bta);
        KminM_bta = mean(Kmin_bta);
        CCminM_bta = mean(CCmin_bta);
        CCmeanNormallM_bta = mean(CCmeanNormall_bta,'omitnan');
        LLNormallM_bta = mean(LLNormall_bta,'omitnan');

        if graphmode == 1
            figure
            hold on
            plot(KminM_bta,'r','displayname','Degree');
            plot(CCminM_bta+1,'b','displayname','ClustCoeff');
            plot(SWmean_bta,'g','displayname','SmallWorld');
            ylim([0.5 2])
            yline(1)
            xlabel('Threshold','fontsize',12)
            ylabel('Value','fontsize',12)
            legend
            box on
            %fig.WindowState = 'maximized';
            savefig([graphpath date '_Plotmin_bta']);
            exportgraphics(gcf,[graphpath 'Plotmin_bta.png'])
            close
        end
        
        clear b D d1 d1R d2 d2R DRR f id idtr indrand isubject itrv lambda lambdaRR meand1 meand2 tr x y ...
            Ai Ai_wta Ai_bta RRi_bta RRmeani_bta ...
            fig histg1 xmin xmax x y pd pg1  histg2 pg2 dim

        if ~isfile([savepath date '_workspace.mat'])
            save([savepath date '_workspace.mat'],'-v7.3')
        else
            save([savepath date '_workspace.mat'],'-append')
        end
        fprintf('Data saved in %s on %s\n', savepath, date);

        %% AUC
        for ich = 1:size(LE_bta,1)
            AUCLE_bta(:,ich) = trapz(squeeze(LE_bta(ich,:,thrrange)),2);
        end
        AUCLEmean_bta = trapz(LEmean_bta(:,thrrange),2);
        for ich = 1:size(LENormallgood_bta,1)
            AUCLENormall_bta(:,ich) = trapz(squeeze(LENormallgood_bta(ich,:,thrrange)),2);
        end
        AUCLEmeanNormall_bta = trapz(LEmeanNormallgood_bta(:,thrrange),2);
        for ich = 1:size(ZLEgood_bta,1)
            AUCZLE_bta(:,ich) = trapz(squeeze(ZLEgood_bta(ich,:,thrrange)),2);
        end
        AUCZLEmean_bta = trapz(ZLEmeangood_bta(:,thrrange),2);
        AUCGE_bta = trapz(GE_bta(:,thrrange),2);
        AUCGENormall_bta = trapz(GENormallgood_bta(:,thrrange),2);
        AUCZGE_bta = trapz(ZGEgood_bta(:,thrrange),2);
        AUCLL_bta = trapz(LLgood_bta(:,thrrange),2);
        AUCLLNormall_bta = trapz(LLNormallgood_bta(:,thrrange),2);
        AUCZLL_bta = trapz(ZLLgood_bta(:,thrrange),2);
        for ich = 1:size(CC_bta,1)
            AUCCC_bta(:,ich) = trapz(squeeze(CC_bta(ich,:,thrrange)),2);
        end
        AUCCCmean_bta = trapz(CCmean_bta(:,thrrange),2);
        for ich = 1:size(CCNormallgood_bta,1)
            AUCCCNormall_bta(:,ich) = trapz(squeeze(CCNormallgood_bta(ich,:,thrrange)),2);
        end
        AUCCCmeanNormall_bta = trapz(CCmeanNormallgood_bta(:,thrrange),2);
        for ich = 1:size(ZCCgood_bta,1)
            AUCZCC_bta(:,ich) = trapz(squeeze(ZCCgood_bta(ich,:,thrrange)),2);
        end
        AUCZCCmean_bta = trapz(ZCCmeangood_bta(:,thrrange),2);
        AUCSW_bta = trapz(SWgood_bta(:,thrrange),2);
        AUCSWNormall_bta = trapz(SWNormallgood_bta(:,thrrange),2);
        AUCZSW_bta = trapz(ZSWgood_bta(:,thrrange),2);
        AUCKmean_bta = trapz(Kmean_bta(:,thrrange),2);
        
        clear ich

        if ~isfile([savepath date '_workspace.mat'])
            save([savepath date '_workspace.mat'],'-v7.3')
        else
            save([savepath date '_workspace.mat'],'-append')
        end
        fprintf('Data saved in %s on %s\n', savepath, date);
    end

    %% Proportional threshold & weighted
    if proportionalthresh ==1 && weightedmode ==1 %WTP
        if computegraphmode
            disp('Computing proportional threshold & weighted metrics, running')
    
            LE_wtp = zeros(size(MATall,1),size(MATall,3), numel(itrp));
            LERRall_wtp = zeros(size(MATall,1), nrand, size(MATall,3), numel(itrp),'single');
            GE_wtp = zeros(size(MATall,3), numel(itrp));
            GERRall_wtp = zeros(nrand, size(MATall,3), numel(itrp),'single');
            LL_wtp = zeros(size(MATall,3), numel(itrp));
            LLRRall_wtp = zeros(nrand, size(MATall,3), numel(itrp),'single');
            CC_wtp = zeros(size(MATall,1),size(MATall,3),numel(itrp));
            CCRRall_wtp = zeros(size(MATall,1), nrand, size(MATall,3),numel(itrp),'single');
            K_wtp = zeros(size(MATall,1),size(MATall,3),numel(itrp));
            KRRall_wtp = zeros(size(MATall,1), nrand, size(MATall,3),numel(itrp),'single');
            RRi_wtp = zeros(size(MATall,1),size(MATall,2),nrand,'single');
            if savemode == 1
                A_wtp = zeros(size(MATall,1),size(MATall,2),numel(itrp),size(MATall,3),'single');
                RR_wtp = zeros(size(MATall,1),size(MATall,2),nrand,numel(itrp),size(MATall,3),'single'); %RR_wtp = zeros(size(MATall,1),size(MATall,2),100);
                D_wtp = zeros(size(MATall,1),size(MATall,2),size(MATall,3),numel(itrp),'single');
                DRRall_wtp = zeros(size(MATall,1),size(MATall,2),nrand,size(MATall,3),numel(itrp),'single');
            end
    
            for isubject = 1:size(MATall,3)
                fprintf('\t Computing subject %.0f %s\n',isubject, part(isubject));
                fprintf(1,'\t \t Computing threshold:     \n');
                for idtr = 1:numel(itrp)
                    itrv = itrp(idtr);
                    str = num2str(itrv);
                    if length(str) < 4
                        if length(str) == 2
                            str = append(str,'  ');
                        elseif length(str) == 3
                            str = append(str,' ');
                        end
                    end
                    fprintf(1,'\b\b\b\b%s',str);
                    Ai = squeeze(A(:,:,isubject)); %prendre la matrice de chaque participant%
                    Ai_wtp = threshold_proportional(Ai,itrv); %A_wtp %Thresholder la matrice (weigthed)%
                    if savemode == 1
                        A_wtp(:,:,idtr,isubject) = Ai_wtp;
                        %A_wtp(isnan(A_wtp)) = 0;  %% Mettre les NAN à 0
                    end
    
                    for indrand = 1:nrand %% For each subject, 100 random networks that preserved the degree distribution of the original network were estimated%
                        if ~any(Ai_wtp(:))
                            fprintf(1,'\n \t \t Matrix %s is all zeros. Rdm matrix set as NaN.\n', num2str(indrand))
                            RRi_wtp(:,:,:) = NaN;
                            if itrv~=itrp(end)
                                fprintf(1,'\t Computing threshold:     ')
                            end
                            if savemode == 1
                                RR_wtp(:,:,:,idtr,isubject) = nan(size(RRi_wtp));
                            end
                            break
                        end
                        RRi_wtp(:,:,indrand) = randmio_und(Ai_wtp,1); %RR_wtp(:,:,indrand) %% random matrix
                        if isnan(RRi_wtp(:,:,indrand)) %RR_wtp(:,:,indrand)
                            fprintf(1,'\n \t \t Rdm matrix # %s wont converge. Set as NaN.\n', num2str(indrand))
                            RRi_wtp(:,:,:) = NaN;
                            if itrv~=itrp(end)
                                fprintf(1,'\t Computing threshold:     ')
                            end
                            if savemode == 1
                                RR_wtp(:,:,:,idtr,isubject) = nan(size(RRi_wtp));
                            end
                            break
                        end
                        if savemode == 1
                            RR_wtp(:,:,indrand,idtr,isubject) = RRi_wtp(:,:,indrand); %RR_wtp(:,:,indrand) %% random matrix
                        end
                    end
    
                    RRi_wtp(isnan(RRi_wtp)) = 0; %Mettre les NAN à 0
                    K_wtp(:,isubject, idtr) = degrees_und(Ai_wtp);
                    LE_wtp(:,isubject, idtr) = efficiency_wei(Ai_wtp,2); %calculer le LE pour chaque canal, de chaque matrice binarisée de chaque participant%
                    GE_wtp(isubject, idtr) = efficiency_wei(Ai_wtp,0); %calculer le GE pour chaque matrice de chaque participant%
                    D = distance_wei(weight_conversion(Ai_wtp,'lengths')); %calculer la matrice de distance% %%Pas certaine du bon input pour la fonction de distance
                    [lambda,~,~,~,~] = charpath(D,0,0);
                    if savemode == 1
                        D_wtp(:,:,isubject,idtr) = D;
                    end
                    LL_wtp(isubject, idtr) = lambda; %calculer le LL pour chaque matrice de distance binarisée de chaque participant%
                    CC_wtp(:,isubject, idtr) = clustering_coef_wu(Ai_wtp); %calculer le CC pour chaque canal, de chaque matrice de chaque participant%
    
                    for indrand = 1:nrand
                        %RRi_wtp(:,:,indrand) = squeeze(RR_wtp(:,:,indrand,idtr,isubject));
                        KRRall_wtp(:,indrand,isubject,idtr) = degrees_und(RRi_wtp(:,:,indrand));
                        LERRall_wtp(:,indrand,isubject,idtr) = efficiency_wei(RRi_wtp(:,:,indrand),2); %calculer le LE pour chaque canal, de chaque matrice de chaque participant%
                        GERRall_wtp(indrand,isubject,idtr) = efficiency_wei(RRi_wtp(:,:,indrand)); %calculer le GE pour chaque matrice de chaque participant%
                        DRR = distance_wei(weight_conversion(RRi_wtp(:,:,indrand),'lengths'));
                        [lambdaRR,~,~,~,~] = charpath(DRR,0,0);
                        if savemode == 1
                        DRRall_wtp (:,:,indrand, isubject, idtr) = DRR;
                        end
                        LLRRall_wtp(indrand, isubject,idtr) = lambdaRR; %calculer le LL pour chaque matrice de distance de chaque participant%
                        CCRRall_wtp(:,indrand,isubject,idtr) = clustering_coef_wu(RRi_wtp(:,:,indrand)); %calculer le CC pour chaque canal, de chaque matrice de chaque participant%
                    end
                end
                fprintf(1,'\n')
            end
            clear RRi_wtp
    
            if savemode == 1
                if ~isfile([savepath date '_RR.mat'])
                    save([savepath date '_RR.mat'],'RR_wtp','-v7.3')
                else
                    save([savepath date '_RR.mat'],'RR_wtp','-append')
                end
                if ~isfile([savepath date '_A.mat'])
                    save([savepath date '_A.mat'],'A_wtp','-v7.3')
                else
                    save([savepath date '_A.mat'],'A_wtp','-append')
                end
                if ~isfile([savepath date '_D.mat'])
                    save([savepath date '_D.mat'],'D_wtp','DRRall_wtp','-v7.3')
                else
                    save([savepath date '_D.mat'],'D_wtp','DRRall_wtp','-append')
                end
                clear RR_wtp A_wtp D_wtp DRRall_wtp
            end
            disp('Computing proportional threshold & weighted metrics, done')
        end

        %% Metriques
        %Local Efficiency Input Distance matrix, une valeur par canal, par threshold, par participant%
        fprintf('\t Local Efficiency_wtp, running\n')

        LEmean_wtp = squeeze(mean(LE_wtp,1)); %faire la moyenne du LE de tous les canaux, pour chaque matrice de chaque participant%
        LERRmean_wtp = squeeze(mean(LERRall_wtp,2));
        LERRstd_wtp = squeeze(std(LERRall_wtp,0,2));
        LEmeanRRall_wtp = squeeze(mean(LERRmean_wtp,1));

        LENormall_wtp = LE_wtp./mean(LERRmean_wtp,2,'omitnan');
        LENormall_wtp(isinf(LENormall_wtp)) = NaN;
        LEmeanNormall_wtp = squeeze(mean(LENormall_wtp,1,'omitnan'));
        LENormallgood_wtp = LENormall_wtp;
        LENormallgood_wtp(isnan(LENormallgood_wtp)) = 0;
        LEmeanNormallgood_wtp = LEmeanNormall_wtp;
        LEmeanNormallgood_wtp(isnan(LEmeanNormallgood_wtp)) = 0;
        
        ZLE_wtp = (LE_wtp-LERRmean_wtp)./LERRstd_wtp;
        ZLE_wtp(isinf(ZLE_wtp)) = NaN;
        ZLEmean_wtp = squeeze(mean(ZLE_wtp,1,'omitnan'));
        ZLEgood_wtp = ZLE_wtp;
        ZLEgood_wtp(isnan(ZLEgood_wtp)) = 0;
        ZLEmeangood_wtp = ZLEmean_wtp;
        ZLEmeangood_wtp(isnan(ZLEmeangood_wtp)) = 0;

        fprintf('\t Local Efficiency_wtp, done\n')

        %%Global efficiency, une valeur par threshold par participant
        fprintf('\t Global Efficiency_wtp, running\n')

        GERRmean_wtp = squeeze(mean(GERRall_wtp,1));
        GERRstd_wtp = squeeze(std(GERRall_wtp,0,1));

        GENormall_wtp = GE_wtp./mean(GERRmean_wtp,2,'omitnan');
        GENormall_wtp(isinf(GENormall_wtp)) = NaN;
        GENormallgood_wtp = GENormall_wtp;
        GENormallgood_wtp(isnan(GENormallgood_wtp)) = 0;

        ZGE_wtp = (GE_wtp-GERRmean_wtp)./GERRstd_wtp;
        ZGE_wtp(isinf(ZGE_wtp)) = NaN;
        ZGEgood_wtp = ZGE_wtp;
        ZGEgood_wtp(isnan(ZGEgood_wtp)) = 0;

        fprintf('\t Global Efficiency_wtp, done\n')

        %%Characteristic pathlenght charpath(D)(BU, BD, WU, WD networks), une valeur par threshold par participant.
        fprintf('\t Characteristic path length_wtp, running\n')

        LLRRmean_wtp = squeeze(mean(LLRRall_wtp,1));
        LLRRstd_wtp = squeeze(std(LLRRall_wtp,0,1));

        LLNormall_wtp = LL_wtp./mean(LLRRmean_wtp,2,'omitnan');
        LLNormall_wtp(isinf(LLNormall_wtp)) = NaN;
        LLNormallgood_wtp = LLNormall_wtp;
        LLNormallgood_wtp(isnan(LLNormallgood_wtp)) = 0;

        ZLL_wtp = (LL_wtp-LLRRmean_wtp)./LLRRstd_wtp;
        ZLL_wtp(isinf(ZLL_wtp)) = NaN;
        ZLLgood_wtp = ZLL_wtp;
        ZLLgood_wtp(isnan(ZLLgood_wtp)) = 0;
        
        fprintf('\tCharacteristic path length_wtp, done\n')

        %%Clustering Coefficient clustering_coef_bu(A), une valeur par canal, par threshold par participant
        fprintf('\t Clustering coefficient_wtp, running\n')

        CCmean_wtp = squeeze(mean(CC_wtp,1)); %faire la moyenne du CC de tous les canaux, pour chaque matrice binarisée de chaque participant%
        CCRRmean_wtp = squeeze(mean(CCRRall_wtp,2));
        CCRRstd_wtp = squeeze(std(CCRRall_wtp,0,2));
        CCmeanRRall_wtp = squeeze(mean(CCRRmean_wtp,1));
        
        CCNormall_wtp = CC_wtp./mean(CCRRmean_wtp,2,'omitnan');
        CCNormall_wtp(isinf(CCNormall_wtp)) = NaN;
        CCmeanNormall_wtp = squeeze(mean(CCNormall_wtp,1,'omitnan'));
        CCNormallgood_wtp = CCNormall_wtp;
        CCNormallgood_wtp(isnan(CCNormallgood_wtp)) = 0;
        CCmeanNormallgood_wtp = CCmeanNormall_wtp;
        CCmeanNormallgood_wtp(isnan(CCmeanNormallgood_wtp)) = 0;

        ZCC_wtp = (CC_wtp-CCRRmean_wtp)./CCRRstd_wtp;
        ZCC_wtp(isinf(ZCC_wtp)) = NaN;
        ZCCmean_wtp = squeeze(mean(ZCC_wtp,1,'omitnan'));
        ZCCgood_wtp = ZCC_wtp;
        ZCCgood_wtp(isnan(ZCCgood_wtp)) = 0;
        ZCCmeangood_wtp = ZCCmean_wtp;
        ZCCmeangood_wtp(isnan(ZCCmeangood_wtp)) = 0;

        fprintf('\t Clustering coefficient_wtp, done\n')

        %%Small worldness index%
        fprintf('\t Small Worldness index_wtp, running\n')

        SW_wtp = CCmean_wtp./LL_wtp;
        SW_wtp(isinf(SW_wtp)) = NaN;
        SWgood_wtp = SW_wtp;
        SWgood_wtp(isnan(SWgood_wtp)) = 0;

        SWNormall_wtp = CCmeanNormall_wtp./LLNormall_wtp;
        SWNormall_wtp(isinf(SWNormall_wtp)) = NaN;
        SWNormallgood_wtp = SWNormall_wtp;
        SWNormallgood_wtp(isnan(SWNormallgood_wtp)) = 0;

        ZSW_wtp = ZCCmean_wtp./ZLL_wtp;
        ZSW_wtp(isinf(ZSW_wtp)) = NaN;
        ZSWgood_wtp = ZSW_wtp;
        ZSWgood_wtp(isnan(ZSWgood_wtp)) = 0;

        SWmean_wtp = mean(SWNormall_wtp,1,'omitnan');
        SWmeanone_wtp = find(SWmean_wtp < 1,1,'first')-1;
        fprintf('\t Small Worldness index_wtp, done\n')

        %%Degree
        fprintf('\t Degree_wtp, running\n')

        Kmean_wtp = squeeze(mean(K_wtp,1));
        KmeanG1_wtp = squeeze(mean(K_wtp(:,idG1,:),1));
        KmeanG2_wtp = squeeze(mean(K_wtp(:,idG2,:),1));

        if graphmode == 1
        fig = figure;
        hold on
        histg1 = histogram(KmeanG1_wtp);
        histg1.Normalization = 'pdf';
        %histg1.BinWidth = 0.05;
        histg1.DisplayName = 'Mean degree G1';
        xmin = floor((min(KmeanG1_wtp,[],'all'))*10)/10;
        xmax = ceil(max(KmeanG1_wtp,[],'all')*10)/10;
        x = xmin:0.05:xmax;
        pd = fitdist(KmeanG1_wtp(:),'Normal');
        y = pdf(pd,x);
        pg1 = plot(x,y,'LineWidth',1);
        pg1.Color = [0 0.4470 0.7410];
        histg2 = histogram(KmeanG2_wtp);
        histg2.Normalization = 'pdf';
        %histg2.BinWidth = 0.05;
        histg2.DisplayName = 'Mean degree G2';
        xmin = floor((min(KmeanG2_wtp,[],'all'))*10)/10;
        xmax = ceil(max(KmeanG2_wtp,[],'all')*10)/10;
        x = xmin:0.05:xmax;
        pd = fitdist(KmeanG2_wtp(:),'Normal');
        y = pdf(pd,x);
        pg2 = plot(x,y,'LineWidth',1);
        pg2.Color = [0.9290 0.6940 0.1250];
        %dim = [.15 .6 .3 .3];
        %str = sprintf('Mean G1 = %.3f\n Std G1 = %.3f\n Mean G2 = %.3f\n Std G2 = %.3f',MmeanchG1p,MstdchG1p,MmeanchG2p,MstdchG2p);
        %annotation('textbox',dim,'String',str,'FitBoxToText','on');
        legend;
        xlabel('Degree Value');
        ylabel('Proportion of the participants (%)');
        title('Histogram of the participants degree values for wtp')
        fig.WindowState = 'maximized';
        savefig([graphpath date '_HistK_wtp']);
        exportgraphics(gcf,[graphpath 'HistK_wtp.png'])
        close
        end

        fprintf('\t Degree_wtp, done\n')

        %%Criteria
        for isubject = 1:size(K_wtp,2)
            goodch = ~ismember(1:size(K_wtp,1), NaNch{1,isubject});
            Kmin_wtp(isubject,:) = min(squeeze(K_wtp(goodch,isubject,:)));
            CCmin_wtp(isubject,:) = min(squeeze(CC_wtp(goodch,isubject,:)));
            try
            Kone_wtp(isubject,:) = find(Kmin_wtp(isubject,:) > 0,1,'first')+1;
            catch
                disp(['Not finding any K>0 for subject ' num2str(isubject)])
               Kone_wtp(isubject,1) = nan;
            end
            try
            CCone_wtp(isubject,:) = find(CCmin_wtp(isubject,:) > 0,1,'first')+1;
            catch
                disp(['Not finding any CC>0 for subject ' num2str(isubject)])
               CCone_wtp(isubject,1) = nan;
            end
        end
        %         for idtr = 1:numel(itra)
        %             Kmin_wtp(:,idtr) = min(K_wtp(:,:,idtr),[],[1 2]);
        %             CCmin_wtp(:,idtr) = min(CC_wtp(:,:,idtr),[],[1 2]);
        %         end

        Kmeanone_wtp = max(Kone_wtp);
        CCmeanone_wtp = max(CCone_wtp);
        Kmeanmin_wtp = min(Kmean_wtp);
        CCmeanmin_wtp = min(CCmean_wtp);
        KminM_wtp = mean(Kmin_wtp);
        CCminM_wtp = mean(CCmin_wtp);
        CCmeanNormallM_wtp = mean(CCmeanNormall_wtp,'omitnan');
        LLNormallM_wtp = mean(LLNormall_wtp,'omitnan');

        if graphmode == 1
            figure
            hold on
            plot(KminM_wtp,'r','displayname','Degree');
            plot(CCminM_wtp+1,'b','displayname','ClustCoeff');
            plot(SWmean_wtp,'g','displayname','SmallWorld');
            ylim([0.5 2])
            yline(1)
            xlabel('Threshold','fontsize',12)
            ylabel('Value','fontsize',12)
            legend
            box on
            %fig.WindowState = 'maximized';
            savefig([graphpath date '_Plotmin_wtp']);
            exportgraphics(gcf,[graphpath 'Plotmin_wtp.png'])
            close
        end
        
        clear b D d1 d1R d2 d2R DRR f id idtr indrand isubject itrv lambda lambdaRR meand1 meand2 tr x y ...
            Ai Ai_wtp RRi_wtp RRmeani_wtp ...
            fig histg1 xmin xmax x y pd pg1  histg2 pg2 dim

        if ~isfile([savepath date '_workspace.mat'])
            save([savepath date '_workspace.mat'],'-v7.3')
        else
            save([savepath date '_workspace.mat'],'-append')
        end
        fprintf('Data saved in %s on %s\n', savepath, date);

        %% AUC
        for ich = 1:size(LE_wtp,1)
            AUCLE_wtp(:,ich) = trapz(squeeze(LE_wtp(ich,:,thrrange)),2);
        end
        AUCLEmean_wtp = trapz(LEmean_wtp(:,thrrange),2);
        for ich = 1:size(LENormallgood_wtp,1)
            AUCLENormall_wtp(:,ich) = trapz(squeeze(LENormallgood_wtp(ich,:,thrrange)),2);
        end
        AUCLEmeanNormall_wtp = trapz(LEmeanNormallgood_wtp(:,thrrange),2);
        for ich = 1:size(ZLEgood_wtp,1)
            AUCZLE_wtp(:,ich) = trapz(squeeze(ZLEgood_wtp(ich,:,thrrange)),2);
        end
        AUCZLEmean_wtp = trapz(ZLEmeangood_wtp(:,thrrange),2);
        AUCGE_wtp = trapz(GE_wtp(:,thrrange),2);
        AUCGENormall_wtp = trapz(GENormallgood_wtp(:,thrrange),2);
        AUCZGE_wtp = trapz(ZGEgood_wtp(:,thrrange),2);
        AUCLL_wtp = trapz(LL_wtp(:,thrrange),2);
        AUCLLNormall_wtp = trapz(LLNormallgood_wtp(:,thrrange),2);
        AUCZLL_wtp = trapz(ZLLgood_wtp(:,thrrange),2);
        for ich = 1:size(CC_wtp,1)
            AUCCC_wtp(:,ich) = trapz(squeeze(CC_wtp(ich,:,thrrange)),2);
        end
        AUCCCmean_wtp = trapz(CCmean_wtp(:,thrrange),2);
        for ich = 1:size(CCNormallgood_wtp,1)
            AUCCCNormall_wtp(:,ich) = trapz(squeeze(CCNormallgood_wtp(ich,:,thrrange)),2);
        end
        AUCCCmeanNormall_wtp = trapz(CCmeanNormallgood_wtp(:,thrrange),2);
        for ich = 1:size(ZCCgood_wtp,1)
            AUCZCC_wtp(:,ich) = trapz(squeeze(ZCCgood_wtp(ich,:,thrrange)),2);
        end
        AUCZCCmean_wtp = trapz(ZCCmeangood_wtp(:,thrrange),2);
        AUCSW_wtp = trapz(SWgood_wtp(:,thrrange),2);
        AUCSWNormall_wtp = trapz(SWNormallgood_wtp(:,thrrange),2);
        AUCZSW_wtp = trapz(ZSWgood_wtp(:,thrrange),2);
        AUCKmean_wtp = trapz(Kmean_wtp(:,thrrange),2);

        clear ich

        if ~isfile([savepath date '_workspace.mat'])
            save([savepath date '_workspace.mat'],'-v7.3')
        else
            save([savepath date '_workspace.mat'],'-append')
        end
        fprintf('Data saved in %s on %s\n', savepath, date);

    end
    %% Proportional threshold & binarized
    if proportionalthresh ==1 && binarizedmode ==1 %BTP
        
        if computegraphmode
            disp('Computing proportional threshold & binarized metrics, running')
    
            LE_btp = zeros(size(MATall,1),size(MATall,3), numel(itrp));
            LERRall_btp = zeros(size(MATall,1), nrand, size(MATall,3), numel(itrp),'single');
            GE_btp = zeros(size(MATall,3), numel(itrp));
            GERRall_btp = zeros(nrand, size(MATall,3), numel(itrp),'single');
            LL_btp = zeros(size(MATall,3), numel(itrp));
            LLRRall_btp = zeros(nrand, size(MATall,3), numel(itrp),'single');
            CC_btp = zeros(size(MATall,1), size(MATall,3), numel(itrp));
            CCRRall_btp = zeros(size(MATall,1), nrand,size(MATall,3), numel(itrp),'single');
            K_btp = zeros(size(MATall,1), size(MATall,3), numel(itrp));
            KRRall_btp = zeros(size(MATall,1), nrand,size(MATall,3), numel(itrp),'single');
            RRi_btp = zeros(size(MATall,1),size(MATall,2),nrand,'single');
            if savemode == 1
                A_btp = zeros(size(MATall,1),size(MATall,2),numel(itrp),size(MATall,3),'single');
                RR_btp = zeros(size(MATall,1),size(MATall,2),nrand,numel(itrp),size(MATall,3),'single'); %RR_btp = zeros(size(MATall,1),size(MATall,2),100);
                D_btp = zeros(size(MATall,1),size(MATall,2),size(MATall,3),numel(itrp),'single');
                DRRall_btp = zeros(size(MATall,1),size(MATall,2),nrand,size(MATall,3),numel(itrp),'single');
            end
    
            for isubject = 1:size(MATall,3)
                fprintf('\t Computing subject %.0f %s\n',isubject, part(isubject));
                fprintf(1,'\t \t Computing threshold:     \n');
                for idtr = 1:numel(itrp)
                    itrv = itrp(idtr);
                    str = num2str(itrv);
                    if length(str) < 4
                        if length(str) == 2
                            str = append(str,'  ');
                        elseif length(str) == 3
                            str = append(str,' ');
                        end
                    end
                    fprintf(1,'\b\b\b\b%s',str);
                    Ai = squeeze(A(:,:,isubject)); %prendre la matrice de chaque participant%
                    Ai_wtp = threshold_proportional(Ai,itrv); %Thresholder la matrice (weigthed)% %Ai_wtp = A_wtp(:,:,idtr,isubject);
                    Ai_btp = weight_conversion(Ai_wtp,'binarize'); %A_btp %Binarizer la matrice (binarized)
                    if savemode == 1
                        A_btp(:,:,idtr,isubject) = Ai_btp;
                        %A_btp(isnan(A_btp)) = 0;  %% Mettre les NAN à 0
                    end
    
                    for indrand = 1:nrand %% For each subject, 100 random networks that preserved the degree distribution of the original network were estimated%
                        if ~any(Ai_btp(:))
                            fprintf(1,'\n \t \t Matrix %s is all zeros. Rdm matrix set as NaN.\n', num2str(indrand))
                            RRi_btp(:,:,:) = NaN;
                            if itrv~=itrp(end)
                                fprintf(1,'\t Computing threshold:     ')
                            end
                            if savemode == 1
                                RR_btp(:,:,:,idtr,isubject) = nan(size(RRi_btp));
                            end
                            break
                        end
                        RRi_btp(:,:,indrand) = randmio_und(Ai_btp,1); %% random matrix %RR_btp(:,:,indrand)
                        if isnan(RRi_btp(:,:,indrand))
                            fprintf(1,'\n \t \t Rdm matrix # %s wont converge. Set as NaN.\n', num2str(indrand))
                            RRi_btp(:,:,:) = NaN;
                            if itrv~=itrp(end)
                                fprintf(1,'\t Computing threshold:     ')
                            end
                            if savemode == 1
                                RR_btp(:,:,:,idtr,isubject) = nan(size(RRi_btp));
                            end
                            break
                        end
                        if savemode == 1
                            RR_btp(:,:,indrand,idtr,isubject) = RRi_btp(:,:,indrand);
                        end
                    end
    
                    RRi_btp(isnan(RRi_btp)) = 0; %Mettre les NAN à 0
                    K_btp(:,isubject, idtr) = degrees_und(Ai_btp);
                    LE_btp(:,isubject, idtr) = efficiency_bin(Ai_btp,1); %calculer le LE pour chaque canal, de chaque matrice binarisée de chaque participant%
                    GE_btp(isubject, idtr) = efficiency_bin(Ai_btp); %calculer le GE pour chaque matrice de chaque participant%
                    D = distance_bin(Ai_btp); %calculer la matrice de distance%
                    [lambda,~,~,~,~] = charpath(D,0,0);
                    if savemode == 1
                        D_btp(:,:,isubject,idtr) = D;
                    end
                    LL_btp(isubject, idtr) = lambda; %calculer le LL pour chaque matrice de distance binarisée de chaque participant%
                    CC_btp(:,isubject, idtr) = clustering_coef_bu(Ai_btp); %calculer le CC pour chaque canal, de chaque matrice de chaque participant%
    
                    for indrand = 1:nrand
                        %RRi_btp(:,:,indrand) = squeeze(RR_btp(:,:,indrand,idtr,isubject));
                        KRRall_btp(:,indrand,isubject,idtr) = degrees_und(RRi_btp(:,:,indrand));
                        LERRall_btp(:,indrand,isubject,idtr) = efficiency_bin(RRi_btp(:,:,indrand),1); %calculer le LE pour chaque canal, de chaque matrice de chaque participant%
                        GERRall_btp(indrand,isubject,idtr) = efficiency_bin(RRi_btp(:,:,indrand)); %calculer le GE pour chaque matrice de chaque participant%
                        DRR = distance_bin(RRi_btp(:,:,indrand));
                        [lambdaRR,~,~,~,~] = charpath(DRR,0,0);
                        if savemode == 1
                            DRRall_btp(:,:,indrand,isubject,idtr) = DRR;
                        end
                        LLRRall_btp(indrand, isubject,idtr) = lambdaRR; %calculer le LL pour chaque matrice de distance de chaque participant%
                        CCRRall_btp(:,indrand,isubject,idtr) = clustering_coef_bu(RRi_btp(:,:,indrand)); %calculer le CC pour chaque canal, de chaque matrice de chaque participant%
                    end
                end
                fprintf(1,'\n')
            end
            clear RRi_btp
    
            if savemode == 1
                if ~isfile([savepath date '_RR.mat'])
                    save([savepath date '_RR.mat'],'RR_btp','-v7.3')
                else
                    save([savepath date '_RR.mat'],'RR_btp','-append')
                end
                if ~isfile([savepath date '_A.mat'])
                    save([savepath date '_A.mat'],'A_btp','-v7.3')
                else
                    save([savepath date '_A.mat'],'A_btp','-append')
                end
                if ~isfile([savepath date '_D.mat'])
                    save([savepath date '_D.mat'],'D_btp','DRRall_btp','-v7.3')
                else
                    save([savepath date '_D.mat'],'D_btp','DRRall_btp','-append')
                end
                clear RR_btp A_btp D_btp DRRall_btp
            end
            disp('Computing proportional threshold & binarized metrics, done')
        end

        %% Metriques
        % Local Efficiency Input Distance matrix, une valeur par canal, par threshold, par participant%
        fprintf('\t Local Efficiency_btp, running\n');

        LEmean_btp = squeeze(mean(LE_btp,1)); %faire la moyenne du LE de tous les canaux, pour chaque matrice de chaque participant%
        LERRmean_btp = squeeze(mean(LERRall_btp,2));
        LERRstd_btp = squeeze(std(LERRall_btp,0,2));
        LEmeanRRall_btp = squeeze(mean(LERRmean_btp,1));

        LENormall_btp = LE_btp./mean(LERRmean_btp,2,'omitnan');
        LENormall_btp(isinf(LENormall_btp)) = NaN;
        LEmeanNormall_btp = squeeze(mean(LENormall_btp,1,'omitnan'));
        LENormallgood_btp = LENormall_btp;
        LENormallgood_btp(isnan(LENormallgood_btp)) = 0;
        LEmeanNormallgood_btp = LEmeanNormall_btp;
        LEmeanNormallgood_btp(isnan(LEmeanNormallgood_btp)) = 0;
        
        ZLE_btp = (LE_btp-LERRmean_btp)./LERRstd_btp;
        ZLE_btp(isinf(ZLE_btp)) = NaN;
        ZLEmean_btp = squeeze(mean(ZLE_btp,1,'omitnan'));
        ZLEgood_btp = ZLE_btp;
        ZLEgood_btp(isnan(ZLEgood_btp)) = 0;
        ZLEmeangood_btp = ZLEmean_btp;
        ZLEmeangood_btp(isnan(ZLEmeangood_btp)) = 0;

        fprintf('\t Local Efficiency_btp, done\n');

        %%Global efficiency, une valeur par threshold par participant
        fprintf('\t Global Efficiency_btp, running\n');

        GERRmean_btp = squeeze(mean(GERRall_btp,1));
        GERRstd_btp = squeeze(std(GERRall_btp,0,1));
        
        GENormall_btp = GE_btp./mean(GERRmean_btp,2,'omitnan');
        GENormall_btp(isinf(GENormall_btp)) = NaN;
        GENormallgood_btp = GENormall_btp;
        GENormallgood_btp (isnan(GENormallgood_btp )) = 0;

        ZGE_btp = (GE_btp-GERRmean_btp)./GERRstd_btp;
        ZGE_btp(isinf(ZGE_btp)) = NaN;
        ZGEgood_btp = ZGE_btp;
        ZGEgood_btp(isnan(ZGEgood_btp)) = 0;

        fprintf('\t Global Efficiency_btp, done\n');

        %%Characteristic pathlenght charpath(D)(BU, BD, WU, WD networks), une valeur par threshold par participant.
        fprintf('\t Characteristic path length_btp, running\n');

        LLRRmean_btp = squeeze(mean(LLRRall_btp,1));
        LLRRstd_btp = squeeze(std(LLRRall_btp,0,1));

        LLNormall_btp = LL_btp./mean(LLRRmean_btp,2,'omitnan');
        LLNormall_btp(isinf(LLNormall_btp)) = NaN;
        LLNormallgood_btp = LLNormall_btp;
        LLNormallgood_btp(isnan(LLNormallgood_btp)) = 0;

        ZLL_btp = (LL_btp-LLRRmean_btp)./LLRRstd_btp;
        ZLL_btp(isinf(ZLL_btp)) = NaN;
        ZLLgood_btp = ZLL_btp;
        ZLLgood_btp(isnan(ZLLgood_btp)) = 0;
 
        fprintf('\t Characteristic path length_btp, done\n');

        %%Clustering Coefficient clustering_coef_bu(A), une valeur par canal, par threshold par participant
        fprintf('\t Clustering coefficient_btp, running\n');

        CCmean_btp = squeeze(mean(CC_btp,1)); %faire la moyenne du CC de tous les canaux, pour chaque matrice binarisée de chaque participant%
        CCRRmean_btp = squeeze(mean(CCRRall_btp,2));
        CCRRstd_btp = squeeze(std(CCRRall_btp,0,2));
        CCmeanRRall_btp = squeeze(mean(CCRRmean_btp,1));
        
        CCNormall_btp = CC_btp./mean(CCRRmean_btp,2,'omitnan');
        CCNormall_btp(isinf(CCNormall_btp)) = NaN;
        CCmeanNormall_btp = squeeze(mean(CCNormall_btp,1,'omitnan'));
        CCNormallgood_btp = CCNormall_btp;
        CCNormallgood_btp(isnan(CCNormallgood_btp)) = 0;
        CCmeanNormallgood_btp = CCmeanNormall_btp;
        CCmeanNormallgood_btp(isnan(CCmeanNormallgood_btp)) = 0;

        ZCC_btp = (CC_btp-CCRRmean_btp)./CCRRstd_btp;
        ZCC_btp(isinf(ZCC_btp)) = NaN;
        ZCCmean_btp = squeeze(mean(ZCC_btp,1,'omitnan'));
        ZCCgood_btp = ZCC_btp;
        ZCCgood_btp(isnan(ZCCgood_btp)) = 0;
        ZCCmeangood_btp = ZCCmean_btp;
        ZCCmeangood_btp(isnan(ZCCmeangood_btp)) = 0;

        fprintf('\t Clustering coefficient_btp, done\n');

        %%Small worldness index%
        fprintf('\t Small Worldness index_btp, running\n');

        SW_btp = CCmean_btp./LL_btp;
        SW_btp(isinf(SW_btp)) = NaN;
        SWgood_btp = SW_btp;
        SWgood_btp(isnan(SWgood_btp)) = 0;

        SWNormall_btp = CCmeanNormall_btp./LLNormall_btp;
        SWNormall_btp(isinf(SWNormall_btp)) = NaN;
        SWNormallgood_btp = SWNormall_btp;
        SWNormallgood_btp(isnan(SWNormallgood_btp)) = 0;
        
        ZSW_btp = ZCCmean_btp./ZLL_btp;
        ZSW_btp(isinf(ZSW_btp)) = NaN;
        ZSWgood_btp = ZSW_btp;
        ZSWgood_btp(isnan(ZSWgood_btp)) = 0;

        SWmean_btp = mean(SWNormall_btp,1,'omitnan');
        SWmeanone_btp = find(SWmean_btp < 1,1,'first')-1;

        fprintf('\t Small Worldness index_btp, done\n');

        %%Degree
        fprintf('\t Degree_btp, running\n')

        Kmean_btp = squeeze(mean(K_btp,1));
        KmeanG1_btp = squeeze(mean(K_btp(:,idG1,:),1));
        KmeanG2_btp = squeeze(mean(K_btp(:,idG2,:),1));

        if graphmode == 1
        fig = figure;
        hold on
        histg1 = histogram(KmeanG1_btp);
        histg1.Normalization = 'pdf';
        %histg1.BinWidth = 0.05;
        histg1.DisplayName = 'Mean degree G1';
        xmin = floor((min(KmeanG1_btp,[],'all'))*10)/10;
        xmax = ceil(max(KmeanG1_btp,[],'all')*10)/10;
        x = xmin:0.05:xmax;
        pd = fitdist(KmeanG1_btp(:),'Normal');
        y = pdf(pd,x);
        pg1 = plot(x,y,'LineWidth',1);
        pg1.Color = [0 0.4470 0.7410];
        histg2 = histogram(KmeanG2_btp);
        histg2.Normalization = 'pdf';
        %histg2.BinWidth = 0.05;
        histg2.DisplayName = 'Mean degree G2';
        xmin = floor((min(KmeanG2_btp,[],'all'))*10)/10;
        xmax = ceil(max(KmeanG2_btp,[],'all')*10)/10;
        x = xmin:0.05:xmax;
        pd = fitdist(KmeanG2_btp(:),'Normal');
        y = pdf(pd,x);
        pg2 = plot(x,y,'LineWidth',1);
        pg2.Color = [0.9290 0.6940 0.1250];
        %dim = [.15 .6 .3 .3];
        %str = sprintf('Mean G1 = %.3f\n Std G1 = %.3f\n Mean G2 = %.3f\n Std G2 = %.3f',MmeanchG1p,MstdchG1p,MmeanchG2p,MstdchG2p);
        %annotation('textbox',dim,'String',str,'FitBoxToText','on');
        legend;
        xlabel('Degree Value');
        ylabel('Proportion of the participants (%)');
        title('Histogram of the participants degree values for btp')
        fig.WindowState = 'maximized';
        savefig([graphpath date '_HistK_btp']);
        exportgraphics(gcf,[graphpath 'HistK_btp.png'])
        close
        end

        fprintf('\t Degree_btp, done\n')

        %%Criteria
        for isubject = 1:size(K_btp,2)
            goodch = ~ismember(1:size(K_btp,1), NaNch{1,isubject});
            Kmin_btp(isubject,:) = min(squeeze(K_btp(goodch,isubject,:)));
            CCmin_btp(isubject,:) = min(squeeze(CC_btp(goodch,isubject,:)));
            try
            Kone_btp(isubject,:) = find(Kmin_btp(isubject,:) > 0,1,'first')+1;
            catch
               disp(['Not finding any K>0 for subject ' num2str(isubject)])
               Kone_btp(isubject,1) = nan;
            end
            try
            CCone_btp(isubject,:) = find(CCmin_btp(isubject,:) > 0,1,'first')+1;
            catch
                disp(['Not finding any CC>0 for subject ' num2str(isubject)])
               CCone_btp(isubject,1) = nan;
            end
        end
        
        Kmeanone_btp = max(Kone_btp);
        CCmeanone_btp = max(CCone_btp);
        Kmeanmin_btp = min(Kmean_btp);
        CCmeanmin_btp = min(CCmean_btp);
        KminM_btp = mean(Kmin_btp);
        CCminM_btp = mean(CCmin_btp);
        CCmeanNormallM_btp = mean(CCmeanNormall_btp,'omitnan');
        LLNormallM_btp = mean(LLNormall_btp,'omitnan');
        
        %         for idtr = 1:numel(itra)
        %             Kmin_btp(:,idtr) = min(K_btp(:,:,idtr),[],[1 2]);
        %             CCmin_btp(:,idtr) = min(CC_btp(:,:,idtr),[],[1 2]);
        %         end

        if graphmode == 1
            figure
            hold on
            plot(KminM_btp,'r','displayname','Degree');
            plot(CCminM_btp+1,'b','displayname','ClustCoeff');
            plot(SWmean_btp,'g','displayname','SmallWorld');
            ylim([0.5 2])
            yline(1)
            xlabel('Threshold','fontsize',12)
            ylabel('Value','fontsize',12)
            legend
            box on
            %fig.WindowState = 'maximized';
            savefig([graphpath date '_Plotmin_btp']);
            exportgraphics(gcf,[graphpath 'Plotmin_btp.png'])
            close
        end

                clear b D d1 d1R d2 d2R DRR f id idtr indrand isubject itrv lambda lambdaRR meand1 meand2 tr x y ...
            Ai Ai_wtp Ai_btp RRi_btp RRmeani_btp ...
            fig histg1 xmin xmax x y pd pg1  histg2 pg2 dim

        %%save%
        if ~isfile([savepath date '_workspace.mat'])
            save([savepath date '_workspace.mat'],'-v7.3')
        else
            save([savepath date '_workspace.mat'],'-append')
        end
        fprintf('Data saved in %s on %s\n', savepath, date);

        %% AUC
        for ich = 1:size(LE_btp,1)
            AUCLE_btp(:,ich) = trapz(squeeze(LE_btp(ich,:,thrrange)),2);
        end
        AUCLEmean_btp = trapz(LEmean_btp(:,thrrange),2);
        for ich = 1:size(LENormallgood_btp,1)
            AUCLENormall_btp(:,ich) = trapz(squeeze(LENormallgood_btp(ich,:,thrrange)),2);
        end
        AUCLEmeanNormall_btp = trapz(LEmeanNormallgood_btp(:,thrrange),2);
        for ich = 1:size(ZLEgood_btp,1)
            AUCZLE_btp(:,ich) = trapz(squeeze(ZLEgood_btp(ich,:,thrrange)),2);
        end
        AUCZLEmean_btp = trapz(ZLEmeangood_btp(:,thrrange),2);
        AUCGE_btp = trapz(GE_btp(:,thrrange),2);
        AUCGENormall_btp = trapz(GENormallgood_btp(:,thrrange),2);
        AUCZGE_btp = trapz(ZGEgood_btp(:,thrrange),2);
        AUCLL_btp = trapz(LL_btp(:,thrrange),2);
        AUCLLNormall_btp = trapz(LLNormallgood_btp(:,thrrange),2);
        AUCZLL_btp = trapz(ZLLgood_btp(:,thrrange),2);
        for ich = 1:size(CC_btp,1)
            AUCCC_btp(:,ich) = trapz(squeeze(CC_btp(ich,:,thrrange)),2);
        end
        AUCCCmean_btp = trapz(CCmean_btp(:,thrrange),2);
        for ich = 1:size(CCNormallgood_btp,1)
            AUCCCNormall_btp(:,ich) = trapz(squeeze(CCNormallgood_btp(ich,:,thrrange)),2);
        end
        AUCCCmeanNormall_btp = trapz(CCmeanNormallgood_btp(:,thrrange),2);
        for ich = 1:size(ZCCgood_btp,1)
            AUCZCC_btp(:,ich) = trapz(squeeze(ZCCgood_btp(ich,:,thrrange)),2);
        end
        AUCZCCmean_btp = trapz(ZCCmeangood_btp(:,thrrange),2);
        AUCSW_btp = trapz(SWgood_btp(:,thrrange),2);
        AUCSWNormall_btp = trapz(SWNormallgood_btp(:,thrrange),2);
        AUCZSW_btp = trapz(ZSWgood_btp(:,thrrange),2);
        AUCKmean_btp = trapz(Kmean_btp(:,thrrange),2);

        %%save%
        if ~isfile([savepath date '_workspace.mat'])
            save([savepath date '_workspace.mat'],'-v7.3')
        else
            save([savepath date '_workspace.mat'],'-append')
        end
        fprintf('Data saved in %s on %s\n', savepath, date);

    end
    %t=toc;
    %disp(datestr(datenum(0,0,0,0,0,t),'HH:MM:SS'))
end

%% GRAPHS
if graphmode == 1

    graphpath = fullfile(savepath,'Graphs', filesep);
    if ~isfolder(graphpath)
        mkdir(graphpath)
    end

    for Meanmode = 0:1
        RRmode = 1;
        label = {'LEmean', {'bta','btp','wta','wtp'}} ;
        Graph_GraphTheory_CO(graphpath, label, itra(thrrange), idG1, idG2, RRmode, Meanmode, LEmean_bta(:,thrrange), LEmeanRRall_bta(:,thrrange), LEmean_btp(:,thrrange), LEmeanRRall_btp(:,thrrange), LEmean_wta(:,thrrange), LEmeanRRall_wta(:,thrrange), LEmean_wtp(:,thrrange), LEmeanRRall_wtp(:,thrrange))
        label = {'GE', {'bta','btp','wta','wtp'}};
        Graph_GraphTheory_CO(graphpath, label, itra(thrrange), idG1, idG2, RRmode, Meanmode, GE_bta(:,thrrange), GERRmean_bta(:,thrrange), GE_btp(:,thrrange), GERRmean_btp(:,thrrange), GE_wta(:,thrrange), GERRmean_wta(:,thrrange), GE_wtp(:,thrrange), GERRmean_wtp(:,thrrange))
        label = {'LL', {'bta','btp','wta','wtp'}};
        Graph_GraphTheory_CO(graphpath, label, itra(thrrange), idG1, idG2, RRmode, Meanmode, LL_bta(:,thrrange), LLRRmean_bta(:,thrrange), LL_btp(:,thrrange), LLRRmean_btp(:,thrrange), LL_wta(:,thrrange), LLRRmean_wta(:,thrrange), LL_wtp(:,thrrange), LLRRmean_wtp(:,thrrange))
        label = {'CCmean', {'bta','btp','wta','wtp'}};
        Graph_GraphTheory_CO(graphpath, label, itra(thrrange), idG1, idG2, RRmode, Meanmode, CCmean_bta(:,thrrange), CCmeanRRall_bta(:,thrrange), CCmean_btp(:,thrrange), CCmeanRRall_btp(:,thrrange), CCmean_wta(:,thrrange), CCmeanRRall_wta(:,thrrange), CCmean_wtp(:,thrrange), CCmeanRRall_wtp(:,thrrange))
        RRmode = 0;
        label = {'LEmeanNorm', {'bta','btp','wta','wtp'}};
        Graph_GraphTheory_CO(graphpath, label, itra(thrrange), idG1, idG2, RRmode, Meanmode, LEmeanNormall_bta(:,thrrange), LEmeanNormall_btp(:,thrrange), LEmeanNormall_wta(:,thrrange), LEmeanNormall_wtp(:,thrrange))
        label = {'GENorm', {'bta','btp','wta','wtp'}};
        Graph_GraphTheory_CO(graphpath, label, itra(thrrange), idG1, idG2, RRmode, Meanmode, GENormall_bta(:,thrrange), GENormall_btp(:,thrrange), GENormall_wta(:,thrrange), GENormall_wtp(:,thrrange))
        label = {'LLNorm', {'bta','btp','wta','wtp'}};
        Graph_GraphTheory_CO(graphpath, label, itra(thrrange), idG1, idG2, RRmode, Meanmode, LLNormall_bta(:,thrrange), LLNormall_btp(:,thrrange), LLNormall_wta(:,thrrange), LLNormall_wtp(:,thrrange))
        label = {'CCmeanNorm', {'bta','btp','wta','wtp'}};
        Graph_GraphTheory_CO(graphpath, label, itra(thrrange), idG1, idG2, RRmode, Meanmode, CCmeanNormall_bta(:,thrrange), CCmeanNormall_btp(:,thrrange), CCmeanNormall_wta(:,thrrange), CCmeanNormall_wtp(:,thrrange))
        label = {'SWNorm', {'bta','btp','wta','wtp'}};
        Graph_GraphTheory_CO(graphpath, label, itra(thrrange), idG1, idG2, RRmode, Meanmode, SWNormall_bta(:,thrrange), SWNormall_btp(:,thrrange), SWNormall_wta(:,thrrange), SWNormall_wtp(:,thrrange))
    end
end

%% ANOVA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if analysemode == 1

    disp(['Computing analyses on ' savepath])

    analysespathi = fullfile(savepath,'Analyses', filesep);
    if ~isfolder(analysespathi)
        mkdir(analysespathi)
    end

    %absolute threshold & weighted
    if absolutethresh ==1 && weightedmode ==1 %WTA
        fprintf('Computing analyses for absolute threshold & weighted metrics, running\n');

        if nointermode
            type = 'main';
            analysespath = fullfile(analysespathi,['ANOVA' type], filesep);
            label = ('wta');
            ANOVA_GraphTheory_CO(LE_wta(:,:,thrrange), LEmean_wta(:,thrrange), GE_wta(:,thrrange),...
                LL_wta(:,thrrange), CC_wta(:,:,thrrange), CCmean_wta(:,thrrange),...
                type, p, itra(thrrange), factors,  labelfactors, label, analysespath, fdrmode, graphmode, SW_wta(:,thrrange),...
                AUCLE_wta, AUCLEmean_wta, AUCGE_wta, AUCLL_wta, AUCCC_wta, AUCCCmean_wta, AUCSW_wta, Kmean_wta(:,thrrange), AUCKmean_wta)

            if normalizemode == 1
                label = ('wtaZ');
                ANOVA_GraphTheory_CO(ZLE_wta(:,:,thrrange), ZLEmean_wta(:,thrrange), ZGE_wta(:,thrrange),...
                    ZLL_wta(:,thrrange), ZCC_wta(:,:,thrrange), ZCCmean_wta(:,thrrange),...
                    type, p, itra(thrrange), factors,  labelfactors, label, analysespath, fdrmode, graphmode, ZSW_wta(:,thrrange), ...
                    AUCZLE_wta, AUCZLEmean_wta, AUCZGE_wta, AUCZLL_wta, AUCZCC_wta, AUCZCCmean_wta, AUCZSW_wta)

                label = ('wtaNormall');
                ANOVA_GraphTheory_CO(LENormall_wta(:,:,thrrange), LEmeanNormall_wta(:,thrrange), GENormall_wta(:,thrrange),...
                    LLNormall_wta(:,thrrange), CCNormall_wta(:,:,thrrange), CCmeanNormall_wta(:,thrrange),...
                    type, p, itra(thrrange), factors,  labelfactors, label, analysespath, fdrmode, graphmode, SWNormall_wta(:,thrrange),...
                    AUCLENormall_wta, AUCLEmeanNormall_wta, AUCGENormall_wta, AUCLLNormall_wta, AUCCCNormall_wta, AUCCCmeanNormall_wta, AUCSWNormall_wta)
            end
        end
        type = 'inter';
        analysespath = fullfile(analysespathi,['ANOVA' type], filesep);
        label = ('wta');
        ANOVA_GraphTheory_CO(LE_wta(:,:,thrrange), LEmean_wta(:,thrrange), GE_wta(:,thrrange),...
            LL_wta(:,thrrange), CC_wta(:,:,thrrange), CCmean_wta(:,thrrange), type,...
            p, itra(thrrange), factors,  labelfactors, label, analysespath, fdrmode, graphmode, SW_wta(:,thrrange),...
            AUCLE_wta, AUCLEmean_wta, AUCGE_wta, AUCLL_wta, AUCCC_wta, AUCCCmean_wta, AUCSW_wta, Kmean_wta(:,thrrange), AUCKmean_wta)

        if normalizemode == 1
            label = ('wtaZ');
            ANOVA_GraphTheory_CO(ZLE_wta(:,:,thrrange), ZLEmean_wta(:,thrrange), ZGE_wta(:,thrrange),...
                ZLL_wta(:,thrrange), ZCC_wta(:,:,thrrange), ZCCmean_wta(:,thrrange),...
                type, p, itra(thrrange), factors,  labelfactors, label, analysespath, fdrmode, graphmode, ZSW_wta(:,thrrange), ...
                AUCZLE_wta, AUCZLEmean_wta, AUCZGE_wta, AUCZLL_wta, AUCZCC_wta, AUCZCCmean_wta, AUCZSW_wta)

            label = ('wtaNormall');
            ANOVA_GraphTheory_CO(LENormall_wta(:,:,thrrange), LEmeanNormall_wta(:,thrrange), GENormall_wta(:,thrrange),...
                LLNormall_wta(:,thrrange), CCNormall_wta(:,:,thrrange), CCmeanNormall_wta(:,thrrange),...
                type, p, itra(thrrange), factors,  labelfactors, label, analysespath, fdrmode, graphmode, SWNormall_wta(:,thrrange),...
                AUCLENormall_wta, AUCLEmeanNormall_wta, AUCGENormall_wta, AUCLLNormall_wta, AUCCCNormall_wta, AUCCCmeanNormall_wta, AUCSWNormall_wta)
        end


        fprintf('Computing analyses for absolute threshold & weighted metrics, done\n\n');
    end

    %%absolute threshold & binarized
    if absolutethresh ==1 && binarizedmode ==1 %BTA
        fprintf('Computing analyses for absolute threshold & binarized metrics, running\n');

        if nointermode
            type = 'main';
            analysespath = fullfile(analysespathi,['ANOVA' type], filesep);
            label = ('bta');
            ANOVA_GraphTheory_CO(LE_bta(:,:,thrrange), LEmean_bta(:,thrrange), GE_bta(:,thrrange),...
                LL_bta(:,thrrange), CC_bta(:,:,thrrange), CCmean_bta(:,thrrange),...
                type, p, itra(thrrange), factors,  labelfactors, label, analysespath, fdrmode, graphmode, SW_bta(:,thrrange),...
                AUCLE_bta, AUCLEmean_bta, AUCGE_bta, AUCLL_bta, AUCCC_bta, AUCCCmean_bta, AUCSW_bta, Kmean_bta(:,thrrange), AUCKmean_bta)

            if normalizemode == 1
                label = ('btaZ');
                ANOVA_GraphTheory_CO(ZLE_bta(:,:,thrrange), ZLEmean_bta(:,thrrange), ZGE_bta(:,thrrange),...
                    ZLL_bta(:,thrrange), ZCC_bta(:,:,thrrange), ZCCmean_bta(:,thrrange),...
                    type, p, itra(thrrange), factors,  labelfactors, label, analysespath, fdrmode, graphmode, ZSW_bta(:,thrrange), ...
                    AUCZLE_bta, AUCZLEmean_bta, AUCZGE_bta, AUCZLL_bta, AUCZCC_bta, AUCZCCmean_bta, AUCZSW_bta)

                label = ('btaNormall');
                ANOVA_GraphTheory_CO(LENormall_bta(:,:,thrrange), LEmeanNormall_bta(:,thrrange), GENormall_bta(:,thrrange),...
                    LLNormall_bta(:,thrrange), CCNormall_bta(:,:,thrrange), CCmeanNormall_bta(:,thrrange),...
                    type, p, itra(thrrange), factors,  labelfactors, label, analysespath, fdrmode, graphmode, SWNormall_bta(:,thrrange),...
                    AUCLENormall_bta, AUCLEmeanNormall_bta, AUCGENormall_bta, AUCLLNormall_bta, AUCCCNormall_bta, AUCCCmeanNormall_bta, AUCSWNormall_bta)
            end
        end
        type = 'inter';
        analysespath = fullfile(analysespathi,['ANOVA' type], filesep);
        label = ('bta');
        ANOVA_GraphTheory_CO(LE_bta(:,:,thrrange), LEmean_bta(:,thrrange), GE_bta(:,thrrange),...
            LL_bta(:,thrrange), CC_bta(:,:,thrrange), CCmean_bta(:,thrrange),...
            type, p, itra(thrrange), factors,  labelfactors, label, analysespath, fdrmode, graphmode, SW_bta(:,thrrange),...
            AUCLE_bta, AUCLEmean_bta, AUCGE_bta, AUCLL_bta, AUCCC_bta, AUCCCmean_bta, AUCSW_bta, Kmean_bta(:,thrrange), AUCKmean_bta)

        if normalizemode == 1
            label = ('btaZ');
            ANOVA_GraphTheory_CO(ZLE_bta(:,:,thrrange), ZLEmean_bta(:,thrrange), ZGE_bta(:,thrrange),...
                ZLL_bta(:,thrrange), ZCC_bta(:,:,thrrange), ZCCmean_bta(:,thrrange),...
                type, p, itra(thrrange), factors,  labelfactors, label, analysespath, fdrmode, graphmode, ZSW_bta(:,thrrange), ...
                AUCZLE_bta, AUCZLEmean_bta, AUCZGE_bta, AUCZLL_bta, AUCZCC_bta, AUCZCCmean_bta, AUCZSW_bta)

            label = ('btaNormall');
            ANOVA_GraphTheory_CO(LENormall_bta(:,:,thrrange), LEmeanNormall_bta(:,thrrange), GENormall_bta(:,thrrange),...
                LLNormall_bta(:,thrrange), CCNormall_bta(:,:,thrrange), CCmeanNormall_bta(:,thrrange),...
                type, p, itra(thrrange), factors,  labelfactors, label, analysespath, fdrmode, graphmode, SWNormall_bta(:,thrrange),...
                AUCLENormall_bta, AUCLEmeanNormall_bta, AUCGENormall_bta, AUCLLNormall_bta, AUCCCNormall_bta, AUCCCmeanNormall_bta, AUCSWNormall_bta)
        end

        fprintf('Computing analyses for absolute threshold & binarized metrics, done\n\n');
    end

    %%proportional threshold & weighted
    if proportionalthresh ==1 && weightedmode ==1 %WTP
        fprintf('Computing analyses for proportional threshold & weighted metrics, running\n');

        if nointermode
            type = 'main';
            analysespath = fullfile(analysespathi,['ANOVA' type], filesep);
            label = ('wtp');
            ANOVA_GraphTheory_CO(LE_wtp(:,:,thrrange), LEmean_wtp(:,thrrange), GE_wtp(:,thrrange),...
                LL_wtp(:,thrrange), CC_wtp(:,:,thrrange), CCmean_wtp(:,thrrange),...
                type, p, itrp(thrrange), factors,  labelfactors, label, analysespath, fdrmode, graphmode, SW_wtp(:,thrrange),...
                AUCLE_wtp, AUCLEmean_wtp, AUCGE_wtp, AUCLL_wtp, AUCCC_wtp, AUCCCmean_wtp, AUCSW_wtp, Kmean_wtp(:,thrrange), AUCKmean_wtp)

            if normalizemode == 1
                label = ('wtpZ');
                ANOVA_GraphTheory_CO(ZLE_wtp(:,:,thrrange), ZLEmean_wtp(:,thrrange), ZGE_wtp(:,thrrange),...
                    ZLL_wtp(:,thrrange), ZCC_wtp(:,:,thrrange), ZCCmean_wtp(:,thrrange),...
                    type, p, itrp(thrrange), factors,  labelfactors, label, analysespath, fdrmode, graphmode, ZSW_wtp(:,thrrange), ...
                    AUCZLE_wtp, AUCZLEmean_wtp, AUCZGE_wtp, AUCZLL_wtp, AUCZCC_wtp, AUCZCCmean_wtp, AUCZSW_wtp)

                label = ('wtpNormall');
                ANOVA_GraphTheory_CO(LENormall_wtp(:,:,thrrange), LEmeanNormall_wtp(:,thrrange), GENormall_wtp(:,thrrange),...
                    LLNormall_wtp(:,thrrange), CCNormall_wtp(:,:,thrrange), CCmeanNormall_wtp(:,thrrange),...
                    type, p, itrp(thrrange), factors,  labelfactors, label, analysespath, fdrmode, graphmode, SWNormall_wtp(:,thrrange),...
                    AUCLENormall_wtp, AUCLEmeanNormall_wtp, AUCGENormall_wtp, AUCLLNormall_wtp, AUCCCNormall_wtp, AUCCCmeanNormall_wtp, AUCSWNormall_wtp)
            end
        end
        type = 'inter';
        analysespath = fullfile(analysespathi,['ANOVA' type], filesep);
        label = ('wtp');
        ANOVA_GraphTheory_CO(LE_wtp(:,:,thrrange), LEmean_wtp(:,thrrange), GE_wtp(:,thrrange),...
            LL_wtp(:,thrrange), CC_wtp(:,:,thrrange), CCmean_wtp(:,thrrange),...
            type, p, itrp(thrrange), factors,  labelfactors, label, analysespath, fdrmode, graphmode, SW_wtp(:,thrrange),...
            AUCLE_wtp, AUCLEmean_wtp, AUCGE_wtp, AUCLL_wtp, AUCCC_wtp, AUCCCmean_wtp, AUCSW_wtp, Kmean_wtp(:,thrrange), AUCKmean_wtp)

        if normalizemode == 1
            label = ('wtpZ');
            ANOVA_GraphTheory_CO(ZLE_wtp(:,:,thrrange), ZLEmean_wtp(:,thrrange), ZGE_wtp(:,thrrange),...
                ZLL_wtp(:,thrrange), ZCC_wtp(:,:,thrrange), ZCCmean_wtp(:,thrrange),...
                type, p, itrp(thrrange), factors,  labelfactors, label, analysespath, fdrmode, graphmode, ZSW_wtp(:,thrrange), ...
                AUCZLE_wtp, AUCZLEmean_wtp, AUCZGE_wtp, AUCZLL_wtp, AUCZCC_wtp, AUCZCCmean_wtp, AUCZSW_wtp)

            label = ('wtpNormall');
            ANOVA_GraphTheory_CO(LENormall_wtp(:,:,thrrange), LEmeanNormall_wtp(:,thrrange), GENormall_wtp(:,thrrange),...
                LLNormall_wtp(:,thrrange), CCNormall_wtp(:,:,thrrange), CCmeanNormall_wtp(:,thrrange),...
                type, p, itrp(thrrange), factors,  labelfactors, label, analysespath, fdrmode, graphmode, SWNormall_wtp(:,thrrange),...
                AUCLENormall_wtp, AUCLEmeanNormall_wtp, AUCGENormall_wtp, AUCLLNormall_wtp, AUCCCNormall_wtp, AUCCCmeanNormall_wtp, AUCSWNormall_wtp)
        end

        fprintf('Computing analyses for proportional threshold & weighted metrics, done\n\n');
    end

    %%proportional threshold & binarized
    if proportionalthresh ==1 && binarizedmode ==1 %BTP
        fprintf('Computing analyses for proportional threshold & binarized metrics, running\n');

        if nointermode
            type = 'main';
            analysespath = fullfile(analysespathi,['ANOVA' type], filesep);
            label = ('btp');
            ANOVA_GraphTheory_CO(LE_btp(:,:,thrrange), LEmean_btp(:,thrrange), GE_btp(:,thrrange),...
                LL_btp(:,thrrange), CC_btp(:,:,thrrange), CCmean_btp(:,thrrange),...
                type, p, itrp(thrrange), factors,  labelfactors, label, analysespath, fdrmode, graphmode, SW_btp(:,thrrange),...
                AUCLE_btp, AUCLEmean_btp, AUCGE_btp, AUCLL_btp, AUCCC_btp, AUCCCmean_btp, AUCSW_btp, Kmean_btp(:,thrrange), AUCKmean_btp)

            if normalizemode == 1
                label = ('btpZ');
                ANOVA_GraphTheory_CO(ZLE_btp(:,:,thrrange), ZLEmean_btp(:,thrrange), ZGE_btp(:,thrrange),...
                    ZLL_btp(:,thrrange), ZCC_btp(:,:,thrrange), ZCCmean_btp(:,thrrange),...
                    type, p, itrp(thrrange), factors,  labelfactors, label, analysespath, fdrmode, graphmode, ZSW_btp(:,thrrange), ...
                    AUCZLE_btp, AUCZLEmean_btp, AUCZGE_btp, AUCZLL_btp, AUCZCC_btp, AUCZCCmean_btp, AUCZSW_btp)

                label = ('btpNormall');
                ANOVA_GraphTheory_CO(LENormall_btp(:,:,thrrange), LEmeanNormall_btp(:,thrrange), GENormall_btp(:,thrrange),...
                    LLNormall_btp(:,thrrange), CCNormall_btp(:,:,thrrange), CCmeanNormall_btp(:,thrrange),...
                    type, p, itrp(thrrange), factors,  labelfactors, label, analysespath, fdrmode, graphmode, SWNormall_btp(:,thrrange),...
                    AUCLENormall_btp, AUCLEmeanNormall_btp, AUCGENormall_btp, AUCLLNormall_btp, AUCCCNormall_btp, AUCCCmeanNormall_btp, AUCSWNormall_btp)
            end
        end

        type = 'inter';
        analysespath = fullfile(analysespathi,['ANOVA' type], filesep);
        label = ('btp');
        ANOVA_GraphTheory_CO(LE_btp(:,:,thrrange), LEmean_btp(:,thrrange), GE_btp(:,thrrange),...
            LL_btp(:,thrrange), CC_btp(:,:,thrrange), CCmean_btp(:,thrrange),...
            type, p, itrp(thrrange), factors,  labelfactors, label, analysespath, fdrmode, graphmode, SW_btp(:,thrrange),...
            AUCLE_btp, AUCLEmean_btp, AUCGE_btp, AUCLL_btp, AUCCC_btp, AUCCCmean_btp, AUCSW_btp, Kmean_btp(:,thrrange), AUCKmean_btp)

        if normalizemode == 1
            label = ('btpZ');
            ANOVA_GraphTheory_CO(ZLE_btp(:,:,thrrange), ZLEmean_btp(:,thrrange), ZGE_btp(:,thrrange),...
                ZLL_btp(:,thrrange), ZCC_btp(:,:,thrrange), ZCCmean_btp(:,thrrange),...
                type, p, itrp(thrrange), factors,  labelfactors, label, analysespath, fdrmode, graphmode, ZSW_btp(:,thrrange), ...
                AUCZLE_btp, AUCZLEmean_btp, AUCZGE_btp, AUCZLL_btp, AUCZCC_btp, AUCZCCmean_btp, AUCZSW_btp)

            label = ('btpNormall');
            ANOVA_GraphTheory_CO(LENormall_btp(:,:,thrrange), LEmeanNormall_btp(:,thrrange), GENormall_btp(:,thrrange),...
                LLNormall_btp(:,thrrange), CCNormall_btp(:,:,thrrange), CCmeanNormall_btp(:,thrrange),...
                type, p, itrp(thrrange), factors,  labelfactors, label, analysespath, fdrmode, graphmode, SWNormall_btp(:,thrrange),...
                AUCLENormall_btp, AUCLEmeanNormall_btp, AUCGENormall_btp, AUCLLNormall_btp, AUCCCNormall_btp, AUCCCmeanNormall_btp, AUCSWNormall_btp)
        end

        fprintf('Computing analyses for proportional threshold & binarized metrics, done\n');
    end
    disp('Computing analyses, done')
end

%% GRAPHS%%%%%%%%
function Graph_GraphTheory_CO(savepath, label, itr, idG1, idG2, RRmode, Meanmode, varargin)

try
if Meanmode
    figure
    t = tiledlayout(2,2,'TileSpacing','Compact');
    if RRmode
        i = 1;
        for n = 1:2:numel(varargin)
            metric = varargin{n};
            metricRR = varargin{n+1};
            nexttile
            hold on
            plot(itr,mean(real(metric(idG1,:)),1,'omitnan'),'Color','#A2142F','displayname','Mal') %'r'
            plot(itr,mean(real(metricRR(idG1,:)),1,'omitnan'),'--','Color','#A2142F','displayname','MalRndm') %'r--'
            plot(itr,mean(real(metric(idG2,:)),1,'omitnan'),'Color','#0072BD','displayname','Ctl') %'b'
            plot(itr,mean(real(metricRR(idG2,:)),1,'omitnan'),'--','Color','#0072BD','displayname','CtlRndm') %'b--'
            xlim([min(itr) max(itr)])
            title([label{1} ' ' label{2}{i}],'fontsize',14)
            hold off
            pbaspect([1 1 1])
            legend
            box on
            i = i + 1;
        end
    else
        for n = 1:numel(varargin)
            metric = varargin{n};
            nexttile
            hold on
            plot(itr,mean(real(metric(idG1,:)),1,'omitnan'),'Color','#A2142F','displayname','Mal')
            plot(itr,mean(real(metric(idG2,:)),1,'omitnan'),'Color','#0072BD','displayname','Ctl')
            xlim([min(itr) max(itr)])
            title([label{1} ' ' label{2}{n}],'fontsize',14)
            hold off
            pbaspect([1 1 1])
            legend
            box on
        end
    end
    xlabel(t,'Threshold','fontsize',14)
    ylabel(t,'Metric value','fontsize',14)
    f = gcf;
    f.Position(3) = 900;
    f.Position(4) = 900;
    movegui(f,'center')
    %f.WindowState = 'maximized';
    savefig([savepath date '_' label{1} '_MEAN.fig'])
    exportgraphics(f,[savepath label{1} '_MEAN.png'])
    close

elseif Meanmode == 0
    figure
    t = tiledlayout(2,2,'TileSpacing','Compact');
    if RRmode
        i = 1;
        for n = 1:2:numel(varargin)
            metric = varargin{n};
            nexttile
            hold on
            plot(itr, metric)
            plot(itr, mean(metric,'omitnan'),'black','lineWidth',1.5)
            xlim([min(itr) max(itr)])
            title([label{1} ' ' label{2}{i}],'fontsize',14)
            hold off
            pbaspect([1 1 1])
            i = i + 1;
        end
        xlabel(t,'Threshold','fontsize',14)
        ylabel(t,'Metric value','fontsize',14)
        f = gcf;
        f.Position(3) = 900;
        f.Position(4) = 900;
        movegui(f,'center')
        %f.WindowState = 'maximized';
        savefig([savepath date '_' label{1} '_ALL.fig'])
        exportgraphics(f,[savepath label{1} '_ALL.png'])
        close

        figure
        t = tiledlayout(2,2,'TileSpacing','Compact');
        i = 1;
        for n = 2:2:numel(varargin)
            metricRR = varargin{n};
            nexttile
            hold on
            plot(itr, metricRR)
            plot(itr, mean(metricRR,'omitnan'),'black','lineWidth',1.5)
            xlim([min(itr) max(itr)])
            title([label{1} 'RR ' label{2}{i}],'fontsize',14)
            hold off
            pbaspect([1 1 1])
            i = i + 1;
        end
        xlabel(t,'Threshold','fontsize',14)
        ylabel(t,'Metric value','fontsize',14)
        f = gcf;
        f.Position(3) = 900;
        f.Position(4) = 900;
        movegui(f,'center')
        %f.WindowState = 'maximized';
        savefig([savepath date '_' label{1} 'RR_ALL.fig'])
        exportgraphics(f,[savepath label{1} 'RR_ALL.png'])
        close
    else
        for n = 1:numel(varargin)
            metric = varargin{n};
            nexttile
            hold on
            plot(itr, metric)
            plot(itr, mean(metric,'omitnan'),'black','lineWidth',1.5)
            xlim([min(itr) max(itr)])
            title([label{1} ' ' label{2}{n}],'fontsize',14)
            hold off
            pbaspect([1 1 1])
        end
        xlabel(t,'Threshold','fontsize',14)
        ylabel(t,'Metric value','fontsize',14)
        f = gcf;
        f.Position(3) = 900;
        f.Position(4) = 900;
        movegui(f,'center')
        %f.WindowState = 'maximized';
        savefig([savepath date '_' label{1} '_ALL.fig'])
        exportgraphics(f,[savepath label{1} '_ALL.png'])
        close
    end
end

catch
    disp('Error during graphic creation')
end
end

%% ANOVA%%%%%%%%%%
function ANOVA_GraphTheory_CO(LE, LEmean, GE, LL, CC, CCmean, type, p, itr, factors, labelfactors, label, savepath, fdrmode, graphmode, varargin)

savemode = 0;

    %try
    if ~isempty(varargin)
        SW = varargin{1};
        AUCLE = varargin{2};
        AUCLEmean = varargin{3};
        AUCGE = varargin{4};
        AUCLL = varargin{5};
        AUCCC = varargin{6};
        AUCCCmean = varargin{7};
        AUCSW = varargin{8};
    end
    if numel(varargin) > 8
        Kmean = varargin{9};
        AUCKmean = varargin{10};
    end
    
    if ~isfolder(savepath)
        mkdir(savepath)
    end
    
    savepathgraph = fullfile(fullfile(savepath,filesep),'Graphs',filesep);
    if ~isfolder(savepathgraph)
        mkdir(savepathgraph)
    end
    
    switch type
        case 'inter'
            fprintf('\t Computing ANOVA with %s and main effects and interactions\n', strjoin(labelfactors,', '))
        case 'main'
            fprintf('\t Computing ANOVA with %s main effects\n', strjoin(labelfactors,', '))
    end
    
    resultsANOVA = struct;
    
    % factors = {gr,ses};
    % labelfactors = {'Group','SES'};
    
    for i = 1:size(LE,1)
        labeldim1{1,i} = ['c' num2str(i)];
    end
    for j = 1:size(LE,3)
        labeldim2{1,j} = ['t' num2str(j)];
    end
    
    labeldata = {'Metric value'; ['LE_' label]};
    [pval, resdata] = ANOVA_job(type, LE, factors, labelfactors, p, labeldata, [{labeldim1};{labeldim2}], fdrmode, graphmode, savepath, savemode);
    resultsANOVA.(labeldata{2}) = {pval; resdata};
    
    if ~isempty(varargin)
        labeldata = {'AUC value'; ['AUCLE_' label]};
        [pval, resdata] = ANOVA_job(type, AUCLE, factors, labelfactors, p, labeldata, labeldim1, fdrmode, graphmode, savepath, savemode);
        resultsANOVA.(labeldata{2}) = {pval; resdata};
    end
    
    labeldata = {'Metric value'; ['LEmean_' label]};
    [pval, resdata] = ANOVA_job(type, LEmean, factors, labelfactors, p, labeldata, labeldim2, fdrmode, graphmode, savepath, savemode);
    resultsANOVA.(labeldata{2}) = {pval; resdata};
    
    if ~isempty(pval)
        f = figure;
        plot(itr, pval{:,:}')
        yline(0.05)
        legend([pval.Properties.RowNames; 'sig'])
        xlabel('Threshold')
        ylabel('Pvalue')
        title(['Pvalue at each threshold for ' strrep(labeldata{2},'_',' ')])
        savefig([savepathgraph 'resultsANOVA' type '_' labeldata{2} '.fig'])
        exportgraphics(f,[savepathgraph 'resultsANOVA' type '_' labeldata{2} '.png'])
        close
    end
    
    if ~isempty(varargin)
        labeldata = {'AUC value'; ['AUCLEmean_' label]};
        [pval, resdata] = ANOVA_job(type, AUCLEmean, factors, labelfactors, p, labeldata, 'NA', fdrmode, graphmode, savepath, savemode);
        resultsANOVA.(labeldata{2}) = {pval; resdata};
    end
    
    labeldata = {'Metric value'; ['GE_' label]};
    [pval, resdata] = ANOVA_job(type, GE, factors, labelfactors, p, labeldata, labeldim2, fdrmode, graphmode, savepath, savemode);
    resultsANOVA.(labeldata{2}) = {pval; resdata};
    
    if ~isempty(pval)
        f = figure;
        plot(itr, pval{:,:}')
        yline(0.05)
        legend([pval.Properties.RowNames; 'sig'])
        xlabel('Threshold')
        ylabel('Pvalue')
        title(['Pvalue at each threshold for  ' strrep(labeldata{2},'_',' ')])
        savefig([savepathgraph 'resultsANOVA' type '_' labeldata{2} '.fig'])
        exportgraphics(f,[savepathgraph 'resultsANOVA' type '_' labeldata{2} '.png'])
        close
    end
    
    if ~isempty(varargin)
        labeldata = {'AUC value'; ['AUCGE_' label]};
        [pval, resdata] = ANOVA_job(type, AUCGE, factors, labelfactors, p, labeldata, 'NA', fdrmode, graphmode, savepath, savemode);
        resultsANOVA.(labeldata{2}) = {pval; resdata};
    end
    
    labeldata = {'Metric value'; ['LL_' label]};
    [pval, resdata] = ANOVA_job(type, LL, factors, labelfactors, p, labeldata, labeldim2, fdrmode, graphmode, savepath, savemode);
    resultsANOVA.(labeldata{2}) = {pval; resdata};
    
    if ~isempty(pval)
        f = figure;
        plot(itr, pval{:,:}')
        yline(0.05)
        legend([pval.Properties.RowNames; 'sig'])
        xlabel('Threshold')
        ylabel('Pvalue')
        title(['Pvalue at each threshold for ' strrep(labeldata{2},'_',' ')])
        savefig([savepathgraph 'resultsANOVA' type '_' labeldata{2} '.fig'])
        exportgraphics(f,[savepathgraph 'resultsANOVA' type '_' labeldata{2} '.png'])
        close
    end
    
    if ~isempty(varargin)
        labeldata = {'AUC value'; ['AUCLL_' label]};
        [pval, resdata] = ANOVA_job(type, AUCLL, factors, labelfactors, p, labeldata, 'NA', fdrmode, graphmode, savepath, savemode);
        resultsANOVA.(labeldata{2}) = {pval; resdata};
    end
    
    labeldata = {'Metric value'; ['CC_' label]};
    [pval, resdata] = ANOVA_job(type, CC, factors, labelfactors, p, labeldata, [{labeldim1};{labeldim2}], fdrmode, graphmode, savepath, savemode);
    resultsANOVA.(labeldata{2}) = {pval; resdata};
    
    if ~isempty(varargin)
        labeldata = {'AUC value'; ['AUCCC_' label]};
        [pval, resdata] = ANOVA_job(type, AUCCC, factors, labelfactors, p, labeldata, labeldim1, fdrmode, graphmode, savepath, savemode);
        resultsANOVA.(labeldata{2}) = {pval; resdata};
    end
    
    labeldata = {'Metric value'; ['CCmean_' label]};
    [pval, resdata] = ANOVA_job(type, CCmean, factors, labelfactors, p, labeldata, labeldim2, fdrmode, graphmode, savepath, savemode);
    resultsANOVA.(labeldata{2}) = {pval; resdata};
    
    if ~isempty(pval)
        f = figure;
        plot(itr, pval{:,:}')
        yline(0.05)
        legend([pval.Properties.RowNames; 'sig'])
        xlabel('Threshold')
        ylabel('Pvalue')
        title(['Pvalue at each threshold for ' strrep(labeldata{2},'_',' ')])
        savefig([savepathgraph 'resultsANOVA' type '_' labeldata{2} '.fig'])
        exportgraphics(f,[savepathgraph 'resultsANOVA' type '_' labeldata{2} '.png'])
        close
    end
    
    if ~isempty(varargin)
        labeldata = {'AUC value'; ['AUCCCmean_' label]};
        [pval, resdata] = ANOVA_job(type, AUCCCmean, factors, labelfactors, p, labeldata, 'NA', fdrmode, graphmode, savepath, savemode);
        resultsANOVA.(labeldata{2}) = {pval; resdata};
    end
    
    if ~isempty(varargin)
        labeldata = {'Metric value'; ['SW_' label]};
        [pval, resdata] = ANOVA_job(type, SW, factors, labelfactors, p, labeldata, labeldim2, fdrmode, graphmode, savepath, savemode);
        resultsANOVA.(labeldata{2}) = {pval; resdata};
    
        if ~isempty(pval)
            f = figure;
            plot(itr, pval{:,:}')
            yline(0.05)
            legend([pval.Properties.RowNames; 'sig'])
            xlabel('Threshold')
            ylabel('Pvalue')
            title(['Pvalue at each threshold for ' strrep(labeldata{2},'_',' ')])
            savefig([savepathgraph 'resultsANOVA' type '_' labeldata{2} '.fig'])
            exportgraphics(f,[savepathgraph 'resultsANOVA' type '_' labeldata{2} '.png'])
            close
        end
    
        labeldata = {'AUC value'; ['AUCSW_' label]};
        [pval, resdata] = ANOVA_job(type, AUCSW, factors, labelfactors, p, labeldata, 'NA', fdrmode, graphmode, savepath, savemode);
        resultsANOVA.(labeldata{2}) = {pval; resdata};
    end
    
    if numel(varargin)>8
        labeldata = {'Metric value'; ['Kmean_' label]};
        [pval, resdata] = ANOVA_job(type, Kmean, factors, labelfactors, p, labeldata, labeldim2, fdrmode, graphmode, savepath, savemode);
        resultsANOVA.(labeldata{2}) = {pval; resdata};
    
        if ~isempty(pval)
            f = figure;
            plot(itr, pval{:,:}')
            yline(0.05)
            legend([pval.Properties.RowNames; 'sig'])
            xlabel('Threshold')
            ylabel('Pvalue')
            title(['Pvalue at each threshold for ' strrep(labeldata{2},'_',' ')])
            savefig([savepathgraph 'resultsANOVA' type '_' labeldata{2} '.fig'])
            exportgraphics(f,[savepathgraph 'resultsANOVA' type '_' labeldata{2} '.png'])
            close
        end
    
        labeldata = {'AUC value'; ['AUCKmean_' label]};
        [pval, resdata] = ANOVA_job(type, AUCKmean, factors, labelfactors, p, labeldata, 'NA', fdrmode, graphmode, savepath, savemode);
        resultsANOVA.(labeldata{2}) = {pval; resdata};
    end
    
    save([fullfile(savepath,filesep) date '_resultsANOVA' type '_' label '.mat'], "resultsANOVA")
    fprintf('\n')
end

%%% RUN For Multiple Datasets (mettre en commentaires le datapath)
% datapath = 'C:\data\Malnutrition\Resting\NIRS\Analyses\Stats\CORR0,01_0,08\PCAPW\ExcludeY\CorrPairC0,25\fisher\';
% run('GraphTheory_CO.m')
% clear
% datapath = 'C:\data\Malnutrition\Resting\NIRS\Analyses\Stats\CORR0,01_0,08\PCAPW\ExcludeY\CorrPairC0\fisher\';
% run('GraphTheory_CO.m')
% clear
% datapath = 'C:\data\Malnutrition\Resting\NIRS\Analyses\Stats\CORR0,01_0,08\PCAPW\ExcludeY\CorrPairC0,5\fisher\';
% run('GraphTheory_CO.m')
% clear