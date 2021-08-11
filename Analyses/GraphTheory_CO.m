%%%%%%%%%%%%%%%%%%%%%%%%%% Graph Theory%%%%%%%%%%%%%%%%%%%%
datapath = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Stats\CORR0,01_0,08\PhysioSatRespEKG\';
%load ([datapath 'workspace.mat'])
load ([datapath 'workspacemat.mat'])

savepath = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Stats\CORR0,01_0,08\PhysioSatRespEKG\GraphTheory\';
if ~isfolder(savepath)
    mkdir(savepath)
end

calculatemode = 0;
analysemode = 1;
normalizemode = 1;
weightedmode = 1;
binarizedmode = 1;
absolutethresh = 1;
proportionalthresh = 1;

%BCT threshold 
%Weigthed = garder les valeurs de corr, Binary = remplacer tout par 0 ou 1.%
itr = 0.01:0.01:1.0;% Seuil variable BCT soustraction
%itrp = 0.01:0.01:0.25;
%p=0.95
fixetr = 0.2; %pour les figures

if calculatemode == 1
    %% Absolute threshold, weighted
    if absolutethresh ==1 && weightedmode ==1 %(WTA)
    %% Local Efficiency Input Distance matrix, une valeur par canal, par threshold, par participant%
    disp('Local Efficiency, running')
    idtr = 1;
    LE_wta = zeros(size(MATall,1),size(MATall,3), numel(itr));
    LERR_wta = zeros(size(MATall,1),size(MATall,3), numel(itr));
    RR_wta = zeros(size(MATall,1),size(MATall,2),100);
    for isubject = 1:size(MATall,3)
        for idtr = 1:numel(itr)
            itrv = itr(idtr);
            A = squeeze(MATall(:,:,isubject)); %prendre la matrice de chaque participant%
            A_wta = threshold_absolute(A, itrv); %Thresholder la matrice (weigthed)%
            A_wta(isnan(A_wta)) = 0;  %% Mettre les NAN à 0
            for indrand = 1:100 %% For each subject, 100 random networks that preserved the degree distribution of the original network were estimated%
                    RR_wta (:,:, indrand) = randomizer_bin_und(A_wta,1); %% random matrix
            end
            RRmean_wta = mean(RR_wta,3); %% average of the 100 estimated random matrices
            LE_wta(:,isubject, idtr)= efficiency_bin(A_wta,1); %calculer le LE pour chaque canal, de chaque matrice binarisée de chaque participant%
            LERR_wta(:,isubject, idtr)= efficiency_bin(RRmean_wta,1); %calculer le LE pour chaque canal, de chaque matrice binarisée de chaque participant%
        end   
    end
    %faire la moyenne du LE de tous les canaux, pour chaque matrice de chaque participant%
    LEmean_wta = squeeze(mean(LE_wta,1));
    %LEmeanG1_wta = nanmean(LEmean_wta(idG1,:),1); %extraire la moyenne de G1
    %LEmeanG2_wta = nanmean(LEmean_wta(idG2,:),1); % extraire la moyenne de G2
    
    %average the random LE of all channels, for each matrix of each participant%.
    LEmeanRR_wta = squeeze(mean(LERR_wta,1));
    %LEmeanRRG1_wta = nanmean(LEmeanRR_wta(idG1,:),1); %extraire la moyenne de G1
    %LEmeanRRG2_wta = nanmean(LEmeanRR_wta(idG2,:),1); % extraire la moyenne de G2

    LENorm_wta = LE_wta./LERR_wta;
    LENormmean_wta = squeeze(mean(LENorm_wta,1));

    disp('Local Efficiency_wta, done')
    save([savepath date '_LE_wta.mat'],'LE_wta','LEmean_wta','LERR_wta','LEmeanRR_wta','LENorm_wta','LENormmean_wta');

    figure
    subplot(1,2,1)
    hold on
    plot(itr,mean(LEmean_wta(idG1,:),1),'b','displayname','mal')
    plot(itr,mean(LEmeanRR_wta(idG1,:),1),'r','displayname','Random mal')
    plot(itr,mean(LEmean_wta(idG2,:),1),'b--','displayname','ctl')
    plot(itr,mean(LEmeanRR_wta(idG2,:),1),'r--','displayname','Random ctl')
    xlabel('Threshold','fontsize',12)
    ylabel('Local Efficiency_wta','fontsize',12)
    legend
    box on

    id = sum(itr< fixetr);
    clear d1 d2 d1R d2R
    d1 = LEmean_wta(idG1,id);
    d2 = LEmean_wta(idG2,id);
    d1R = LEmeanRR_wta(idG1,id);
    d2R = LEmeanRR_wta(idG2,id);

    subplot(1,2,2)
    hold on
    b = bar([mean(d1), mean(d1R); mean(d2), mean(d2R)]);
    b(2).FaceColor = 'r';
    plot(0.85,d1,'x')
    plot(1.15,d1R,'x')
    plot(1.85,d2,'x')
    plot(2.15,d2R,'x')
    title(['Threshold <', num2str(fixetr)])
    xlabel('Groups','fontsize',12)
    ylabel('Local efficiency_wta','fontsize',12)
    xlim([0 , 3])
    legend Real Random
    box on

    savefig([savepath date,'_LocalEfficiency_wta.fig'])
    f=gcf;
    exportgraphics(f,[savepath date '_LocalEfficiency_wta.png'])

    figure
    subplot(1,2,1)
    hold on
    plot(itr,mean(LENormmean_wta(idG1,:),1),'r','displayname','mal')
    plot(itr,mean(LENormmean_wta(idG2,:),1),'b','displayname','ctl')
    xlabel('Threshold','fontsize',12)
    ylabel('Norm Local Efficiency_wta','fontsize',12)

    id = sum(itr< fixetr);
    clear d1 d2
    d1 = LENormmean_wta(idG1,id);
    d2 = LENormmean_wta(idG2,id);

    subplot(2,4,6)
    hold on
    bar(1,mean(d1),'r')
    bar(2,mean(d2),'b')
    plot(1,d1,'x')
    plot(2,d2,'x')
    title(['Threshold <', num2str(fixetr)])
    xlabel('Groups')
    ylabel('Norm Local efficiency_wta')
    xlim([0 , 3])

    savefig([savepath date,'_LocalEfficiency_Norm_wta.fig'])
    f=gcf;
    exportgraphics(f,[savepath date '_LocalEfficiency_Norm_wta.png'])

    %% Global efficiency, une valeur par threshold par participant
    disp('Global Efficiency, running')
    GE_wta = zeros(size(MATall,3), numel(itr)); %%erreur dimension
    GERR_wta = zeros(size(MATall,3), numel(itr));
    RR_wta = zeros(size(MATall,1),size(MATall,2),100);
    idtr = 1;
    for isubject = 1:size(MATall,3) % pour chaque sujet
        for idtr = 1:numel(itr) % pour chaque threshold
            itrv = itr(idtr); %extraire le threshold actuel
            A = squeeze(MATall(:,:,isubject)); %prendre la matrice de chaque participant%
            A_wta = threshold_absolute(A, itrv); %Thresholder la matrice%
            A_wta(isnan(A_wta)) = 0; %% Mettre les NAN à 0
            for indrand = 1:100 %% For each subject, 100 random networks that preserved the degree distribution of the original network were estimated%
                RR_wta (:,:, indrand) = randomizer_bin_und(A_wta,1);
            end
            RRmean_wta = mean(RR_wta,3);
            GE_wta(isubject, idtr) = efficiency_bin(A_wta); %calculer le GE pour chaque matrice de chaque participant%              
            GERR_wta(isubject, idtr) = efficiency_bin(RRmean_wta); %calculer le LE pour chaque canal, de chaque matrice de chaque participant%
        end   
    end
    
    GENorm_wta = GE_wta./GERR_wta;
    
    %average the GE for each matrix of each participant%.
    %GEG1_wta = nanmean(GE_wta(idG1,:),1); %moyenne du G1 pour chaque threshold
    %GEG2_wta = nanmean(GE_wta(idG2,:),1); % moyenne du G2 pour chaque threshold
    %average the random GE for each matrix of each participant%.
    %GERRG1_wta = nanmean(GERR_wta(idG1,:),1); %moyenne du G1 pour chaque threshold
    %GERRG2_wta = nanmean(GERR_wta(idG2,:),1); % moyenne du G2 pour chaque threshold

    disp('Global Efficiency_wta, done')
    save([savepath date '_GE_wta.mat'],'GE_wta', 'GENorm_wta','GERR_wta');

    figure
    subplot(1,2,1);
    hold on;
    plot(itr,mean(GE_wta(idG1,:),1),'b','displayname','mal')
    plot(itr,mean(GERR_wta(idG1,:),1),'r','displayname','Random mal')
    plot(itr,mean(GE_wta(idG2,:),1),'b--','displayname','ctl')
    plot(itr,mean(GERR_wta(idG2,:),1),'r--','displayname','Random ctl')
    xlabel('Threshold','fontsize',12)
    ylabel('Global Efficiency_wta','fontsize',12)
    legend
    box on

    id = sum(itr< fixetr); %trouver la position du threshold sélectionné
    clear d1 d2 d1R d2R
    d1 = GE_wta(idG1,id);
    d2 = GE_wta(idG2,id);
    d1R = GERR_wta(idG1,id);
    d2R = GERR_wta(idG2,id);

    subplot(1,2,2)
    hold on
    b = bar([mean(d1), mean(d1R); mean(d2), mean(d2R)]);
    b(2).FaceColor = 'r';
    plot(0.85,d1,'x')
    plot(1.15,d1R,'x')
    plot(1.85,d2,'x')
    plot(2.15,d2R,'x')
    title(['Threshold <', num2str(fixetr)])
    xlabel('Groups','fontsize',12)
    ylabel('Global Efficiency_wta','fontsize',12)
    xlim([0 , 3])
    legend Real Random
    box on

    savefig([savepath date,'_GlobalEfficiency_wta.fig'])
    f=gcf;
    exportgraphics(f,[savepath date '_GlobalEfficiency_wta.png'])

    figure
    subplot(1,2,1);
    hold on;
    plot(itr,mean(GENorm_wta(idG1,:),1),'r','displayname','mal')
    plot(itr,mean(GENorm_wta(idG2,:),1),'b','displayname','ctl')
    xlabel('Threshold','fontsize',12)
    ylabel('Norm Global Efficiency_wta','fontsize',12)

    id = sum(itr< fixetr); %trouver la position du threshold sélectionné
    d1 = GENorm_wta(idG1,id);
    d2 = GENorm_wta(idG2,id);

    subplot(1,2,2)
    hold on
    bar(1,mean(d1),'r')
    bar(2,mean(d2),'b')
    plot(1,d1,'x')
    plot(2,d2,'x')
    title(['Threshold <', num2str(fixetr)])
    xlabel('Groups')
    ylabel('Norm Global efficiency_wta')
    xlim([0 , 3])

    savefig([savepath date,'_GlobalEfficiency_Norm_wta.fig'])
    f=gcf;
    exportgraphics(f,[savepath date '_GlobalEfficiency_Norm_wta.png'])

    %% Characteristic pathlenght charpath(D)(BU, BD, WU, WD networks), une valeur par threshold par participant.
    disp('Characteristic path length, running')
    idtr = 1;
    LL_wta = zeros(size(MATall,3), numel(itr));
    LLRR_wta = zeros(size(MATall,3), numel(itr));
    RRLL_wta = zeros(size(MATall,1),size(MATall,2),100);
    for isubject = 1:size(MATall,3)
        for idtr = 1:numel(itr)
            itrv = itr(idtr);
            A = squeeze(MATall(:,:,isubject)); %prendre la matrice de chaque participant%
            A_wta = threshold_absolute(A, itrv); %Rendre la matrice binaire en fonction de chaque threshold%
            A_wta(isnan(A_wta)) = 0; %Mettre les NAN à 0
            for indrand = 1:100 %% For each subject, 100 random networks that preserved the degree distribution of the original network were estimated%
                RRLL_wta (:,:, indrand) = randomizer_bin_und(A_wta,1);
            end
            AR_wta = mean(RRLL_wta,3);
            AR_wta(isnan(AR_wta)) = 0; %Mettre les NAN à 0
            RRLLmean_wta = distance_bin(AR_wta);
            D = distance_bin(A_wta); %calculer la matrice de distance%
            [lambda,efficiency,ecc,radius,diameter] = charpath(D,1,0); 
            LL_wta(isubject, idtr) = lambda; %calculer le LL pour chaque matrice de distance binarisée de chaque participant%              
            [lambdaRR,efficiency,ecc,radius,diameter] = charpath(RRLLmean_wta,1,0); 
            LLRR_wta(isubject, idtr) = lambdaRR; %calculer le LL pour chaque matrice de distance binarisée de chaque participant%              
        end   
    end
    %average characteristic pathlenght for each matrix of each participant%.
    %LLG1_wta = nanmean(LL_wta(idG1,:),1);
    %LLG2_wta = nanmean(LL_wta(idG2,:),1);
    %LLmeansubject_wta = nanmean(LL_wta(isubject,:),2);

    %average random characteristic pathlenght for each matrix of each participant%.
    %LLRRG1_wta = nanmean(LLRR_wta(idG1,:),1);
    %LLRRG2_wta = nanmean(LLRR_wta(idG2,:),1);
    %LLRRmeansubject_wta = nanmean(LLRR_wta(isubject,:),2);
    LLNorm_wta = LL_wta./LLRR_wta;
    
    disp('Characteristic path length, done')
    save([savepath date '_LL_wta.mat'],'LL_wta','LLRR_wta','LLNorm_wta');

    figure
    subplot(1,2,1)
    hold on
    plot(itr,mean(LL_wta(idG1,:),1),'b','displayname','mal')
    plot(itr,mean(LLRR_wta(idG1,:),1),'r','displayname','Random mal')
    plot(itr,mean(LL_wta(idG2,:),1),'b--','displayname','ctl')
    plot(itr,mean(LLRR_wta(idG2,:),1),'r--','displayname','Random ctl')
    xlabel('Threshold','fontsize',12)
    ylabel('Characteristic path lenght_wta','fontsize',12)
    legend
    box on

    id = sum(itr< fixetr);
    clear d1 d2 d1R d2R
    d1 = LL_wta(idG1,id);
    d2 = LL_wta(idG2,id);
    d1R = LLRR_wta(idG1,id);
    d2R = LLRR_wta(idG2,id);

    subplot(1,2,2);
    hold on
    b = bar([mean(d1), mean(d1R); mean(d2), mean(d2R)]);
    b(2).FaceColor = 'r';
    plot(0.85,d1,'x')
    plot(1.15,d1R,'x')
    plot(1.85,d2,'x')
    plot(2.15,d2R,'x')
    title(['Threshold <', num2str(fixetr)])
    xlabel('Groups')
    ylabel('Characteristic path lenght_wta')
    xlim([0 , 3])
    legend Real Random
    box on

    savefig([savepath date,'_CharPathLenght_wta.fig'])
    f=gcf;
    exportgraphics(f,[savepath date '_CharPathLenght_wta.png'])
    
    figure
    subplot(1,2,1)
    hold on
    plot(itr,mean(LLNorm_wta(idG1,:),1),'r','displayname','mal')
    plot(itr,mean(LLNorm_wta(idG2,:),1),'b','displayname','ctl')
    xlabel('Threshold','fontsize',12)
    ylabel('Norm Characteristic path lenght_wta','fontsize',12)

    id = sum(itr< fixetr);
    clear d1 d2
    d1 = LLNorm_wta(idG1,id);
    d2 = LLNorm_wta(idG2,id);
    
    subplot(1,2,2)
    hold on
    bar(1,mean(d1),'r')
    bar(2,mean(d2),'b')
    plot(1,d1,'x')
    plot(2,d2,'x')
    xlabel('Groups')
    ylabel('Norm characteristic path length_wta')
    title(['Threshold <', num2str(fixetr)])

    savefig([savepath date,'_CharPathLength_Norm_wta.fig'])
    f=gcf;
    exportgraphics(f,[savepath date '_CharPathLength_Norm_wta.png'])

    %% Clustering Coefficient clustering_coef_bu(A), une valeur par canal, par threshold par participant
    disp('Clustering coefficient, running')
    idtr = 1;
    CC_wta = zeros(size(MATall,1),size(MATall,3), numel(itr));
    CCRR_wta = zeros(size(MATall,1),size(MATall,3), numel(itr));
    RR_wta = zeros(size(MATall,1),size(MATall,2),100);
    for isubject = 1:size(MATall,3)
        for idtr = 1:numel(itr)
            itrv = itr(idtr);
            A = squeeze(MATall(:,:,isubject)); %prendre la matrice de chaque participant%
            A_wta = threshold_absolute(A, itrv); %Rendre la matrice binaire en fonction de chaque threshold%
            A_wta(isnan(A_wta)) = 0;
            for indrand = 1:100 %% For each subject, 100 random networks that preserved the degree distribution of the original network were estimated%
                RR_wta (:,:, indrand) = randomizer_bin_und(A_wta,1);
            end
            RRCCmean_wta = mean(RR_wta,3);
            RRCCmean_wta(isnan(RRmean_wta)) = 0; %Mettre les NAN à 0        
            CC_wta(:,isubject, idtr) = clustering_coef_bu(A_wta); %calculer le CC pour chaque canal, de chaque matrice de chaque participant%
            CCRR_wta(:,isubject, idtr) = clustering_coef_bu(RRmean_wta); %calculer le LE pour chaque canal, de chaque matrice binarisée de chaque participant%
        end
    end

    CCmean_wta = squeeze(mean(CC_wta,1)); %faire la moyenne du CC de tous les canaux, pour chaque matrice binarisée de chaque participant%
    %CCmeanG1 = nanmean(CCmean_wta(idG1,:),1);
    %CCmeanG2 = nanmean(CCmean_wta(idG2,:),1);
    %CCmeansubject = nanmean(CCmean_wta(isubject,:),2);
    
    %%average random CC of all channels, for each binarized matrix of each participant%
    CCmeanRR_wta = squeeze(mean(CCRR_wta,1));
    %CCmeanRRG1 = nanmean(CCmeanRR(idG1,:),1); %extraire la moyenne de G1
    %CCmeanRRG2 = nanmean(CCmeanRR(idG2,:),1); % extraire la moyenne de G2
    %CCmeanRRsubject = nanmean(CCmean_wta(isubject,:),2);
    CCmeanNorm_wta = CCmean_wta./CCmeanRR_wta; 

    disp('Clustering coefficient, done')
    save([savepath date '_CC.mat'],'CC_wta','CCmean_wta','CCmeanRR_wta','CCmeanNorm_wta');

    figure
    subplot(1,2,1)
    hold on
    plot(itr,mean(CCmean_wta(idG1,:),1),'b','displayname','mal')
    plot(itr,mean(CCmeanRR_wta(idG1,:),1),'r','displayname','Random mal')
    plot(itr,mean(CCmean_wta(idG2,:),1),'b--','displayname','ctl')
    plot(itr,mean(CCmeanRR_wta(idG2,:),1),'r--','displayname','Random ctl')
    xlabel('Threshold','fontsize',12)
    ylabel('Clustering Coefficient_wta','fontsize',12)
    legend

    id = sum(itr < fixetr);
    clear d1 d2 d1R d2R
    d1 = CCmean_wta(idG1,id);
    d2 = CCmean_wta(idG2,id);
    d1R = CCmeanRR_wta(idG1,id);
    d2R = CCmeanRR_wta(idG2,id);

    subplot(1,2,2)
    hold on
    b = bar([mean(d1), mean(d1R); mean(d2), mean(d2R)]);
    b(2).FaceColor = 'r';
    plot(0.85,d1,'x')
    plot(1.15,d1R,'x')
    plot(1.85,d2,'x')
    plot(2.15,d2R,'x')
    xlabel('Groups')
    ylabel('Clustering Coefficient_wta')
    title(['Threshold <', num2str(fixetr)])
    legend Real Random
    box on
    xlim([0 , 3])

    savefig([savepath date,'_ClusteringCoeff_wta.fig'])
    f=gcf;
    exportgraphics(f,[savepath date '_ClusteringCoeff_wta.png'])

    figure
    subplot(1,2,1)
    hold on
    plot(itr,mean(CCmeanNorm_wta(idG1,:),1),'r','displayname','mal')
    plot(itr,mean(CCmeanNorm_wta(idG2,:),1),'b','displayname','ctl')
    xlabel('Threshold','fontsize',12)
    ylabel('Norm Clustering Coefficient_wta','fontsize',12)
    legend

    id = sum(itr < fixetr);
    clear d1 d2
    d1 = CCmeanNorm_wta(idG1,id);
    d2 = CCmeanNorm_wta(idG2,id);
    
    subplot(1,2,2)
    hold on
    bar(1,nanmean(d1),'r')
    bar(2,nanmean(d2),'b')
    plot(1,d1,'x')
    plot(2,d2,'x')
    xlabel('Groups')
    ylabel('Norm Clustering Coefficient_wta')
    title(['Threshold <', num2str(fixetr)])
    
    savefig([savepath date,'_ClusteringCoeff_Norm_wta.fig'])
    f=gcf;
    exportgraphics(f,[savepath date '_ClusteringCoeff_Norm_wta.png'])
    
    %savefig([savepath date 'Metrics.fig'])
    %f=gcf;
    %exportgraphics(f,[savepath date 'Metrics.png'])

    %% Small worldness index%
    disp('Small Worldness index, running')

    SW_wta = CCmean_wta./LL_wta;
    %SWG1 = nanmean(SW_wta(idG1,:),1);
    %SWG2 = nanmean(SW_wta(idG2,:),1);

    SWNorm_wta = CCmeanNorm_wta./LLNorm_wta;
    %SWNormG1 = nanmean(SWNorm_wta(idG1,:),1);
    %SWNormG2 = nanmean(SWNorm_wta(idG2,:),1);
    
    figure
    subplot(1,2,1)
    hold on
    plot(itr,squeeze(nanmean(SW_wta(idG1,:),1)),'r','displayname','ctl')
    plot(itr,squeeze(nanmean(SW_wta(idG2,:),1)),'b','displayname','mal')
    xlabel('Threshold','fontsize',16)
    ylabel('Clustering/Path length','fontsize',16)
    title('Small World Index_wta')

    subplot(1,2,2)
    meand1 = mean(SW_wta(idG1,id));
    meand2 = mean(SW_wta(idG2,id));
    x = [1 2];
    y = [meand1 meand2];
    b = bar(x,y,'r');
    b.FaceColor = 'flat';
    b.CData(2,:) = [0 0 1];
    xlabel('Groups')
    ylabel('Small World Index_wta')
    title(['Threshold <', num2str(fixetr)])
    box on
    xlim([0 , 3])

    savefig([savepath date,'_SmallWorld_wta.fig'])
    f=gcf;
    exportgraphics(f,[savepath date 'SmallWorld_wta.png'])

    disp('Small Worldness index, done')
    save([savepath date '_SW_wta.mat'],'SW_wta','SWNorm_wta');

    figure
    subplot(1,2,1)
    hold on
    plot(itr,squeeze(nanmean(SWNorm_wta(idG1,:),1)),'r','displayname','ctl')
    plot(itr,squeeze(nanmean(SWNorm_wta(idG2,:),1)),'b','displayname','mal')
    xlabel('Threshold','fontsize',16)
    ylabel('Norm Clustering/Norm Path length','fontsize',16)
    title('Small World Index_wta')

    subplot(1,2,2)
    meand1 = mean(SWNorm_wta(idG1,id));
    meand2 = mean(SWNorm_wta(idG2,id));
    x = [1 2];
    y = [meand1 meand2];
    b = bar(x,y,'r');
    b.FaceColor = 'flat';
    b.CData(2,:) = [0 0 1];
    xlabel('Groups')
    ylabel('Small World Index Norm_wta')
    title(['Threshold <', num2str(fixetr)])
    box on
    xlim([0 , 3])

    savefig([savepath date,'_SmallWorldNorm_wta.fig'])
    f=gcf;
    exportgraphics(f,[savepath date 'SmallWorldNorm_wta.png'])

    % figure %ajout solène
    % meand1 = mean(d1)
    % meand2 = mean(d2)
    % x = [1 2]
    % y = [meand1 meand2]
    % b = bar(x,y,'r')
    % b.FaceColor = 'flat'
    % b.CData(2,:) = [0 0 1]
    % xlabel('Groups')
    % ylabel('Mean characteristic path lenght')
    % title(['Mean L by groups'])
    % clear y
    % clear x
    %savefig([savepath date,'Mean metrics.fig'])
    %exportgraphics(b,[savepath date 'Mean metrics.png'])
    end
    
    %% Absolute threshold, binarized
    if absolutethresh == 1 && binarizedmode ==1 %BTA
       %% Local Efficiency Input Distance matrix, une valeur par canal, par threshold, par participant%
        disp('Local Efficiency, running')
        idtr = 1;
        LE_bta = zeros(size(MATall,1),size(MATall,3), numel(itr));
        LERR_bta = zeros(size(MATall,1),size(MATall,3), numel(itr));
        RR_bta = zeros(size(MATall,1),size(MATall,2),100);
        for isubject = 1:size(MATall,3)
            for idtr = 1:numel(itr)
                itrv = itr(idtr);
                A = squeeze(MATall(:,:,isubject)); %prendre la matrice de chaque participant%
                A_wta = threshold_absolute(A, itrv); %Thresholder la matrice (weigthed)%
                A_bta = weight_conversion(A_wta,'binarize'); %Binarizer la matrice (binarized)
                A_bta(isnan(A_bta)) = 0;  %% Mettre les NAN à 0
                for indrand = 1:100 %% For each subject, 100 random networks that preserved the degree distribution of the original network were estimated%
                        RR_bta (:,:, indrand) = randomizer_bin_und(A_bta,1); %% random matrix
                end
                RRmean_bta = mean(RR_bta,3); %% average of the 100 estimated random matrices
                LE_bta(:,isubject, idtr)= efficiency_bin(A_bta,1); %calculer le LE pour chaque canal, de chaque matrice binarisée de chaque participant%
                LERR_bta(:,isubject, idtr)= efficiency_bin(RRmean_bta,1); %calculer le LE pour chaque canal, de chaque matrice binarisée de chaque participant%
            end   
        end
        %faire la moyenne du LE de tous les canaux, pour chaque matrice de chaque participant%
        LEmean_bta = squeeze(mean(LE_bta,1));
        %LEmeanG1_bta = nanmean(LEmean_bta(idG1,:),1); %extraire la moyenne de G1
        %LEmeanG2_bta = nanmean(LEmean_bta(idG2,:),1); % extraire la moyenne de G2

        %average the random LE of all channels, for each matrix of each participant%.
        LEmeanRR_bta = squeeze(mean(LERR_bta,1));
        %LEmeanRRG1_bta = nanmean(LEmeanRR_bta(idG1,:),1); %extraire la moyenne de G1
        %LEmeanRRG2_bta = nanmean(LEmeanRR_bta(idG2,:),1); % extraire la moyenne de G2

        LENorm_bta = LE_bta./LERR_bta;
        LENormmean_bta = squeeze(mean(LENorm_bta,1));

        disp('Local Efficiency_bta, done')
        save([savepath date '_LE_bta.mat'],'LE_bta','LEmean_bta','LERR_bta','LEmeanRR_bta','LENorm_bta','LENormmean_bta');

        figure
        subplot(1,2,1)
        hold on
        plot(itr,mean(LEmean_bta(idG1,:),1),'b','displayname','mal')
        plot(itr,mean(LEmeanRR_bta(idG1,:),1),'r','displayname','Random mal')
        plot(itr,mean(LEmean_bta(idG2,:),1),'b--','displayname','ctl')
        plot(itr,mean(LEmeanRR_bta(idG2,:),1),'r--','displayname','Random ctl')
        xlabel('Threshold','fontsize',12)
        ylabel('Local Efficiency_bta','fontsize',12)
        legend
        box on

        id = sum(itr< fixetr);
        clear d1 d2 d1R d2R
        d1 = LEmean_bta(idG1,id);
        d2 = LEmean_bta(idG2,id);
        d1R = LEmeanRR_bta(idG1,id);
        d2R = LEmeanRR_bta(idG2,id);

        subplot(1,2,2)
        hold on
        b = bar([mean(d1), mean(d1R); mean(d2), mean(d2R)]);
        b(2).FaceColor = 'r';
        plot(0.85,d1,'x')
        plot(1.15,d1R,'x')
        plot(1.85,d2,'x')
        plot(2.15,d2R,'x')
        title(['Threshold <', num2str(fixetr)])
        xlabel('Groups','fontsize',12)
        ylabel('Local efficiency_bta','fontsize',12)
        xlim([0 , 3])
        legend Real Random
        box on

        savefig([savepath date,'_LocalEfficiency_bta.fig'])
        f=gcf;
        exportgraphics(f,[savepath date '_LocalEfficiency_bta.png'])

        figure
        subplot(1,2,1)
        hold on
        plot(itr,mean(LENormmean_bta(idG1,:),1),'r','displayname','mal')
        plot(itr,mean(LENormmean_bta(idG2,:),1),'b','displayname','ctl')
        xlabel('Threshold','fontsize',12)
        ylabel('Norm Local Efficiency_bta','fontsize',12)

        id = sum(itr< fixetr);
        clear d1 d2
        d1 = LENormmean_bta(idG1,id);
        d2 = LENormmean_bta(idG2,id);

        subplot(2,4,6)
        hold on
        bar(1,mean(d1),'r')
        bar(2,mean(d2),'b')
        plot(1,d1,'x')
        plot(2,d2,'x')
        title(['Threshold <', num2str(fixetr)])
        xlabel('Groups')
        ylabel('Norm Local efficiency_bta')
        xlim([0 , 3])

        savefig([savepath date,'_LocalEfficiency_Norm_bta.fig'])
        f=gcf;
        exportgraphics(f,[savepath date '_LocalEfficiency_Norm_bta.png'])

        %% Global efficiency, une valeur par threshold par participant
        disp('Global Efficiency, running')
        GE_bta = zeros(size(MATall,3), numel(itr)); %%erreur dimension
        GERR_bta = zeros(size(MATall,3), numel(itr));
        RR_bta = zeros(size(MATall,1),size(MATall,2),100);
        idtr = 1;
        for isubject = 1:size(MATall,3) % pour chaque sujet
            for idtr = 1:numel(itr) % pour chaque threshold
                itrv = itr(idtr); %extraire le threshold actuel
                A = squeeze(MATall(:,:,isubject)); %prendre la matrice de chaque participant%
                A_wta = threshold_absolute(A, itrv); %Thresholder la matrice%
                A_bta = weight_conversion(A_wta,'binarize'); %Binarizer la matrice (binarized)
                A_bta(isnan(A_bta)) = 0; %% Mettre les NAN à 0
                for indrand = 1:100 %% For each subject, 100 random networks that preserved the degree distribution of the original network were estimated%
                    RR_bta (:,:, indrand) = randomizer_bin_und(A_bta,1);
                end
                RRmean_bta = mean(RR_bta,3);
                GE_bta(isubject, idtr) = efficiency_bin(A_bta); %calculer le GE pour chaque matrice de chaque participant%              
                GERR_bta(isubject, idtr) = efficiency_bin(RRmean_bta); %calculer le LE pour chaque canal, de chaque matrice de chaque participant%
            end   
        end

        GENorm_bta = GE_bta./GERR_bta;

        %average the GE for each matrix of each participant%.
        %GEG1_bta = nanmean(GE_bta(idG1,:),1); %moyenne du G1 pour chaque threshold
        %GEG2_bta = nanmean(GE_bta(idG2,:),1); % moyenne du G2 pour chaque threshold
        %average the random GE for each matrix of each participant%.
        %GERRG1_bta = nanmean(GERR_bta(idG1,:),1); %moyenne du G1 pour chaque threshold
        %GERRG2_bta = nanmean(GERR_bta(idG2,:),1); % moyenne du G2 pour chaque threshold

        disp('Global Efficiency_bta, done')
        save([savepath date '_GE_bta.mat'],'GE_bta', 'GENorm_bta','GERR_bta');

        figure
        subplot(1,2,1);
        hold on;
        plot(itr,mean(GE_bta(idG1,:),1),'b','displayname','mal')
        plot(itr,mean(GERR_bta(idG1,:),1),'r','displayname','Random mal')
        plot(itr,mean(GE_bta(idG2,:),1),'b--','displayname','ctl')
        plot(itr,mean(GERR_bta(idG2,:),1),'r--','displayname','Random ctl')
        xlabel('Threshold','fontsize',12)
        ylabel('Global Efficiency_bta','fontsize',12)
        legend
        box on

        id = sum(itr< fixetr); %trouver la position du threshold sélectionné
        clear d1 d2 d1R d2R
        d1 = GE_bta(idG1,id);
        d2 = GE_bta(idG2,id);
        d1R = GERR_bta(idG1,id);
        d2R = GERR_bta(idG2,id);

        subplot(1,2,2)
        hold on
        b = bar([mean(d1), mean(d1R); mean(d2), mean(d2R)]);
        b(2).FaceColor = 'r';
        plot(0.85,d1,'x')
        plot(1.15,d1R,'x')
        plot(1.85,d2,'x')
        plot(2.15,d2R,'x')
        title(['Threshold <', num2str(fixetr)])
        xlabel('Groups','fontsize',12)
        ylabel('Global Efficiency_bta','fontsize',12)
        xlim([0 , 3])
        legend Real Random
        box on

        savefig([savepath date,'_GlobalEfficiency_bta.fig'])
        f=gcf;
        exportgraphics(f,[savepath date '_GlobalEfficiency_bta.png'])

        figure
        subplot(1,2,1);
        hold on;
        plot(itr,mean(GENorm_bta(idG1,:),1),'r','displayname','mal')
        plot(itr,mean(GENorm_bta(idG2,:),1),'b','displayname','ctl')
        xlabel('Threshold','fontsize',12)
        ylabel('Norm Global Efficiency_bta','fontsize',12)

        id = sum(itr< fixetr); %trouver la position du threshold sélectionné
        d1 = GENorm_bta(idG1,id);
        d2 = GENorm_bta(idG2,id);

        subplot(1,2,2)
        hold on
        bar(1,mean(d1),'r')
        bar(2,mean(d2),'b')
        plot(1,d1,'x')
        plot(2,d2,'x')
        title(['Threshold <', num2str(fixetr)])
        xlabel('Groups')
        ylabel('Norm Global efficiency_bta')
        xlim([0 , 3])

        savefig([savepath date,'_GlobalEfficiency_Norm_bta.fig'])
        f=gcf;
        exportgraphics(f,[savepath date '_GlobalEfficiency_Norm_bta.png'])

        %% Characteristic pathlenght charpath(D)(BU, BD, WU, WD networks), une valeur par threshold par participant.
        disp('Characteristic path length, running')
        idtr = 1;
        LL_bta = zeros(size(MATall,3), numel(itr));
        LLRR_bta = zeros(size(MATall,3), numel(itr));
        RRLL_bta = zeros(size(MATall,1),size(MATall,2),100);
        for isubject = 1:size(MATall,3)
            for idtr = 1:numel(itr)
                itrv = itr(idtr);
                A = squeeze(MATall(:,:,isubject)); %prendre la matrice de chaque participant%
                A_wta = threshold_absolute(A, itrv); %Rendre la matrice binaire en fonction de chaque threshold%
                A_bta = weight_conversion(A_wta,'binarize'); %Binarizer la matrice (binarized)
                A_bta(isnan(A_bta)) = 0; %Mettre les NAN à 0
                for indrand = 1:100 %% For each subject, 100 random networks that preserved the degree distribution of the original network were estimated%
                    RRLL_bta (:,:, indrand) = randomizer_bin_und(A_bta,1);
                end
                AR_bta = mean(RRLL_bta,3);
                AR_bta(isnan(AR_bta)) = 0; %Mettre les NAN à 0
                RRLLmean_bta = distance_bin(AR_bta);
                D = distance_bin(A_bta); %calculer la matrice de distance%
                [lambda,efficiency,ecc,radius,diameter] = charpath(D,1,0); 
                LL_bta(isubject, idtr) = lambda; %calculer le LL pour chaque matrice de distance binarisée de chaque participant%              
                [lambdaRR,efficiency,ecc,radius,diameter] = charpath(RRLLmean_bta,1,0); 
                LLRR_bta(isubject, idtr) = lambdaRR; %calculer le LL pour chaque matrice de distance binarisée de chaque participant%              
            end   
        end
        %average characteristic pathlenght for each matrix of each participant%.
        %LLG1_bta = nanmean(LL_bta(idG1,:),1);
        %LLG2_bta = nanmean(LL_bta(idG2,:),1);
        %LLmeansubject_bta = nanmean(LL_bta(isubject,:),2);

        %average random characteristic pathlenght for each matrix of each participant%.
        %LLRRG1_bta = nanmean(LLRR_bta(idG1,:),1);
        %LLRRG2_bta = nanmean(LLRR_bta(idG2,:),1);
        %LLRRmeansubject_bta = nanmean(LLRR_bta(isubject,:),2);
        LLNorm_bta = LL_bta./LLRR_bta;

        disp('Characteristic path length, done')
        save([savepath date '_LL_bta.mat'],'LL_bta','LLRR_bta','LLNorm_bta');

        figure
        subplot(1,2,1)
        hold on
        plot(itr,mean(LL_bta(idG1,:),1),'b','displayname','mal')
        plot(itr,mean(LLRR_bta(idG1,:),1),'r','displayname','Random mal')
        plot(itr,mean(LL_bta(idG2,:),1),'b--','displayname','ctl')
        plot(itr,mean(LLRR_bta(idG2,:),1),'r--','displayname','Random ctl')
        xlabel('Threshold','fontsize',12)
        ylabel('Characteristic path lenght_bta','fontsize',12)
        legend
        box on

        id = sum(itr< fixetr);
        clear d1 d2 d1R d2R
        d1 = LL_bta(idG1,id);
        d2 = LL_bta(idG2,id);
        d1R = LLRR_bta(idG1,id);
        d2R = LLRR_bta(idG2,id);

        subplot(1,2,2);
        hold on
        b = bar([mean(d1), mean(d1R); mean(d2), mean(d2R)]);
        b(2).FaceColor = 'r';
        plot(0.85,d1,'x')
        plot(1.15,d1R,'x')
        plot(1.85,d2,'x')
        plot(2.15,d2R,'x')
        title(['Threshold <', num2str(fixetr)])
        xlabel('Groups')
        ylabel('Characteristic path lenght_bta')
        xlim([0 , 3])
        legend Real Random
        box on

        savefig([savepath date,'_CharPathLenght_bta.fig'])
        f=gcf;
        exportgraphics(f,[savepath date '_CharPathLenght_bta.png'])

        figure
        subplot(1,2,1)
        hold on
        plot(itr,mean(LLNorm_bta(idG1,:),1),'r','displayname','mal')
        plot(itr,mean(LLNorm_bta(idG2,:),1),'b','displayname','ctl')
        xlabel('Threshold','fontsize',12)
        ylabel('Norm Characteristic path lenght_bta','fontsize',12)

        id = sum(itr< fixetr);
        clear d1 d2
        d1 = LLNorm_bta(idG1,id);
        d2 = LLNorm_bta(idG2,id);

        subplot(1,2,2)
        hold on
        bar(1,mean(d1),'r')
        bar(2,mean(d2),'b')
        plot(1,d1,'x')
        plot(2,d2,'x')
        xlabel('Groups')
        ylabel('Norm characteristic path length_bta')
        title(['Threshold <', num2str(fixetr)])

        savefig([savepath date,'_CharPathLength_Norm_bta.fig'])
        f=gcf;
        exportgraphics(f,[savepath date '_CharPathLength_Norm_bta.png'])

        %% Clustering Coefficient clustering_coef_bu(A), une valeur par canal, par threshold par participant
        disp('Clustering coefficient, running')
        idtr = 1;
        CC_bta = zeros(size(MATall,1),size(MATall,3), numel(itr));
        CCRR_bta = zeros(size(MATall,1),size(MATall,3), numel(itr));
        RR_bta = zeros(size(MATall,1),size(MATall,2),100);
        for isubject = 1:size(MATall,3)
            for idtr = 1:numel(itr)
                itrv = itr(idtr);
                A = squeeze(MATall(:,:,isubject)); %prendre la matrice de chaque participant%
                A_wta = threshold_absolute(A, itrv); %Rendre la matrice binaire en fonction de chaque threshold%
                A_bta = weight_conversion(A_wta,'binarize'); %Binarizer la matrice (binarized)
                A_bta(isnan(A_bta)) = 0;
                for indrand = 1:100 %% For each subject, 100 random networks that preserved the degree distribution of the original network were estimated%
                    RR_bta (:,:, indrand) = randomizer_bin_und(A_bta,1);
                end
                RRCCmean_bta = mean(RR_bta,3);
                RRCCmean_bta(isnan(RRmean_bta)) = 0; %Mettre les NAN à 0        
                CC_bta(:,isubject, idtr) = clustering_coef_bu(A_bta); %calculer le CC pour chaque canal, de chaque matrice de chaque participant%
                CCRR_bta(:,isubject, idtr) = clustering_coef_bu(RRmean_bta); %calculer le LE pour chaque canal, de chaque matrice binarisée de chaque participant%
            end
        end

        CCmean_bta = squeeze(mean(CC_bta,1)); %faire la moyenne du CC de tous les canaux, pour chaque matrice binarisée de chaque participant%
        %CCmeanG1 = nanmean(CCmean_bta(idG1,:),1);
        %CCmeanG2 = nanmean(CCmean_bta(idG2,:),1);
        %CCmeansubject = nanmean(CCmean_bta(isubject,:),2);

        %%average random CC of all channels, for each binarized matrix of each participant%
        CCmeanRR_bta = squeeze(mean(CCRR_bta,1));
        %CCmeanRRG1 = nanmean(CCmeanRR(idG1,:),1); %extraire la moyenne de G1
        %CCmeanRRG2 = nanmean(CCmeanRR(idG2,:),1); % extraire la moyenne de G2
        %CCmeanRRsubject = nanmean(CCmean_bta(isubject,:),2);
        CCmeanNorm_bta = CCmean_bta./CCmeanRR_bta; 

        disp('Clustering coefficient, done')
        save([savepath date '_CC.mat'],'CC_bta','CCmean_bta','CCmeanRR_bta','CCmeanNorm_bta');

        figure
        subplot(1,2,1)
        hold on
        plot(itr,mean(CCmean_bta(idG1,:),1),'b','displayname','mal')
        plot(itr,mean(CCmeanRR_bta(idG1,:),1),'r','displayname','Random mal')
        plot(itr,mean(CCmean_bta(idG2,:),1),'b--','displayname','ctl')
        plot(itr,mean(CCmeanRR_bta(idG2,:),1),'r--','displayname','Random ctl')
        xlabel('Threshold','fontsize',12)
        ylabel('Clustering Coefficient_bta','fontsize',12)
        legend

        id = sum(itr < fixetr);
        clear d1 d2 d1R d2R
        d1 = CCmean_bta(idG1,id);
        d2 = CCmean_bta(idG2,id);
        d1R = CCmeanRR_bta(idG1,id);
        d2R = CCmeanRR_bta(idG2,id);

        subplot(1,2,2)
        hold on
        b = bar([mean(d1), mean(d1R); mean(d2), mean(d2R)]);
        b(2).FaceColor = 'r';
        plot(0.85,d1,'x')
        plot(1.15,d1R,'x')
        plot(1.85,d2,'x')
        plot(2.15,d2R,'x')
        xlabel('Groups')
        ylabel('Clustering Coefficient_bta')
        title(['Threshold <', num2str(fixetr)])
        legend Real Random
        box on
        xlim([0 , 3])

        savefig([savepath date,'_ClusteringCoeff_bta.fig'])
        f=gcf;
        exportgraphics(f,[savepath date '_ClusteringCoeff_bta.png'])

        figure
        subplot(1,2,1)
        hold on
        plot(itr,mean(CCmeanNorm_bta(idG1,:),1),'r','displayname','mal')
        plot(itr,mean(CCmeanNorm_bta(idG2,:),1),'b','displayname','ctl')
        xlabel('Threshold','fontsize',12)
        ylabel('Norm Clustering Coefficient_bta','fontsize',12)
        legend

        id = sum(itr < fixetr);
        clear d1 d2
        d1 = CCmeanNorm_bta(idG1,id);
        d2 = CCmeanNorm_bta(idG2,id);

        subplot(1,2,2)
        hold on
        bar(1,nanmean(d1),'r')
        bar(2,nanmean(d2),'b')
        plot(1,d1,'x')
        plot(2,d2,'x')
        xlabel('Groups')
        ylabel('Norm Clustering Coefficient_bta')
        title(['Threshold <', num2str(fixetr)])

        savefig([savepath date,'_ClusteringCoeff_Norm_bta.fig'])
        f=gcf;
        exportgraphics(f,[savepath date '_ClusteringCoeff_Norm_bta.png'])

        %savefig([savepath date 'Metrics.fig'])
        %f=gcf;
        %exportgraphics(f,[savepath date 'Metrics.png'])

        %% Small worldness index%
        disp('Small Worldness index, running')

        SW_bta = CCmean_bta./LL_bta;
        %SWG1 = nanmean(SW_bta(idG1,:),1);
        %SWG2 = nanmean(SW_bta(idG2,:),1);

        SWNorm_bta = CCmeanNorm_bta./LLNorm_bta;
        %SWNormG1 = nanmean(SWNorm_bta(idG1,:),1);
        %SWNormG2 = nanmean(SWNorm_bta(idG2,:),1);

        figure
        subplot(1,2,1)
        hold on
        plot(itr,squeeze(nanmean(SW_bta(idG1,:),1)),'r','displayname','ctl')
        plot(itr,squeeze(nanmean(SW_bta(idG2,:),1)),'b','displayname','mal')
        xlabel('Threshold','fontsize',16)
        ylabel('Clustering/Path length','fontsize',16)
        title('Small World Index_bta')

        subplot(1,2,2)
        meand1 = mean(SW_bta(idG1,id));
        meand2 = mean(SW_bta(idG2,id));
        x = [1 2];
        y = [meand1 meand2];
        b = bar(x,y,'r');
        b.FaceColor = 'flat';
        b.CData(2,:) = [0 0 1];
        xlabel('Groups')
        ylabel('Small World Index_bta')
        title(['Threshold <', num2str(fixetr)])
        box on
        xlim([0 , 3])

        savefig([savepath date,'_SmallWorld_bta.fig'])
        f=gcf;
        exportgraphics(f,[savepath date 'SmallWorld_bta.png'])

        disp('Small Worldness index, done')
        save([savepath date '_SW_bta.mat'],'SW_bta','SWNorm_bta');

        figure
        subplot(1,2,1)
        hold on
        plot(itr,squeeze(nanmean(SWNorm_bta(idG1,:),1)),'r','displayname','ctl')
        plot(itr,squeeze(nanmean(SWNorm_bta(idG2,:),1)),'b','displayname','mal')
        xlabel('Threshold','fontsize',16)
        ylabel('Norm Clustering/Norm Path length','fontsize',16)
        title('Small World Index_bta')

        subplot(1,2,2)
        meand1 = mean(SWNorm_bta(idG1,id));
        meand2 = mean(SWNorm_bta(idG2,id));
        x = [1 2];
        y = [meand1 meand2];
        b = bar(x,y,'r');
        b.FaceColor = 'flat';
        b.CData(2,:) = [0 0 1];
        xlabel('Groups')
        ylabel('Small World Index Norm_bta')
        title(['Threshold <', num2str(fixetr)])
        box on
        xlim([0 , 3])

        savefig([savepath date,'_SmallWorldNorm_bta.fig'])
        f=gcf;
        exportgraphics(f,[savepath date 'SmallWorldNorm_bta.png'])
    end
    %% Proportional threshold & weighted
    if proportionalthresh ==1 && weightedmode ==1 %WTP
       %% Local Efficiency Input Distance matrix, une valeur par canal, par threshold, par participant%
        disp('Local Efficiency, running')
        idtr = 1;
        LE_wtp = zeros(size(MATall,1),size(MATall,3), numel(itr));
        LERR_wtp = zeros(size(MATall,1),size(MATall,3), numel(itr));
        RR_wtp = zeros(size(MATall,1),size(MATall,2),100);
        for isubject = 1:size(MATall,3)
            for idtr = 1:numel(itr)
                itrv = itr(idtr);
                A = squeeze(MATall(:,:,isubject)); %prendre la matrice de chaque participant%
                A_wtp = threshold_proportional(A,itrv); %Thresholder la matrice (weigthed)%
                A_wtp(isnan(A_wtp)) = 0;  %% Mettre les NAN à 0
                for indrand = 1:100 %% For each subject, 100 random networks that preserved the degree distribution of the original network were estimated%
                        RR_wtp (:,:, indrand) = randomizer_bin_und(A_wtp,1); %% random matrix
                end
                RRmean_wtp = mean(RR_wtp,3); %% average of the 100 estimated random matrices
                LE_wtp(:,isubject, idtr)= efficiency_bin(A_wtp,1); %calculer le LE pour chaque canal, de chaque matrice binarisée de chaque participant%
                LERR_wtp(:,isubject, idtr)= efficiency_bin(RRmean_wtp,1); %calculer le LE pour chaque canal, de chaque matrice binarisée de chaque participant%
            end   
        end
        %faire la moyenne du LE de tous les canaux, pour chaque matrice de chaque participant%
        LEmean_wtp = squeeze(mean(LE_wtp,1));
        %LEmeanG1_wtp = nanmean(LEmean_wtp(idG1,:),1); %extraire la moyenne de G1
        %LEmeanG2_wtp = nanmean(LEmean_wtp(idG2,:),1); % extraire la moyenne de G2

        %average the random LE of all channels, for each matrix of each participant%.
        LEmeanRR_wtp = squeeze(mean(LERR_wtp,1));
        %LEmeanRRG1_wtp = nanmean(LEmeanRR_wtp(idG1,:),1); %extraire la moyenne de G1
        %LEmeanRRG2_wtp = nanmean(LEmeanRR_wtp(idG2,:),1); % extraire la moyenne de G2

        LENorm_wtp = LE_wtp./LERR_wtp;
        LENormmean_wtp = squeeze(mean(LENorm_wtp,1));

        disp('Local Efficiency_wtp, done')
        save([savepath date '_LE_wtp.mat'],'LE_wtp','LEmean_wtp','LERR_wtp','LEmeanRR_wtp','LENorm_wtp','LENormmean_wtp');

        figure
        subplot(1,2,1)
        hold on
        plot(itr,mean(LEmean_wtp(idG1,:),1),'b','displayname','mal')
        plot(itr,mean(LEmeanRR_wtp(idG1,:),1),'r','displayname','Random mal')
        plot(itr,mean(LEmean_wtp(idG2,:),1),'b--','displayname','ctl')
        plot(itr,mean(LEmeanRR_wtp(idG2,:),1),'r--','displayname','Random ctl')
        xlabel('Threshold','fontsize',12)
        ylabel('Local Efficiency_wtp','fontsize',12)
        legend
        box on

        id = sum(itr< fixetr);
        clear d1 d2 d1R d2R
        d1 = LEmean_wtp(idG1,id);
        d2 = LEmean_wtp(idG2,id);
        d1R = LEmeanRR_wtp(idG1,id);
        d2R = LEmeanRR_wtp(idG2,id);

        subplot(1,2,2)
        hold on
        b = bar([mean(d1), mean(d1R); mean(d2), mean(d2R)]);
        b(2).FaceColor = 'r';
        plot(0.85,d1,'x')
        plot(1.15,d1R,'x')
        plot(1.85,d2,'x')
        plot(2.15,d2R,'x')
        title(['Threshold <', num2str(fixetr)])
        xlabel('Groups','fontsize',12)
        ylabel('Local efficiency_wtp','fontsize',12)
        xlim([0 , 3])
        legend Real Random
        box on

        savefig([savepath date,'_LocalEfficiency_wtp.fig'])
        f=gcf;
        exportgraphics(f,[savepath date '_LocalEfficiency_wtp.png'])

        figure
        subplot(1,2,1)
        hold on
        plot(itr,mean(LENormmean_wtp(idG1,:),1),'r','displayname','mal')
        plot(itr,mean(LENormmean_wtp(idG2,:),1),'b','displayname','ctl')
        xlabel('Threshold','fontsize',12)
        ylabel('Norm Local Efficiency_wtp','fontsize',12)

        id = sum(itr< fixetr);
        clear d1 d2
        d1 = LENormmean_wtp(idG1,id);
        d2 = LENormmean_wtp(idG2,id);

        subplot(2,4,6)
        hold on
        bar(1,mean(d1),'r')
        bar(2,mean(d2),'b')
        plot(1,d1,'x')
        plot(2,d2,'x')
        title(['Threshold <', num2str(fixetr)])
        xlabel('Groups')
        ylabel('Norm Local efficiency_wtp')
        xlim([0 , 3])

        savefig([savepath date,'_LocalEfficiency_Norm_wtp.fig'])
        f=gcf;
        exportgraphics(f,[savepath date '_LocalEfficiency_Norm_wtp.png'])

        %% Global efficiency, une valeur par threshold par participant
        disp('Global Efficiency, running')
        GE_wtp = zeros(size(MATall,3), numel(itr)); %%erreur dimension
        GERR_wtp = zeros(size(MATall,3), numel(itr));
        RR_wtp = zeros(size(MATall,1),size(MATall,2),100);
        idtr = 1;
        for isubject = 1:size(MATall,3) % pour chaque sujet
            for idtr = 1:numel(itr) % pour chaque threshold
                itrv = itr(idtr); %extraire le threshold actuel
                A = squeeze(MATall(:,:,isubject)); %prendre la matrice de chaque participant%
                A_wtp = threshold_proportional(A,itrv); %Thresholder la matrice (weigthed)%
                A_wtp(isnan(A_wtp)) = 0; %% Mettre les NAN à 0
                for indrand = 1:100 %% For each subject, 100 random networks that preserved the degree distribution of the original network were estimated%
                    RR_wtp (:,:, indrand) = randomizer_bin_und(A_wtp,1);
                end
                RRmean_wtp = mean(RR_wtp,3);
                GE_wtp(isubject, idtr) = efficiency_bin(A_wtp); %calculer le GE pour chaque matrice de chaque participant%              
                GERR_wtp(isubject, idtr) = efficiency_bin(RRmean_wtp); %calculer le LE pour chaque canal, de chaque matrice de chaque participant%
            end   
        end

        GENorm_wtp = GE_wtp./GERR_wtp;

        %average the GE for each matrix of each participant%.
        %GEG1_wtp = nanmean(GE_wtp(idG1,:),1); %moyenne du G1 pour chaque threshold
        %GEG2_wtp = nanmean(GE_wtp(idG2,:),1); % moyenne du G2 pour chaque threshold
        %average the random GE for each matrix of each participant%.
        %GERRG1_wtp = nanmean(GERR_wtp(idG1,:),1); %moyenne du G1 pour chaque threshold
        %GERRG2_wtp = nanmean(GERR_wtp(idG2,:),1); % moyenne du G2 pour chaque threshold

        disp('Global Efficiency_wtp, done')
        save([savepath date '_GE_wtp.mat'],'GE_wtp', 'GENorm_wtp','GERR_wtp');

        figure
        subplot(1,2,1);
        hold on;
        plot(itr,mean(GE_wtp(idG1,:),1),'b','displayname','mal')
        plot(itr,mean(GERR_wtp(idG1,:),1),'r','displayname','Random mal')
        plot(itr,mean(GE_wtp(idG2,:),1),'b--','displayname','ctl')
        plot(itr,mean(GERR_wtp(idG2,:),1),'r--','displayname','Random ctl')
        xlabel('Threshold','fontsize',12)
        ylabel('Global Efficiency_wtp','fontsize',12)
        legend
        box on

        id = sum(itr< fixetr); %trouver la position du threshold sélectionné
        clear d1 d2 d1R d2R
        d1 = GE_wtp(idG1,id);
        d2 = GE_wtp(idG2,id);
        d1R = GERR_wtp(idG1,id);
        d2R = GERR_wtp(idG2,id);

        subplot(1,2,2)
        hold on
        b = bar([mean(d1), mean(d1R); mean(d2), mean(d2R)]);
        b(2).FaceColor = 'r';
        plot(0.85,d1,'x')
        plot(1.15,d1R,'x')
        plot(1.85,d2,'x')
        plot(2.15,d2R,'x')
        title(['Threshold <', num2str(fixetr)])
        xlabel('Groups','fontsize',12)
        ylabel('Global Efficiency_wtp','fontsize',12)
        xlim([0 , 3])
        legend Real Random
        box on

        savefig([savepath date,'_GlobalEfficiency_wtp.fig'])
        f=gcf;
        exportgraphics(f,[savepath date '_GlobalEfficiency_wtp.png'])

        figure
        subplot(1,2,1);
        hold on;
        plot(itr,mean(GENorm_wtp(idG1,:),1),'r','displayname','mal')
        plot(itr,mean(GENorm_wtp(idG2,:),1),'b','displayname','ctl')
        xlabel('Threshold','fontsize',12)
        ylabel('Norm Global Efficiency_wtp','fontsize',12)

        id = sum(itr< fixetr); %trouver la position du threshold sélectionné
        d1 = GENorm_wtp(idG1,id);
        d2 = GENorm_wtp(idG2,id);

        subplot(1,2,2)
        hold on
        bar(1,mean(d1),'r')
        bar(2,mean(d2),'b')
        plot(1,d1,'x')
        plot(2,d2,'x')
        title(['Threshold <', num2str(fixetr)])
        xlabel('Groups')
        ylabel('Norm Global efficiency_wtp')
        xlim([0 , 3])

        savefig([savepath date,'_GlobalEfficiency_Norm_wtp.fig'])
        f=gcf;
        exportgraphics(f,[savepath date '_GlobalEfficiency_Norm_wtp.png'])

        %% Characteristic pathlenght charpath(D)(BU, BD, WU, WD networks), une valeur par threshold par participant.
        disp('Characteristic path length, running')
        idtr = 1;
        LL_wtp = zeros(size(MATall,3), numel(itr));
        LLRR_wtp = zeros(size(MATall,3), numel(itr));
        RRLL_wtp = zeros(size(MATall,1),size(MATall,2),100);
        for isubject = 1:size(MATall,3)
            for idtr = 1:numel(itr)
                itrv = itr(idtr);
                A = squeeze(MATall(:,:,isubject)); %prendre la matrice de chaque participant%
                A_wtp = threshold_proportional(A,itrv); %Thresholder la matrice (weigthed)%
                A_wtp(isnan(A_wtp)) = 0; %Mettre les NAN à 0
                for indrand = 1:100 %% For each subject, 100 random networks that preserved the degree distribution of the original network were estimated%
                    RRLL_wtp (:,:, indrand) = randomizer_bin_und(A_wtp,1);
                end
                AR_wtp = mean(RRLL_wtp,3);
                AR_wtp(isnan(AR_wtp)) = 0; %Mettre les NAN à 0
                RRLLmean_wtp = distance_bin(AR_wtp);
                D = distance_bin(A_wtp); %calculer la matrice de distance%
                [lambda,efficiency,ecc,radius,diameter] = charpath(D,1,0); 
                LL_wtp(isubject, idtr) = lambda; %calculer le LL pour chaque matrice de distance binarisée de chaque participant%              
                [lambdaRR,efficiency,ecc,radius,diameter] = charpath(RRLLmean_wtp,1,0); 
                LLRR_wtp(isubject, idtr) = lambdaRR; %calculer le LL pour chaque matrice de distance binarisée de chaque participant%              
            end   
        end
        %average characteristic pathlenght for each matrix of each participant%.
        %LLG1_wtp = nanmean(LL_wtp(idG1,:),1);
        %LLG2_wtp = nanmean(LL_wtp(idG2,:),1);
        %LLmeansubject_wtp = nanmean(LL_wtp(isubject,:),2);

        %average random characteristic pathlenght for each matrix of each participant%.
        %LLRRG1_wtp = nanmean(LLRR_wtp(idG1,:),1);
        %LLRRG2_wtp = nanmean(LLRR_wtp(idG2,:),1);
        %LLRRmeansubject_wtp = nanmean(LLRR_wtp(isubject,:),2);
        LLNorm_wtp = LL_wtp./LLRR_wtp;

        disp('Characteristic path length, done')
        save([savepath date '_LL_wtp.mat'],'LL_wtp','LLRR_wtp','LLNorm_wtp');

        figure
        subplot(1,2,1)
        hold on
        plot(itr,mean(LL_wtp(idG1,:),1),'b','displayname','mal')
        plot(itr,mean(LLRR_wtp(idG1,:),1),'r','displayname','Random mal')
        plot(itr,mean(LL_wtp(idG2,:),1),'b--','displayname','ctl')
        plot(itr,mean(LLRR_wtp(idG2,:),1),'r--','displayname','Random ctl')
        xlabel('Threshold','fontsize',12)
        ylabel('Characteristic path lenght_wtp','fontsize',12)
        legend
        box on

        id = sum(itr< fixetr);
        clear d1 d2 d1R d2R
        d1 = LL_wtp(idG1,id);
        d2 = LL_wtp(idG2,id);
        d1R = LLRR_wtp(idG1,id);
        d2R = LLRR_wtp(idG2,id);

        subplot(1,2,2);
        hold on
        b = bar([mean(d1), mean(d1R); mean(d2), mean(d2R)]);
        b(2).FaceColor = 'r';
        plot(0.85,d1,'x')
        plot(1.15,d1R,'x')
        plot(1.85,d2,'x')
        plot(2.15,d2R,'x')
        title(['Threshold <', num2str(fixetr)])
        xlabel('Groups')
        ylabel('Characteristic path lenght_wtp')
        xlim([0 , 3])
        legend Real Random
        box on

        savefig([savepath date,'_CharPathLenght_wtp.fig'])
        f=gcf;
        exportgraphics(f,[savepath date '_CharPathLenght_wtp.png'])

        figure
        subplot(1,2,1)
        hold on
        plot(itr,mean(LLNorm_wtp(idG1,:),1),'r','displayname','mal')
        plot(itr,mean(LLNorm_wtp(idG2,:),1),'b','displayname','ctl')
        xlabel('Threshold','fontsize',12)
        ylabel('Norm Characteristic path lenght_wtp','fontsize',12)

        id = sum(itr< fixetr);
        clear d1 d2
        d1 = LLNorm_wtp(idG1,id);
        d2 = LLNorm_wtp(idG2,id);

        subplot(1,2,2)
        hold on
        bar(1,mean(d1),'r')
        bar(2,mean(d2),'b')
        plot(1,d1,'x')
        plot(2,d2,'x')
        xlabel('Groups')
        ylabel('Norm characteristic path length_wtp')
        title(['Threshold <', num2str(fixetr)])

        savefig([savepath date,'_CharPathLength_Norm_wtp.fig'])
        f=gcf;
        exportgraphics(f,[savepath date '_CharPathLength_Norm_wtp.png'])

        %% Clustering Coefficient clustering_coef_bu(A), une valeur par canal, par threshold par participant
        disp('Clustering coefficient, running')
        idtr = 1;
        CC_wtp = zeros(size(MATall,1),size(MATall,3), numel(itr));
        CCRR_wtp = zeros(size(MATall,1),size(MATall,3), numel(itr));
        RR_wtp = zeros(size(MATall,1),size(MATall,2),100);
        for isubject = 1:size(MATall,3)
            for idtr = 1:numel(itr)
                itrv = itr(idtr);
                A = squeeze(MATall(:,:,isubject)); %prendre la matrice de chaque participant%
                A_wtp = threshold_proportional(A,itrv); %Thresholder la matrice (weigthed)%
                A_wtp(isnan(A_wtp)) = 0;
                for indrand = 1:100 %% For each subject, 100 random networks that preserved the degree distribution of the original network were estimated%
                    RR_wtp (:,:, indrand) = randomizer_bin_und(A_wtp,1);
                end
                RRCCmean_wtp = mean(RR_wtp,3);
                RRCCmean_wtp(isnan(RRmean_wtp)) = 0; %Mettre les NAN à 0        
                CC_wtp(:,isubject, idtr) = clustering_coef_bu(A_wtp); %calculer le CC pour chaque canal, de chaque matrice de chaque participant%
                CCRR_wtp(:,isubject, idtr) = clustering_coef_bu(RRmean_wtp); %calculer le LE pour chaque canal, de chaque matrice binarisée de chaque participant%
            end
        end

        CCmean_wtp = squeeze(mean(CC_wtp,1)); %faire la moyenne du CC de tous les canaux, pour chaque matrice binarisée de chaque participant%
        %CCmeanG1 = nanmean(CCmean_wtp(idG1,:),1);
        %CCmeanG2 = nanmean(CCmean_wtp(idG2,:),1);
        %CCmeansubject = nanmean(CCmean_wtp(isubject,:),2);

        %%average random CC of all channels, for each binarized matrix of each participant%
        CCmeanRR_wtp = squeeze(mean(CCRR_wtp,1));
        %CCmeanRRG1 = nanmean(CCmeanRR(idG1,:),1); %extraire la moyenne de G1
        %CCmeanRRG2 = nanmean(CCmeanRR(idG2,:),1); % extraire la moyenne de G2
        %CCmeanRRsubject = nanmean(CCmean_wtp(isubject,:),2);
        CCmeanNorm_wtp = CCmean_wtp./CCmeanRR_wtp; 

        disp('Clustering coefficient, done')
        save([savepath date '_CC.mat'],'CC_wtp','CCmean_wtp','CCmeanRR_wtp','CCmeanNorm_wtp');

        figure
        subplot(1,2,1)
        hold on
        plot(itr,mean(CCmean_wtp(idG1,:),1),'b','displayname','mal')
        plot(itr,mean(CCmeanRR_wtp(idG1,:),1),'r','displayname','Random mal')
        plot(itr,mean(CCmean_wtp(idG2,:),1),'b--','displayname','ctl')
        plot(itr,mean(CCmeanRR_wtp(idG2,:),1),'r--','displayname','Random ctl')
        xlabel('Threshold','fontsize',12)
        ylabel('Clustering Coefficient_wtp','fontsize',12)
        legend

        id = sum(itr < fixetr);
        clear d1 d2 d1R d2R
        d1 = CCmean_wtp(idG1,id);
        d2 = CCmean_wtp(idG2,id);
        d1R = CCmeanRR_wtp(idG1,id);
        d2R = CCmeanRR_wtp(idG2,id);

        subplot(1,2,2)
        hold on
        b = bar([mean(d1), mean(d1R); mean(d2), mean(d2R)]);
        b(2).FaceColor = 'r';
        plot(0.85,d1,'x')
        plot(1.15,d1R,'x')
        plot(1.85,d2,'x')
        plot(2.15,d2R,'x')
        xlabel('Groups')
        ylabel('Clustering Coefficient_wtp')
        title(['Threshold <', num2str(fixetr)])
        legend Real Random
        box on
        xlim([0 , 3])

        savefig([savepath date,'_ClusteringCoeff_wtp.fig'])
        f=gcf;
        exportgraphics(f,[savepath date '_ClusteringCoeff_wtp.png'])

        figure
        subplot(1,2,1)
        hold on
        plot(itr,mean(CCmeanNorm_wtp(idG1,:),1),'r','displayname','mal')
        plot(itr,mean(CCmeanNorm_wtp(idG2,:),1),'b','displayname','ctl')
        xlabel('Threshold','fontsize',12)
        ylabel('Norm Clustering Coefficient_wtp','fontsize',12)
        legend

        id = sum(itr < fixetr);
        clear d1 d2
        d1 = CCmeanNorm_wtp(idG1,id);
        d2 = CCmeanNorm_wtp(idG2,id);

        subplot(1,2,2)
        hold on
        bar(1,nanmean(d1),'r')
        bar(2,nanmean(d2),'b')
        plot(1,d1,'x')
        plot(2,d2,'x')
        xlabel('Groups')
        ylabel('Norm Clustering Coefficient_wtp')
        title(['Threshold <', num2str(fixetr)])

        savefig([savepath date,'_ClusteringCoeff_Norm_wtp.fig'])
        f=gcf;
        exportgraphics(f,[savepath date '_ClusteringCoeff_Norm_wtp.png'])

        %savefig([savepath date 'Metrics.fig'])
        %f=gcf;
        %exportgraphics(f,[savepath date 'Metrics.png'])

        %% Small worldness index%
        disp('Small Worldness index, running')

        SW_wtp = CCmean_wtp./LL_wtp;
        %SWG1 = nanmean(SW_wtp(idG1,:),1);
        %SWG2 = nanmean(SW_wtp(idG2,:),1);

        SWNorm_wtp = CCmeanNorm_wtp./LLNorm_wtp;
        %SWNormG1 = nanmean(SWNorm_wtp(idG1,:),1);
        %SWNormG2 = nanmean(SWNorm_wtp(idG2,:),1);

        figure
        subplot(1,2,1)
        hold on
        plot(itr,squeeze(nanmean(SW_wtp(idG1,:),1)),'r','displayname','ctl')
        plot(itr,squeeze(nanmean(SW_wtp(idG2,:),1)),'b','displayname','mal')
        xlabel('Threshold','fontsize',16)
        ylabel('Clustering/Path length','fontsize',16)
        title('Small World Index_wtp')

        subplot(1,2,2)
        meand1 = mean(SW_wtp(idG1,id));
        meand2 = mean(SW_wtp(idG2,id));
        x = [1 2];
        y = [meand1 meand2];
        b = bar(x,y,'r');
        b.FaceColor = 'flat';
        b.CData(2,:) = [0 0 1];
        xlabel('Groups')
        ylabel('Small World Index_wtp')
        title(['Threshold <', num2str(fixetr)])
        box on
        xlim([0 , 3])

        savefig([savepath date,'_SmallWorld_wtp.fig'])
        f=gcf;
        exportgraphics(f,[savepath date 'SmallWorld_wtp.png'])

        disp('Small Worldness index, done')
        save([savepath date '_SW_wtp.mat'],'SW_wtp','SWNorm_wtp');

        figure
        subplot(1,2,1)
        hold on
        plot(itr,squeeze(nanmean(SWNorm_wtp(idG1,:),1)),'r','displayname','ctl')
        plot(itr,squeeze(nanmean(SWNorm_wtp(idG2,:),1)),'b','displayname','mal')
        xlabel('Threshold','fontsize',16)
        ylabel('Norm Clustering/Norm Path length','fontsize',16)
        title('Small World Index_wtp')

        subplot(1,2,2)
        meand1 = mean(SWNorm_wtp(idG1,id));
        meand2 = mean(SWNorm_wtp(idG2,id));
        x = [1 2];
        y = [meand1 meand2];
        b = bar(x,y,'r');
        b.FaceColor = 'flat';
        b.CData(2,:) = [0 0 1];
        xlabel('Groups')
        ylabel('Small World Index Norm_wtp')
        title(['Threshold <', num2str(fixetr)])
        box on
        xlim([0 , 3])

        savefig([savepath date,'_SmallWorldNorm_wtp.fig'])
        f=gcf;
        exportgraphics(f,[savepath date 'SmallWorldNorm_wtp.png'])
    end
    %% Proportional threshold & binarized
    if proportionalthresh ==1 && binarizedmode ==1 %BTP
         %% Local Efficiency Input Distance matrix, une valeur par canal, par threshold, par participant%
        disp('Local Efficiency, running')
        idtr = 1;
        LE_btp = zeros(size(MATall,1),size(MATall,3), numel(itr));
        LERR_btp = zeros(size(MATall,1),size(MATall,3), numel(itr));
        RR_btp = zeros(size(MATall,1),size(MATall,2),100);
        for isubject = 1:size(MATall,3)
            for idtr = 1:numel(itr)
                itrv = itr(idtr);
                A = squeeze(MATall(:,:,isubject)); %prendre la matrice de chaque participant%
                A_wtp = threshold_proportional(A,itrv); %Thresholder la matrice (weigthed)%
                A_btp = weight_conversion(A_wtp,'binarize'); %Binarizer la matrice (binarized)
                A_btp(isnan(A_btp)) = 0;  %% Mettre les NAN à 0
                for indrand = 1:100 %% For each subject, 100 random networks that preserved the degree distribution of the original network were estimated%
                        RR_btp (:,:, indrand) = randomizer_bin_und(A_btp,1); %% random matrix
                end
                RRmean_btp = mean(RR_btp,3); %% average of the 100 estimated random matrices
                LE_btp(:,isubject, idtr)= efficiency_bin(A_btp,1); %calculer le LE pour chaque canal, de chaque matrice binarisée de chaque participant%
                LERR_btp(:,isubject, idtr)= efficiency_bin(RRmean_btp,1); %calculer le LE pour chaque canal, de chaque matrice binarisée de chaque participant%
            end   
        end
        %faire la moyenne du LE de tous les canaux, pour chaque matrice de chaque participant%
        LEmean_btp = squeeze(mean(LE_btp,1));
        %LEmeanG1_btp = nanmean(LEmean_btp(idG1,:),1); %extraire la moyenne de G1
        %LEmeanG2_btp = nanmean(LEmean_btp(idG2,:),1); % extraire la moyenne de G2

        %average the random LE of all channels, for each matrix of each participant%.
        LEmeanRR_btp = squeeze(mean(LERR_btp,1));
        %LEmeanRRG1_btp = nanmean(LEmeanRR_btp(idG1,:),1); %extraire la moyenne de G1
        %LEmeanRRG2_btp = nanmean(LEmeanRR_btp(idG2,:),1); % extraire la moyenne de G2

        LENorm_btp = LE_btp./LERR_btp;
        LENormmean_btp = squeeze(mean(LENorm_btp,1));

        disp('Local Efficiency_btp, done')
        save([savepath date '_LE_btp.mat'],'LE_btp','LEmean_btp','LERR_btp','LEmeanRR_btp','LENorm_btp','LENormmean_btp');

        figure
        subplot(1,2,1)
        hold on
        plot(itr,mean(LEmean_btp(idG1,:),1),'b','displayname','mal')
        plot(itr,mean(LEmeanRR_btp(idG1,:),1),'r','displayname','Random mal')
        plot(itr,mean(LEmean_btp(idG2,:),1),'b--','displayname','ctl')
        plot(itr,mean(LEmeanRR_btp(idG2,:),1),'r--','displayname','Random ctl')
        xlabel('Threshold','fontsize',12)
        ylabel('Local Efficiency_btp','fontsize',12)
        legend
        box on

        id = sum(itr< fixetr);
        clear d1 d2 d1R d2R
        d1 = LEmean_btp(idG1,id);
        d2 = LEmean_btp(idG2,id);
        d1R = LEmeanRR_btp(idG1,id);
        d2R = LEmeanRR_btp(idG2,id);

        subplot(1,2,2)
        hold on
        b = bar([mean(d1), mean(d1R); mean(d2), mean(d2R)]);
        b(2).FaceColor = 'r';
        plot(0.85,d1,'x')
        plot(1.15,d1R,'x')
        plot(1.85,d2,'x')
        plot(2.15,d2R,'x')
        title(['Threshold <', num2str(fixetr)])
        xlabel('Groups','fontsize',12)
        ylabel('Local efficiency_btp','fontsize',12)
        xlim([0 , 3])
        legend Real Random
        box on

        savefig([savepath date,'_LocalEfficiency_btp.fig'])
        f=gcf;
        exportgraphics(f,[savepath date '_LocalEfficiency_btp.png'])

        figure
        subplot(1,2,1)
        hold on
        plot(itr,mean(LENormmean_btp(idG1,:),1),'r','displayname','mal')
        plot(itr,mean(LENormmean_btp(idG2,:),1),'b','displayname','ctl')
        xlabel('Threshold','fontsize',12)
        ylabel('Norm Local Efficiency_btp','fontsize',12)

        id = sum(itr< fixetr);
        clear d1 d2
        d1 = LENormmean_btp(idG1,id);
        d2 = LENormmean_btp(idG2,id);

        subplot(2,4,6)
        hold on
        bar(1,mean(d1),'r')
        bar(2,mean(d2),'b')
        plot(1,d1,'x')
        plot(2,d2,'x')
        title(['Threshold <', num2str(fixetr)])
        xlabel('Groups')
        ylabel('Norm Local efficiency_btp')
        xlim([0 , 3])

        savefig([savepath date,'_LocalEfficiency_Norm_btp.fig'])
        f=gcf;
        exportgraphics(f,[savepath date '_LocalEfficiency_Norm_btp.png'])

        %% Global efficiency, une valeur par threshold par participant
        disp('Global Efficiency, running')
        GE_btp = zeros(size(MATall,3), numel(itr)); %%erreur dimension
        GERR_btp = zeros(size(MATall,3), numel(itr));
        RR_btp = zeros(size(MATall,1),size(MATall,2),100);
        idtr = 1;
        for isubject = 1:size(MATall,3) % pour chaque sujet
            for idtr = 1:numel(itr) % pour chaque threshold
                itrv = itr(idtr); %extraire le threshold actuel
                A = squeeze(MATall(:,:,isubject)); %prendre la matrice de chaque participant%
                A_wtp = threshold_proportional(A,itrv); %Thresholder la matrice (weigthed)%
                A_btp = weight_conversion(A_wtp,'binarize'); %Binarizer la matrice (binarized)
                A_btp(isnan(A_btp)) = 0; %% Mettre les NAN à 0
                for indrand = 1:100 %% For each subject, 100 random networks that preserved the degree distribution of the original network were estimated%
                    RR_btp (:,:, indrand) = randomizer_bin_und(A_btp,1);
                end
                RRmean_btp = mean(RR_btp,3);
                GE_btp(isubject, idtr) = efficiency_bin(A_btp); %calculer le GE pour chaque matrice de chaque participant%              
                GERR_btp(isubject, idtr) = efficiency_bin(RRmean_btp); %calculer le LE pour chaque canal, de chaque matrice de chaque participant%
            end   
        end

        GENorm_btp = GE_btp./GERR_btp;

        %average the GE for each matrix of each participant%.
        %GEG1_btp = nanmean(GE_btp(idG1,:),1); %moyenne du G1 pour chaque threshold
        %GEG2_btp = nanmean(GE_btp(idG2,:),1); % moyenne du G2 pour chaque threshold
        %average the random GE for each matrix of each participant%.
        %GERRG1_btp = nanmean(GERR_btp(idG1,:),1); %moyenne du G1 pour chaque threshold
        %GERRG2_btp = nanmean(GERR_btp(idG2,:),1); % moyenne du G2 pour chaque threshold

        disp('Global Efficiency_btp, done')
        save([savepath date '_GE_btp.mat'],'GE_btp', 'GENorm_btp','GERR_btp');

        figure
        subplot(1,2,1);
        hold on;
        plot(itr,mean(GE_btp(idG1,:),1),'b','displayname','mal')
        plot(itr,mean(GERR_btp(idG1,:),1),'r','displayname','Random mal')
        plot(itr,mean(GE_btp(idG2,:),1),'b--','displayname','ctl')
        plot(itr,mean(GERR_btp(idG2,:),1),'r--','displayname','Random ctl')
        xlabel('Threshold','fontsize',12)
        ylabel('Global Efficiency_btp','fontsize',12)
        legend
        box on

        id = sum(itr< fixetr); %trouver la position du threshold sélectionné
        clear d1 d2 d1R d2R
        d1 = GE_btp(idG1,id);
        d2 = GE_btp(idG2,id);
        d1R = GERR_btp(idG1,id);
        d2R = GERR_btp(idG2,id);

        subplot(1,2,2)
        hold on
        b = bar([mean(d1), mean(d1R); mean(d2), mean(d2R)]);
        b(2).FaceColor = 'r';
        plot(0.85,d1,'x')
        plot(1.15,d1R,'x')
        plot(1.85,d2,'x')
        plot(2.15,d2R,'x')
        title(['Threshold <', num2str(fixetr)])
        xlabel('Groups','fontsize',12)
        ylabel('Global Efficiency_btp','fontsize',12)
        xlim([0 , 3])
        legend Real Random
        box on

        savefig([savepath date,'_GlobalEfficiency_btp.fig'])
        f=gcf;
        exportgraphics(f,[savepath date '_GlobalEfficiency_btp.png'])

        figure
        subplot(1,2,1);
        hold on;
        plot(itr,mean(GENorm_btp(idG1,:),1),'r','displayname','mal')
        plot(itr,mean(GENorm_btp(idG2,:),1),'b','displayname','ctl')
        xlabel('Threshold','fontsize',12)
        ylabel('Norm Global Efficiency_btp','fontsize',12)

        id = sum(itr< fixetr); %trouver la position du threshold sélectionné
        d1 = GENorm_btp(idG1,id);
        d2 = GENorm_btp(idG2,id);

        subplot(1,2,2)
        hold on
        bar(1,mean(d1),'r')
        bar(2,mean(d2),'b')
        plot(1,d1,'x')
        plot(2,d2,'x')
        title(['Threshold <', num2str(fixetr)])
        xlabel('Groups')
        ylabel('Norm Global efficiency_btp')
        xlim([0 , 3])

        savefig([savepath date,'_GlobalEfficiency_Norm_btp.fig'])
        f=gcf;
        exportgraphics(f,[savepath date '_GlobalEfficiency_Norm_btp.png'])

        %% Characteristic pathlenght charpath(D)(BU, BD, WU, WD networks), une valeur par threshold par participant.
        disp('Characteristic path length, running')
        idtr = 1;
        LL_btp = zeros(size(MATall,3), numel(itr));
        LLRR_btp = zeros(size(MATall,3), numel(itr));
        RRLL_btp = zeros(size(MATall,1),size(MATall,2),100);
        for isubject = 1:size(MATall,3)
            for idtr = 1:numel(itr)
                itrv = itr(idtr);
                A = squeeze(MATall(:,:,isubject)); %prendre la matrice de chaque participant%
                A_wtp = threshold_proportional(A,itrv); %Thresholder la matrice (weigthed)%
                A_btp = weight_conversion(A_wtp,'binarize'); %Binarizer la matrice (binarized)
                A_btp(isnan(A_btp)) = 0; %Mettre les NAN à 0
                for indrand = 1:100 %% For each subject, 100 random networks that preserved the degree distribution of the original network were estimated%
                    RRLL_btp (:,:, indrand) = randomizer_bin_und(A_btp,1);
                end
                AR_btp = mean(RRLL_btp,3);
                AR_btp(isnan(AR_btp)) = 0; %Mettre les NAN à 0
                RRLLmean_btp = distance_bin(AR_btp);
                D = distance_bin(A_btp); %calculer la matrice de distance%
                [lambda,efficiency,ecc,radius,diameter] = charpath(D,1,0); 
                LL_btp(isubject, idtr) = lambda; %calculer le LL pour chaque matrice de distance binarisée de chaque participant%              
                [lambdaRR,efficiency,ecc,radius,diameter] = charpath(RRLLmean_btp,1,0); 
                LLRR_btp(isubject, idtr) = lambdaRR; %calculer le LL pour chaque matrice de distance binarisée de chaque participant%              
            end   
        end
        %average characteristic pathlenght for each matrix of each participant%.
        %LLG1_btp = nanmean(LL_btp(idG1,:),1);
        %LLG2_btp = nanmean(LL_btp(idG2,:),1);
        %LLmeansubject_btp = nanmean(LL_btp(isubject,:),2);

        %average random characteristic pathlenght for each matrix of each participant%.
        %LLRRG1_btp = nanmean(LLRR_btp(idG1,:),1);
        %LLRRG2_btp = nanmean(LLRR_btp(idG2,:),1);
        %LLRRmeansubject_btp = nanmean(LLRR_btp(isubject,:),2);
        LLNorm_btp = LL_btp./LLRR_btp;

        disp('Characteristic path length, done')
        save([savepath date '_LL_btp.mat'],'LL_btp','LLRR_btp','LLNorm_btp');

        figure
        subplot(1,2,1)
        hold on
        plot(itr,mean(LL_btp(idG1,:),1),'b','displayname','mal')
        plot(itr,mean(LLRR_btp(idG1,:),1),'r','displayname','Random mal')
        plot(itr,mean(LL_btp(idG2,:),1),'b--','displayname','ctl')
        plot(itr,mean(LLRR_btp(idG2,:),1),'r--','displayname','Random ctl')
        xlabel('Threshold','fontsize',12)
        ylabel('Characteristic path lenght_btp','fontsize',12)
        legend
        box on

        id = sum(itr< fixetr);
        clear d1 d2 d1R d2R
        d1 = LL_btp(idG1,id);
        d2 = LL_btp(idG2,id);
        d1R = LLRR_btp(idG1,id);
        d2R = LLRR_btp(idG2,id);

        subplot(1,2,2);
        hold on
        b = bar([mean(d1), mean(d1R); mean(d2), mean(d2R)]);
        b(2).FaceColor = 'r';
        plot(0.85,d1,'x')
        plot(1.15,d1R,'x')
        plot(1.85,d2,'x')
        plot(2.15,d2R,'x')
        title(['Threshold <', num2str(fixetr)])
        xlabel('Groups')
        ylabel('Characteristic path lenght_btp')
        xlim([0 , 3])
        legend Real Random
        box on

        savefig([savepath date,'_CharPathLenght_btp.fig'])
        f=gcf;
        exportgraphics(f,[savepath date '_CharPathLenght_btp.png'])

        figure
        subplot(1,2,1)
        hold on
        plot(itr,mean(LLNorm_btp(idG1,:),1),'r','displayname','mal')
        plot(itr,mean(LLNorm_btp(idG2,:),1),'b','displayname','ctl')
        xlabel('Threshold','fontsize',12)
        ylabel('Norm Characteristic path lenght_btp','fontsize',12)

        id = sum(itr< fixetr);
        clear d1 d2
        d1 = LLNorm_btp(idG1,id);
        d2 = LLNorm_btp(idG2,id);

        subplot(1,2,2)
        hold on
        bar(1,mean(d1),'r')
        bar(2,mean(d2),'b')
        plot(1,d1,'x')
        plot(2,d2,'x')
        xlabel('Groups')
        ylabel('Norm characteristic path length_btp')
        title(['Threshold <', num2str(fixetr)])

        savefig([savepath date,'_CharPathLength_Norm_btp.fig'])
        f=gcf;
        exportgraphics(f,[savepath date '_CharPathLength_Norm_btp.png'])

        %% Clustering Coefficient clustering_coef_bu(A), une valeur par canal, par threshold par participant
        disp('Clustering coefficient, running')
        idtr = 1;
        CC_btp = zeros(size(MATall,1),size(MATall,3), numel(itr));
        CCRR_btp = zeros(size(MATall,1),size(MATall,3), numel(itr));
        RR_btp = zeros(size(MATall,1),size(MATall,2),100);
        for isubject = 1:size(MATall,3)
            for idtr = 1:numel(itr)
                itrv = itr(idtr);
                A = squeeze(MATall(:,:,isubject)); %prendre la matrice de chaque participant%
                A_wtp = threshold_proportional(A,itrv); %Thresholder la matrice (weigthed)%
                A_btp = weight_conversion(A_wtp,'binarize'); %Binarizer la matrice (binarized)
                A_btp(isnan(A_btp)) = 0;
                for indrand = 1:100 %% For each subject, 100 random networks that preserved the degree distribution of the original network were estimated%
                    RR_btp (:,:, indrand) = randomizer_bin_und(A_btp,1);
                end
                RRCCmean_btp = mean(RR_btp,3);
                RRCCmean_btp(isnan(RRmean_btp)) = 0; %Mettre les NAN à 0        
                CC_btp(:,isubject, idtr) = clustering_coef_bu(A_btp); %calculer le CC pour chaque canal, de chaque matrice de chaque participant%
                CCRR_btp(:,isubject, idtr) = clustering_coef_bu(RRmean_btp); %calculer le LE pour chaque canal, de chaque matrice binarisée de chaque participant%
            end
        end

        CCmean_btp = squeeze(mean(CC_btp,1)); %faire la moyenne du CC de tous les canaux, pour chaque matrice binarisée de chaque participant%
        %CCmeanG1 = nanmean(CCmean_btp(idG1,:),1);
        %CCmeanG2 = nanmean(CCmean_btp(idG2,:),1);
        %CCmeansubject = nanmean(CCmean_btp(isubject,:),2);

        %%average random CC of all channels, for each binarized matrix of each participant%
        CCmeanRR_btp = squeeze(mean(CCRR_btp,1));
        %CCmeanRRG1 = nanmean(CCmeanRR(idG1,:),1); %extraire la moyenne de G1
        %CCmeanRRG2 = nanmean(CCmeanRR(idG2,:),1); % extraire la moyenne de G2
        %CCmeanRRsubject = nanmean(CCmean_btp(isubject,:),2);
        CCmeanNorm_btp = CCmean_btp./CCmeanRR_btp; 

        disp('Clustering coefficient, done')
        save([savepath date '_CC.mat'],'CC_btp','CCmean_btp','CCmeanRR_btp','CCmeanNorm_btp');

        figure
        subplot(1,2,1)
        hold on
        plot(itr,mean(CCmean_btp(idG1,:),1),'b','displayname','mal')
        plot(itr,mean(CCmeanRR_btp(idG1,:),1),'r','displayname','Random mal')
        plot(itr,mean(CCmean_btp(idG2,:),1),'b--','displayname','ctl')
        plot(itr,mean(CCmeanRR_btp(idG2,:),1),'r--','displayname','Random ctl')
        xlabel('Threshold','fontsize',12)
        ylabel('Clustering Coefficient_btp','fontsize',12)
        legend

        id = sum(itr < fixetr);
        clear d1 d2 d1R d2R
        d1 = CCmean_btp(idG1,id);
        d2 = CCmean_btp(idG2,id);
        d1R = CCmeanRR_btp(idG1,id);
        d2R = CCmeanRR_btp(idG2,id);

        subplot(1,2,2)
        hold on
        b = bar([mean(d1), mean(d1R); mean(d2), mean(d2R)]);
        b(2).FaceColor = 'r';
        plot(0.85,d1,'x')
        plot(1.15,d1R,'x')
        plot(1.85,d2,'x')
        plot(2.15,d2R,'x')
        xlabel('Groups')
        ylabel('Clustering Coefficient_btp')
        title(['Threshold <', num2str(fixetr)])
        legend Real Random
        box on
        xlim([0 , 3])

        savefig([savepath date,'_ClusteringCoeff_btp.fig'])
        f=gcf;
        exportgraphics(f,[savepath date '_ClusteringCoeff_btp.png'])

        figure
        subplot(1,2,1)
        hold on
        plot(itr,mean(CCmeanNorm_btp(idG1,:),1),'r','displayname','mal')
        plot(itr,mean(CCmeanNorm_btp(idG2,:),1),'b','displayname','ctl')
        xlabel('Threshold','fontsize',12)
        ylabel('Norm Clustering Coefficient_btp','fontsize',12)
        legend

        id = sum(itr < fixetr);
        clear d1 d2
        d1 = CCmeanNorm_btp(idG1,id);
        d2 = CCmeanNorm_btp(idG2,id);

        subplot(1,2,2)
        hold on
        bar(1,nanmean(d1),'r')
        bar(2,nanmean(d2),'b')
        plot(1,d1,'x')
        plot(2,d2,'x')
        xlabel('Groups')
        ylabel('Norm Clustering Coefficient_btp')
        title(['Threshold <', num2str(fixetr)])

        savefig([savepath date,'_ClusteringCoeff_Norm_btp.fig'])
        f=gcf;
        exportgraphics(f,[savepath date '_ClusteringCoeff_Norm_btp.png'])

        %savefig([savepath date 'Metrics.fig'])
        %f=gcf;
        %exportgraphics(f,[savepath date 'Metrics.png'])

        %% Small worldness index%
        disp('Small Worldness index, running')

        SW_btp = CCmean_btp./LL_btp;
        %SWG1 = nanmean(SW_btp(idG1,:),1);
        %SWG2 = nanmean(SW_btp(idG2,:),1);

        SWNorm_btp = CCmeanNorm_btp./LLNorm_btp;
        %SWNormG1 = nanmean(SWNorm_btp(idG1,:),1);
        %SWNormG2 = nanmean(SWNorm_btp(idG2,:),1);

        figure
        subplot(1,2,1)
        hold on
        plot(itr,squeeze(nanmean(SW_btp(idG1,:),1)),'r','displayname','ctl')
        plot(itr,squeeze(nanmean(SW_btp(idG2,:),1)),'b','displayname','mal')
        xlabel('Threshold','fontsize',16)
        ylabel('Clustering/Path length','fontsize',16)
        title('Small World Index_btp')

        subplot(1,2,2)
        meand1 = mean(SW_btp(idG1,id));
        meand2 = mean(SW_btp(idG2,id));
        x = [1 2];
        y = [meand1 meand2];
        b = bar(x,y,'r');
        b.FaceColor = 'flat';
        b.CData(2,:) = [0 0 1];
        xlabel('Groups')
        ylabel('Small World Index_btp')
        title(['Threshold <', num2str(fixetr)])
        box on
        xlim([0 , 3])

        savefig([savepath date,'_SmallWorld_btp.fig'])
        f=gcf;
        exportgraphics(f,[savepath date 'SmallWorld_btp.png'])

        disp('Small Worldness index, done')
        save([savepath date '_SW_btp.mat'],'SW_btp','SWNorm_btp');

        figure
        subplot(1,2,1)
        hold on
        plot(itr,squeeze(nanmean(SWNorm_btp(idG1,:),1)),'r','displayname','ctl')
        plot(itr,squeeze(nanmean(SWNorm_btp(idG2,:),1)),'b','displayname','mal')
        xlabel('Threshold','fontsize',16)
        ylabel('Norm Clustering/Norm Path length','fontsize',16)
        title('Small World Index_btp')

        subplot(1,2,2)
        meand1 = mean(SWNorm_btp(idG1,id));
        meand2 = mean(SWNorm_btp(idG2,id));
        x = [1 2];
        y = [meand1 meand2];
        b = bar(x,y,'r');
        b.FaceColor = 'flat';
        b.CData(2,:) = [0 0 1];
        xlabel('Groups')
        ylabel('Small World Index Norm_btp')
        title(['Threshold <', num2str(fixetr)])
        box on
        xlim([0 , 3])

        savefig([savepath date,'_SmallWorldNorm_btp.fig'])
        f=gcf;
        exportgraphics(f,[savepath date 'SmallWorldNorm_btp.png'])
    end
end

if analysemode == 1
    %% ANOVA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p = 0.07;
    
    %% absolute threshold & weighted
    if absolutethresh ==1 && weightedmode ==1 %WTA
        %%%%%%%Effets principaux de Gr et SSE%%%
        x = 1;
        resGE_wta = [];
        for idtr = 1:numel(itr)
            tr = itr(idtr);
            [pGE_wta(:,x),res] = anovan(GE_wta(:,idtr),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
            close hidden;
            resGE_wta = [resGE_wta, res];
            x = x + 1;
        end

        n_sig = sum(pGE_wta(1,:) <= p);
        fprintf('%d Group effects are significant p<=%.2f without correction for GE_wta\n',n_sig, p)

        save([savepath date '_resultsANOVA_GE_wta.mat'],'pGE_wta', 'resGE_wta');
        clear n_sig res

        x = 1;
        resLL_wta = [];
        for idtr = 1:numel(itr)
            tr = itr(idtr);
            [pLL_wta(:,x),res] = anovan(LL_wta(:,idtr),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
            close hidden;
            resLL_wta = [resLL_wta, res];
            x = x + 1;
        end

        n_sig = sum(pLL_wta(1,:) <= p);
        fprintf('%d Group effects are significant p<=%.2f without correction for LL_wta\n',n_sig, p)

        save([savepath date '_resultsANOVA_LL_wta.mat'],'pLL_wta', 'resLL_wta');
        clear n_sig res

        x = 1;
        resSW_wta = [];
        for idtr = 1:numel(itr)-1
            tr = itr(idtr);
            [pSW_wta(:,x),res] = anovan(SW_wta(:,idtr),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
            close hidden;
            resSW_wta = [resSW_wta, res];
            x = x + 1;
        end

        n_sig = sum(pSW_wta(1,:) <= p);
        fprintf('%d Group effects are significant p<=%.2f without correction for SW_wta\n',n_sig, p)

        save([savepath date '_resultsANOVA_SW_wta.mat'],'pSW_wta', 'resSW_wta');
        clear n_sig res
        clear resGE_wta pGE_wta resLL_wta pLL_wta resSW_wta pSW_wta 

        %%%%%%%Ajout de l'interaction%%%
        x = 1;
        resGE_wta = [];
        for idtr = 1:numel(itr)
            tr = itr(idtr);
            [pGE_wta(:,x),res] = anovan(GE_wta(:,idtr),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
            close hidden;
            resGE_wta = [resGE_wta, res];
            x = x + 1;
        end

        n_sig = sum(pGE_wta(1,:) <= p);
        fprintf('%d Group effects are significant p<=%.2f without correction for GE_wta\n',n_sig, p)

        n_sig = sum(pGE_wta(3,:) <= p);
        fprintf('%d Group*SES effects are significant p<=%.2f without correction for GE_wta\n',n_sig, p)

        save([savepath date '_resultsANOVAinter_GE_wta.mat'],'pGE_wta', 'resGE_wta');
        clear n_sig res

        x = 1;
        resLL_wta = [];
        for idtr = 1:numel(itr)
            tr = itr(idtr);
            [pLL_wta(:,x),res] = anovan(LL_wta(:,idtr),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
            close hidden;
            resLL_wta = [resLL_wta, res];
            x = x + 1;
        end

        n_sig = sum(pLL_wta(1,:) <= p);
        fprintf('%d Group effects are significant p<=%.2f without correction for LL_wta\n',n_sig, p)

        n_sig = sum(pLL_wta(3,:) <= p);
        fprintf('%d Group*SES effects are significant p<=%.2f without correction for LL_wta\n',n_sig, p)

        save([savepath date '_resultsANOVAinter_LL_wta.mat'],'pLL_wta', 'resLL_wta');
        clear n_sig res

        x = 1;
        resSW_wta = [];
        for idtr = 1:numel(itr)-1
            tr = itr(idtr);
            [pSW_wta(:,x),res] = anovan(SW_wta(:,idtr),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
            close hidden;
            resSW_wta = [resSW_wta, res];
            x = x + 1;
        end

        n_sig = sum(pSW_wta(1,:) <= p);
        fprintf('%d Group effects are significant p<=%.2f without correction for SW_wta\n',n_sig, p)

        n_sig = sum(pSW_wta(3,:) <= p);
        fprintf('%d Group*SES effects are significant p<=%.2f without correction for SW_wta\n',n_sig, p)

        save([savepath date '_resultsANOVAinter_SW_wta.mat'],'pSW_wta', 'resSW_wta');
        clear n_sig res
        clear resGE_wta pGE_wta resLL_wta pLL_wta resSW_wta pSW_wta 

        %%%%%%%Sans covarier pour SES%%%
        x = 1;
        resGE_wta = [];
        for idtr = 1:numel(itr)
            tr = itr(idtr);
            [pGE_wta(:,x),res] = anovan(GE_wta(:,idtr),{gr},'varnames',{'Group'});
            close hidden;
            resGE_wta = [resGE_wta, res];
            x = x + 1;
        end

        n_sig = sum(pGE_wta(1,:) <= p);
        fprintf('%d Group effects are significant p<=%.2f without correction for GE_wta\n',n_sig, p)

        save([savepath date '_resultsANOVApascov_GE_wta.mat'],'pGE_wta', 'resGE_wta');
        clear n_sig res

        x = 1;
        resLL_wta = [];
        for idtr = 1:numel(itr)
            tr = itr(idtr);
            [pLL_wta(:,x),res] = anovan(LL_wta(:,idtr),{gr},'varnames',{'Group'});
            close hidden;
            resLL_wta = [resLL_wta, res];
            x = x + 1;
        end

        n_sig = sum(pLL_wta(1,:) <= p);
        fprintf('%d Group effects are significant p<=%.2f without correction for LL_wta\n',n_sig, p)

        save([savepath date '_resultsANOVApascov_LL_wta.mat'],'pLL_wta', 'resLL_wta');
        clear n_sig res
        
        x = 1;
        resSW_wta = [];
        for idtr = 1:numel(itr)-1
            tr = itr(idtr);
            [pSW_wta(:,x),res] = anovan(SW_wta(:,idtr),{gr},'varnames',{'Group'});
            close hidden;
            resSW_wta = [resSW_wta, res];
            x = x + 1;
        end
        
        n_sig = sum(pSW_wta(1,:) <= p);
        fprintf('%d Group effects are significant p<=%.2f without correction for SW_wta\n',n_sig, p)
        
        save([savepath date '_resultsANOVApascov_SW_wta.mat'],'pSW_wta', 'resSW_wta');
        clear n_sig res
        clear resGE_wta pGE_wta resLL_wta pLL_wta resSW_wta pSW_wta
        
        if normalizemode == 1
            %%%%%%%Effets principaux de Gr et SSE%%%
            x = 1;
            resGENorm_wta = [];
            for idtr = 1:numel(itr)
                tr = itr(idtr);
                [pGENorm_wta(:,x),res] = anovan(GENorm_wta(:,idtr),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
                close hidden;
                resGENorm_wta = [resGENorm_wta, res];
                x = x + 1;
            end
            
            n_sig = sum(pGENorm_wta(1,:) <= p);
            fprintf('%d Group effects are significant p<=%.2f without correction for GENorm_wta\n',n_sig, p)
            
            save([savepath date '_resultsANOVA_GENorm_wta.mat'],'pGENorm_wta', 'resGENorm_wta');
            clear n_sig res
            
            x = 1;
            resLLNorm_wta = [];
            for idtr = 1:numel(itr)
                tr = itr(idtr);
                [pLLNorm_wta(:,x),res] = anovan(LLNorm_wta(:,idtr),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
                close hidden;
                resLLNorm_wta = [resLLNorm_wta, res];
                x = x + 1;
            end
            
            n_sig = sum(pLLNorm_wta(1,:) <= p);
            fprintf('%d Group effects are significant p<=%.2f without correction for LLNorm_wta\n',n_sig, p)
            
            save([savepath date '_resultsANOVA_LLNorm_wta.mat'],'pLLNorm_wta', 'resLLNorm_wta');
            clear n_sig res
            
            x = 1;
            resSWNorm_wta = [];
            for idtr = 1:numel(itr)-1
                tr = itr(idtr);
                [pSWNorm_wta(:,x),res] = anovan(SWNorm_wta(:,idtr),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
                close hidden;
                resSWNorm_wta = [resSWNorm_wta, res];
                x = x + 1;
            end
            
            n_sig = sum(pSWNorm_wta(1,:) <= p);
            fprintf('%d Group effects are significant p<=%.2f without correction for SWNorm_wta\n',n_sig, p)
            
            save([savepath date '_resultsANOVA_SWNorm_wta.mat'],'pSWNorm_wta', 'resSWNorm_wta');
            clear n_sig res
            clear resGENorm_wta pGENorm_wta resLLNorm_wta pLLNorm_wta resSWNorm_wta pSWNorm_wta
            
            %%%%%%%Ajout de l'interaction%%%
            x = 1;
            resGENorm_wta = [];
            for idtr = 1:numel(itr)
                tr = itr(idtr);
                [pGENorm_wta(:,x),res] = anovan(GENorm_wta(:,idtr),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
                close hidden;
                resGENorm_wta = [resGENorm_wta, res];
                x = x + 1;
            end
            
            n_sig = sum(pGENorm_wta(1,:) <= p);
            fprintf('%d Group effects are significant p<=%.2f without correction for GENorm_wta\n',n_sig, p)
            
            n_sig = sum(pGENorm_wta(3,:) <= p);
            fprintf('%d Group*SES effects are significant p<=%.2f without correction for GENorm_wta\n',n_sig, p)
            
            save([savepath date '_resultsANOVAinter_GENorm_wta.mat'],'pGENorm_wta', 'resGENorm_wta');
            clear n_sig res
            
            x = 1;
            resLLNorm_wta = [];
            for idtr = 1:numel(itr)
                tr = itr(idtr);
                [pLLNorm_wta(:,x),res] = anovan(LLNorm_wta(:,idtr),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
                close hidden;
                resLLNorm_wta = [resLLNorm_wta, res];
                x = x + 1;
            end
            
            n_sig = sum(pLLNorm_wta(1,:) <= p);
            fprintf('%d Group effects are significant p<=%.2f without correction for LLNorm_wta\n',n_sig, p)
            
            n_sig = sum(pLLNorm_wta(3,:) <= p);
            fprintf('%d Group*SES effects are significant p<=%.2f without correction for LLNorm_wta\n',n_sig, p)
            
            save([savepath date '_resultsANOVAinter_LLNorm_wta.mat'],'pLLNorm_wta', 'resLLNorm_wta');
            clear n_sig res
            
            x = 1;
            resSWNorm_wta = [];
            for idtr = 1:numel(itr)-1
                tr = itr(idtr);
                [pSWNorm_wta(:,x),res] = anovan(SWNorm_wta(:,idtr),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
                close hidden;
                resSWNorm_wta = [resSWNorm_wta, res];
                x = x + 1;
            end
            
            n_sig = sum(pSWNorm_wta(1,:) <= p);
            fprintf('%d Group effects are significant p<=%.2f without correction for SWNorm_wta\n',n_sig, p)
            
            n_sig = sum(pSWNorm_wta(3,:) <= p);
            fprintf('%d Group*SES effects are significant p<=%.2f without correction for SWNorm_wta\n',n_sig, p)
            
            save([savepath date '_resultsANOVAinter_SWNorm_wta.mat'],'pSWNorm_wta', 'resSWNorm_wta');
            clear n_sig res
            clear resGENorm_wta pGENorm_wta resLLNorm_wta pLLNorm_wta resSWNorm_wta pSWNorm_wta
            
            %%%%%%%Sans covarier pour SES%%%
            x = 1;
            resGENorm_wta = [];
            for idtr = 1:numel(itr)
                tr = itr(idtr);
                [pGENorm_wta(:,x),res] = anovan(GENorm_wta(:,idtr),{gr},'varnames',{'Group'});
                close hidden;
                resGENorm_wta = [resGENorm_wta, res];
                x = x + 1;
            end
            
            n_sig = sum(pGENorm_wta(1,:) <= p);
            fprintf('%d Group effects are significant p<=%.2f without correction for GENorm_wta\n',n_sig, p)
            
            save([savepath date '_resultsANOVApascov_GENorm_wta.mat'],'pGENorm_wta', 'resGENorm_wta');
            clear n_sig res
            
            x = 1;
            resLLNorm_wta = [];
            for idtr = 1:numel(itr)
                tr = itr(idtr);
                [pLLNorm_wta(:,x),res] = anovan(LLNorm_wta(:,idtr),{gr},'varnames',{'Group'});
                close hidden;
                resLLNorm_wta = [resLLNorm_wta, res];
                x = x + 1;
            end
            
            n_sig = sum(pLLNorm_wta(1,:) <= p);
            fprintf('%d Group effects are significant p<=%.2f without correction for LLNorm_wta\n',n_sig, p)
            
            save([savepath date '_resultsANOVApascov_LLNorm_wta.mat'],'pLLNorm_wta', 'resLLNorm_wta');
            clear n_sig res
            
            x = 1;
            resSWNorm_wta = [];
            for idtr = 1:numel(itr)-1
                tr = itr(idtr);
                [pSWNorm_wta(:,x),res] = anovan(SWNorm_wta(:,idtr),{gr},'varnames',{'Group'});
                close hidden;
                resSWNorm_wta = [resSWNorm_wta, res];
                x = x + 1;
            end
            
            n_sig = sum(pSWNorm_wta(1,:) <= p);
            fprintf('%d Group effects are significant p<=%.2f without correction for SWNorm_wta\n',n_sig, p)
            
            save([savepath date '_resultsANOVApascov_SWNorm_wta.mat'],'pSWNorm_wta', 'resSWNorm_wta');
            clear n_sig res
            clear resGENorm_wta pGENorm_wta resLLNorm_wta pLLNorm_wta resSWNorm_wta pSWNorm_wta
        end
    end
    
    %% absolute threshold & binarized
   if absolutethresh ==1 && binarizedmode ==1 %BTA
      %%%%%%%Effets principaux de Gr et SSE%%%
        x = 1;
        resGE_bta = [];
        for idtr = 1:numel(itr)
            tr = itr(idtr);
            [pGE_bta(:,x),res] = anovan(GE_bta(:,idtr),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
            close hidden;
            resGE_bta = [resGE_bta, res];
            x = x + 1;
        end

        n_sig = sum(pGE_bta(1,:) <= p);
        fprintf('%d Group effects are significant p<=%.2f without correction for GE_bta\n',n_sig, p)

        save([savepath date '_resultsANOVA_GE_bta.mat'],'pGE_bta', 'resGE_bta');
        clear n_sig res

        x = 1;
        resLL_bta = [];
        for idtr = 1:numel(itr)
            tr = itr(idtr);
            [pLL_bta(:,x),res] = anovan(LL_bta(:,idtr),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
            close hidden;
            resLL_bta = [resLL_bta, res];
            x = x + 1;
        end

        n_sig = sum(pLL_bta(1,:) <= p);
        fprintf('%d Group effects are significant p<=%.2f without correction for LL_bta\n',n_sig, p)

        save([savepath date '_resultsANOVA_LL_bta.mat'],'pLL_bta', 'resLL_bta');
        clear n_sig res

        x = 1;
        resSW_bta = [];
        for idtr = 1:numel(itr)-1
            tr = itr(idtr);
            [pSW_bta(:,x),res] = anovan(SW_bta(:,idtr),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
            close hidden;
            resSW_bta = [resSW_bta, res];
            x = x + 1;
        end

        n_sig = sum(pSW_bta(1,:) <= p);
        fprintf('%d Group effects are significant p<=%.2f without correction for SW_bta\n',n_sig, p)

        save([savepath date '_resultsANOVA_SW_bta.mat'],'pSW_bta', 'resSW_bta');
        clear n_sig res
        clear resGE_bta pGE_bta resLL_bta pLL_bta resSW_bta pSW_bta 

        %%%%%%%Ajout de l'interaction%%%
        x = 1;
        resGE_bta = [];
        for idtr = 1:numel(itr)
            tr = itr(idtr);
            [pGE_bta(:,x),res] = anovan(GE_bta(:,idtr),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
            close hidden;
            resGE_bta = [resGE_bta, res];
            x = x + 1;
        end

        n_sig = sum(pGE_bta(1,:) <= p);
        fprintf('%d Group effects are significant p<=%.2f without correction for GE_bta\n',n_sig, p)

        n_sig = sum(pGE_bta(3,:) <= p);
        fprintf('%d Group*SES effects are significant p<=%.2f without correction for GE_bta\n',n_sig, p)

        save([savepath date '_resultsANOVAinter_GE_bta.mat'],'pGE_bta', 'resGE_bta');
        clear n_sig res

        x = 1;
        resLL_bta = [];
        for idtr = 1:numel(itr)
            tr = itr(idtr);
            [pLL_bta(:,x),res] = anovan(LL_bta(:,idtr),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
            close hidden;
            resLL_bta = [resLL_bta, res];
            x = x + 1;
        end

        n_sig = sum(pLL_bta(1,:) <= p);
        fprintf('%d Group effects are significant p<=%.2f without correction for LL_bta\n',n_sig, p)

        n_sig = sum(pLL_bta(3,:) <= p);
        fprintf('%d Group*SES effects are significant p<=%.2f without correction for LL_bta\n',n_sig, p)

        save([savepath date '_resultsANOVAinter_LL_bta.mat'],'pLL_bta', 'resLL_bta');
        clear n_sig res

        x = 1;
        resSW_bta = [];
        for idtr = 1:numel(itr)-1
            tr = itr(idtr);
            [pSW_bta(:,x),res] = anovan(SW_bta(:,idtr),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
            close hidden;
            resSW_bta = [resSW_bta, res];
            x = x + 1;
        end

        n_sig = sum(pSW_bta(1,:) <= p);
        fprintf('%d Group effects are significant p<=%.2f without correction for SW_bta\n',n_sig, p)

        n_sig = sum(pSW_bta(3,:) <= p);
        fprintf('%d Group*SES effects are significant p<=%.2f without correction for SW_bta\n',n_sig, p)

        save([savepath date '_resultsANOVAinter_SW_bta.mat'],'pSW_bta', 'resSW_bta');
        clear n_sig res
        clear resGE_bta pGE_bta resLL_bta pLL_bta resSW_bta pSW_bta 

        %%%%%%%Sans covarier pour SES%%%
        x = 1;
        resGE_bta = [];
        for idtr = 1:numel(itr)
            tr = itr(idtr);
            [pGE_bta(:,x),res] = anovan(GE_bta(:,idtr),{gr},'varnames',{'Group'});
            close hidden;
            resGE_bta = [resGE_bta, res];
            x = x + 1;
        end

        n_sig = sum(pGE_bta(1,:) <= p);
        fprintf('%d Group effects are significant p<=%.2f without correction for GE_bta\n',n_sig, p)

        save([savepath date '_resultsANOVApascov_GE_bta.mat'],'pGE_bta', 'resGE_bta');
        clear n_sig res

        x = 1;
        resLL_bta = [];
        for idtr = 1:numel(itr)
            tr = itr(idtr);
            [pLL_bta(:,x),res] = anovan(LL_bta(:,idtr),{gr},'varnames',{'Group'});
            close hidden;
            resLL_bta = [resLL_bta, res];
            x = x + 1;
        end

        n_sig = sum(pLL_bta(1,:) <= p);
        fprintf('%d Group effects are significant p<=%.2f without correction for LL_bta\n',n_sig, p)

        save([savepath date '_resultsANOVApascov_LL_bta.mat'],'pLL_bta', 'resLL_bta');
        clear n_sig res

        x = 1;
        resSW_bta = [];
        for idtr = 1:numel(itr)-1
            tr = itr(idtr);
            [pSW_bta(:,x),res] = anovan(SW_bta(:,idtr),{gr},'varnames',{'Group'});
            close hidden;
            resSW_bta = [resSW_bta, res];
            x = x + 1;
        end

        n_sig = sum(pSW_bta(1,:) <= p);
        fprintf('%d Group effects are significant p<=%.2f without correction for SW_bta\n',n_sig, p)

        save([savepath date '_resultsANOVApascov_SW_bta.mat'],'pSW_bta', 'resSW_bta');
        clear n_sig res
        clear resGE_bta pGE_bta resLL_bta pLL_bta resSW_bta pSW_bta
        
        if normalizemode == 1
            %%%%%%%Effets principaux de Gr et SSE%%%
            x = 1;
            resGENorm_bta = [];
            for idtr = 1:numel(itr)
                tr = itr(idtr);
                [pGENorm_bta(:,x),res] = anovan(GENorm_bta(:,idtr),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
                close hidden;
                resGENorm_bta = [resGENorm_bta, res];
                x = x + 1;
            end
            
            n_sig = sum(pGENorm_bta(1,:) <= p);
            fprintf('%d Group effects are significant p<=%.2f without correction for GENorm_bta\n',n_sig, p)
            
            save([savepath date '_resultsANOVA_GENorm_bta.mat'],'pGENorm_bta', 'resGENorm_bta');
            clear n_sig res
            
            x = 1;
            resLLNorm_bta = [];
            for idtr = 1:numel(itr)
                tr = itr(idtr);
                [pLLNorm_bta(:,x),res] = anovan(LLNorm_bta(:,idtr),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
                close hidden;
                resLLNorm_bta = [resLLNorm_bta, res];
                x = x + 1;
            end
            
            n_sig = sum(pLLNorm_bta(1,:) <= p);
            fprintf('%d Group effects are significant p<=%.2f without correction for LLNorm_bta\n',n_sig, p)
            
            save([savepath date '_resultsANOVA_LLNorm_bta.mat'],'pLLNorm_bta', 'resLLNorm_bta');
            clear n_sig res
            
            x = 1;
            resSWNorm_bta = [];
            for idtr = 1:numel(itr)-1
                tr = itr(idtr);
                [pSWNorm_bta(:,x),res] = anovan(SWNorm_bta(:,idtr),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
                close hidden;
                resSWNorm_bta = [resSWNorm_bta, res];
                x = x + 1;
            end
            
            n_sig = sum(pSWNorm_bta(1,:) <= p);
            fprintf('%d Group effects are significant p<=%.2f without correction for SWNorm_bta\n',n_sig, p)
            
            save([savepath date '_resultsANOVA_SWNorm_bta.mat'],'pSWNorm_bta', 'resSWNorm_bta');
            clear n_sig res
            clear resGENorm_bta pGENorm_bta resLLNorm_bta pLLNorm_bta resSWNorm_bta pSWNorm_bta
            
            %%%%%%%Ajout de l'interaction%%%
            x = 1;
            resGENorm_bta = [];
            for idtr = 1:numel(itr)
                tr = itr(idtr);
                [pGENorm_bta(:,x),res] = anovan(GENorm_bta(:,idtr),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
                close hidden;
                resGENorm_bta = [resGENorm_bta, res];
                x = x + 1;
            end
            
            n_sig = sum(pGENorm_bta(1,:) <= p);
            fprintf('%d Group effects are significant p<=%.2f without correction for GENorm_bta\n',n_sig, p)
            
            n_sig = sum(pGENorm_bta(3,:) <= p);
            fprintf('%d Group*SES effects are significant p<=%.2f without correction for GENorm_bta\n',n_sig, p)
            
            save([savepath date '_resultsANOVAinter_GENorm_bta.mat'],'pGENorm_bta', 'resGENorm_bta');
            clear n_sig res
            
            x = 1;
            resLLNorm_bta = [];
            for idtr = 1:numel(itr)
                tr = itr(idtr);
                [pLLNorm_bta(:,x),res] = anovan(LLNorm_bta(:,idtr),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
                close hidden;
                resLLNorm_bta = [resLLNorm_bta, res];
                x = x + 1;
            end
            
            n_sig = sum(pLLNorm_bta(1,:) <= p);
            fprintf('%d Group effects are significant p<=%.2f without correction for LLNorm_bta\n',n_sig, p)
            
            n_sig = sum(pLLNorm_bta(3,:) <= p);
            fprintf('%d Group*SES effects are significant p<=%.2f without correction for LLNorm_bta\n',n_sig, p)
            
            save([savepath date '_resultsANOVAinter_LLNorm_bta.mat'],'pLLNorm_bta', 'resLLNorm_bta');
            clear n_sig res
            
            x = 1;
            resSWNorm_bta = [];
            for idtr = 1:numel(itr)-1
                tr = itr(idtr);
                [pSWNorm_bta(:,x),res] = anovan(SWNorm_bta(:,idtr),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
                close hidden;
                resSWNorm_bta = [resSWNorm_bta, res];
                x = x + 1;
            end
            
            n_sig = sum(pSWNorm_bta(1,:) <= p);
            fprintf('%d Group effects are significant p<=%.2f without correction for SWNorm_bta\n',n_sig, p)
            
            n_sig = sum(pSWNorm_bta(3,:) <= p);
            fprintf('%d Group*SES effects are significant p<=%.2f without correction for SWNorm_bta\n',n_sig, p)
            
            save([savepath date '_resultsANOVAinter_SWNorm_bta.mat'],'pSWNorm_bta', 'resSWNorm_bta');
            clear n_sig res
            clear resGENorm_bta pGENorm_bta resLLNorm_bta pLLNorm_bta resSWNorm_bta pSWNorm_bta
            
            %%%%%%%Sans covarier pour SES%%%
            x = 1;
            resGENorm_bta = [];
            for idtr = 1:numel(itr)
                tr = itr(idtr);
                [pGENorm_bta(:,x),res] = anovan(GENorm_bta(:,idtr),{gr},'varnames',{'Group'});
                close hidden;
                resGENorm_bta = [resGENorm_bta, res];
                x = x + 1;
            end
            
            n_sig = sum(pGENorm_bta(1,:) <= p);
            fprintf('%d Group effects are significant p<=%.2f without correction for GENorm_bta\n',n_sig, p)
            
            save([savepath date '_resultsANOVApascov_GENorm_bta.mat'],'pGENorm_bta', 'resGENorm_bta');
            clear n_sig res
            
            x = 1;
            resLLNorm_bta = [];
            for idtr = 1:numel(itr)
                tr = itr(idtr);
                [pLLNorm_bta(:,x),res] = anovan(LLNorm_bta(:,idtr),{gr},'varnames',{'Group'});
                close hidden;
                resLLNorm_bta = [resLLNorm_bta, res];
                x = x + 1;
            end
            
            n_sig = sum(pLLNorm_bta(1,:) <= p);
            fprintf('%d Group effects are significant p<=%.2f without correction for LLNorm_bta\n',n_sig, p)
            
            save([savepath date '_resultsANOVApascov_LLNorm_bta.mat'],'pLLNorm_bta', 'resLLNorm_bta');
            clear n_sig res
            
            x = 1;
            resSWNorm_bta = [];
            for idtr = 1:numel(itr)-1
                tr = itr(idtr);
                [pSWNorm_bta(:,x),res] = anovan(SWNorm_bta(:,idtr),{gr},'varnames',{'Group'});
                close hidden;
                resSWNorm_bta = [resSWNorm_bta, res];
                x = x + 1;
            end
            
            n_sig = sum(pSWNorm_bta(1,:) <= p);
            fprintf('%d Group effects are significant p<=%.2f without correction for SWNorm_bta\n',n_sig, p)
            
            save([savepath date '_resultsANOVApascov_SWNorm_bta.mat'],'pSWNorm_bta', 'resSWNorm_bta');
            clear n_sig res
            clear resGENorm_bta pGENorm_bta resLLNorm_bta pLLNorm_bta resSWNorm_bta pSWNorm_bta
        end
   end
   
   %% proportional threshold & weighted
   if proportionalthresh ==1 && weightedmode ==1 %WTP
      %%%%%%%Effets principaux de Gr et SSE%%%
        x = 1;
        resGE_wtp = [];
        for idtr = 1:numel(itr)
            tr = itr(idtr);
            [pGE_wtp(:,x),res] = anovan(GE_wtp(:,idtr),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
            close hidden;
            resGE_wtp = [resGE_wtp, res];
            x = x + 1;
        end

        n_sig = sum(pGE_wtp(1,:) <= p);
        fprintf('%d Group effects are significant p<=%.2f without correction for GE_wtp\n',n_sig, p)

        save([savepath date '_resultsANOVA_GE_wtp.mat'],'pGE_wtp', 'resGE_wtp');
        clear n_sig res

        x = 1;
        resLL_wtp = [];
        for idtr = 1:numel(itr)
            tr = itr(idtr);
            [pLL_wtp(:,x),res] = anovan(LL_wtp(:,idtr),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
            close hidden;
            resLL_wtp = [resLL_wtp, res];
            x = x + 1;
        end

        n_sig = sum(pLL_wtp(1,:) <= p);
        fprintf('%d Group effects are significant p<=%.2f without correction for LL_wtp\n',n_sig, p)

        save([savepath date '_resultsANOVA_LL_wtp.mat'],'pLL_wtp', 'resLL_wtp');
        clear n_sig res

        x = 1;
        resSW_wtp = [];
        for idtr = 1:numel(itr)-1
            tr = itr(idtr);
            [pSW_wtp(:,x),res] = anovan(SW_wtp(:,idtr),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
            close hidden;
            resSW_wtp = [resSW_wtp, res];
            x = x + 1;
        end

        n_sig = sum(pSW_wtp(1,:) <= p);
        fprintf('%d Group effects are significant p<=%.2f without correction for SW_wtp\n',n_sig, p)

        save([savepath date '_resultsANOVA_SW_wtp.mat'],'pSW_wtp', 'resSW_wtp');
        clear n_sig res
        clear resGE_wtp pGE_wtp resLL_wtp pLL_wtp resSW_wtp pSW_wtp 

        %%%%%%%Ajout de l'interaction%%%
        x = 1;
        resGE_wtp = [];
        for idtr = 1:numel(itr)
            tr = itr(idtr);
            [pGE_wtp(:,x),res] = anovan(GE_wtp(:,idtr),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
            close hidden;
            resGE_wtp = [resGE_wtp, res];
            x = x + 1;
        end

        n_sig = sum(pGE_wtp(1,:) <= p);
        fprintf('%d Group effects are significant p<=%.2f without correction for GE_wtp\n',n_sig, p)

        n_sig = sum(pGE_wtp(3,:) <= p);
        fprintf('%d Group*SES effects are significant p<=%.2f without correction for GE_wtp\n',n_sig, p)

        save([savepath date '_resultsANOVAinter_GE_wtp.mat'],'pGE_wtp', 'resGE_wtp');
        clear n_sig res

        x = 1;
        resLL_wtp = [];
        for idtr = 1:numel(itr)
            tr = itr(idtr);
            [pLL_wtp(:,x),res] = anovan(LL_wtp(:,idtr),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
            close hidden;
            resLL_wtp = [resLL_wtp, res];
            x = x + 1;
        end

        n_sig = sum(pLL_wtp(1,:) <= p);
        fprintf('%d Group effects are significant p<=%.2f without correction for LL_wtp\n',n_sig, p)

        n_sig = sum(pLL_wtp(3,:) <= p);
        fprintf('%d Group*SES effects are significant p<=%.2f without correction for LL_wtp\n',n_sig, p)

        save([savepath date '_resultsANOVAinter_LL_wtp.mat'],'pLL_wtp', 'resLL_wtp');
        clear n_sig res

        x = 1;
        resSW_wtp = [];
        for idtr = 1:numel(itr)-1
            tr = itr(idtr);
            [pSW_wtp(:,x),res] = anovan(SW_wtp(:,idtr),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
            close hidden;
            resSW_wtp = [resSW_wtp, res];
            x = x + 1;
        end

        n_sig = sum(pSW_wtp(1,:) <= p);
        fprintf('%d Group effects are significant p<=%.2f without correction for SW_wtp\n',n_sig, p)

        n_sig = sum(pSW_wtp(3,:) <= p);
        fprintf('%d Group*SES effects are significant p<=%.2f without correction for SW_wtp\n',n_sig, p)

        save([savepath date '_resultsANOVAinter_SW_wtp.mat'],'pSW_wtp', 'resSW_wtp');
        clear n_sig res
        clear resGE_wtp pGE_wtp resLL_wtp pLL_wtp resSW_wtp pSW_wtp 

        %%%%%%%Sans covarier pour SES%%%
        x = 1;
        resGE_wtp = [];
        for idtr = 1:numel(itr)
            tr = itr(idtr);
            [pGE_wtp(:,x),res] = anovan(GE_wtp(:,idtr),{gr},'varnames',{'Group'});
            close hidden;
            resGE_wtp = [resGE_wtp, res];
            x = x + 1;
        end

        n_sig = sum(pGE_wtp(1,:) <= p);
        fprintf('%d Group effects are significant p<=%.2f without correction for GE_wtp\n',n_sig, p)

        save([savepath date '_resultsANOVApascov_GE_wtp.mat'],'pGE_wtp', 'resGE_wtp');
        clear n_sig res

        x = 1;
        resLL_wtp = [];
        for idtr = 1:numel(itr)
            tr = itr(idtr);
            [pLL_wtp(:,x),res] = anovan(LL_wtp(:,idtr),{gr},'varnames',{'Group'});
            close hidden;
            resLL_wtp = [resLL_wtp, res];
            x = x + 1;
        end

        n_sig = sum(pLL_wtp(1,:) <= p);
        fprintf('%d Group effects are significant p<=%.2f without correction for LL_wtp\n',n_sig, p)

        save([savepath date '_resultsANOVApascov_LL_wtp.mat'],'pLL_wtp', 'resLL_wtp');
        clear n_sig res
        
        x = 1;
        resSW_wtp = [];
        for idtr = 1:numel(itr)-1
            tr = itr(idtr);
            [pSW_wtp(:,x),res] = anovan(SW_wtp(:,idtr),{gr},'varnames',{'Group'});
            close hidden;
            resSW_wtp = [resSW_wtp, res];
            x = x + 1;
        end
        
        n_sig = sum(pSW_wtp(1,:) <= p);
        fprintf('%d Group effects are significant p<=%.2f without correction for SW_wtp\n',n_sig, p)
        
        save([savepath date '_resultsANOVApascov_SW_wtp.mat'],'pSW_wtp', 'resSW_wtp');
        clear n_sig res
        clear resGE_wtp pGE_wtp resLL_wtp pLL_wtp resSW_wtp pSW_wtp
        
        if normalizemode == 1
            %%%%%%%Effets principaux de Gr et SSE%%%
            x = 1;
            resGENorm_wtp = [];
            for idtr = 1:numel(itr)
                tr = itr(idtr);
                [pGENorm_wtp(:,x),res] = anovan(GENorm_wtp(:,idtr),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
                close hidden;
                resGENorm_wtp = [resGENorm_wtp, res];
                x = x + 1;
            end
            
            n_sig = sum(pGENorm_wtp(1,:) <= p);
            fprintf('%d Group effects are significant p<=%.2f without correction for GENorm_wtp\n',n_sig, p)
            
            save([savepath date '_resultsANOVA_GENorm_wtp.mat'],'pGENorm_wtp', 'resGENorm_wtp');
            clear n_sig res
            
            x = 1;
            resLLNorm_wtp = [];
            for idtr = 1:numel(itr)
                tr = itr(idtr);
                [pLLNorm_wtp(:,x),res] = anovan(LLNorm_wtp(:,idtr),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
                close hidden;
                resLLNorm_wtp = [resLLNorm_wtp, res];
                x = x + 1;
            end
            
            n_sig = sum(pLLNorm_wtp(1,:) <= p);
            fprintf('%d Group effects are significant p<=%.2f without correction for LLNorm_wtp\n',n_sig, p)
            
            save([savepath date '_resultsANOVA_LLNorm_wtp.mat'],'pLLNorm_wtp', 'resLLNorm_wtp');
            clear n_sig res
            
            x = 1;
            resSWNorm_wtp = [];
            for idtr = 1:numel(itr)-1
                tr = itr(idtr);
                [pSWNorm_wtp(:,x),res] = anovan(SWNorm_wtp(:,idtr),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
                close hidden;
                resSWNorm_wtp = [resSWNorm_wtp, res];
                x = x + 1;
            end
            
            n_sig = sum(pSWNorm_wtp(1,:) <= p);
            fprintf('%d Group effects are significant p<=%.2f without correction for SWNorm_wtp\n',n_sig, p)
            
            save([savepath date '_resultsANOVA_SWNorm_wtp.mat'],'pSWNorm_wtp', 'resSWNorm_wtp');
            clear n_sig res
            clear resGENorm_wtp pGENorm_wtp resLLNorm_wtp pLLNorm_wtp resSWNorm_wtp pSWNorm_wtp
            
            %%%%%%%Ajout de l'interaction%%%
            x = 1;
            resGENorm_wtp = [];
            for idtr = 1:numel(itr)
                tr = itr(idtr);
                [pGENorm_wtp(:,x),res] = anovan(GENorm_wtp(:,idtr),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
                close hidden;
                resGENorm_wtp = [resGENorm_wtp, res];
                x = x + 1;
            end
            
            n_sig = sum(pGENorm_wtp(1,:) <= p);
            fprintf('%d Group effects are significant p<=%.2f without correction for GENorm_wtp\n',n_sig, p)
            
            n_sig = sum(pGENorm_wtp(3,:) <= p);
            fprintf('%d Group*SES effects are significant p<=%.2f without correction for GENorm_wtp\n',n_sig, p)
            
            save([savepath date '_resultsANOVAinter_GENorm_wtp.mat'],'pGENorm_wtp', 'resGENorm_wtp');
            clear n_sig res
            
            x = 1;
            resLLNorm_wtp = [];
            for idtr = 1:numel(itr)
                tr = itr(idtr);
                [pLLNorm_wtp(:,x),res] = anovan(LLNorm_wtp(:,idtr),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
                close hidden;
                resLLNorm_wtp = [resLLNorm_wtp, res];
                x = x + 1;
            end
            
            n_sig = sum(pLLNorm_wtp(1,:) <= p);
            fprintf('%d Group effects are significant p<=%.2f without correction for LLNorm_wtp\n',n_sig, p)
            
            n_sig = sum(pLLNorm_wtp(3,:) <= p);
            fprintf('%d Group*SES effects are significant p<=%.2f without correction for LLNorm_wtp\n',n_sig, p)
            
            save([savepath date '_resultsANOVAinter_LLNorm_wtp.mat'],'pLLNorm_wtp', 'resLLNorm_wtp');
            clear n_sig res
            
            x = 1;
            resSWNorm_wtp = [];
            for idtr = 1:numel(itr)-1
                tr = itr(idtr);
                [pSWNorm_wtp(:,x),res] = anovan(SWNorm_wtp(:,idtr),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
                close hidden;
                resSWNorm_wtp = [resSWNorm_wtp, res];
                x = x + 1;
            end
            
            n_sig = sum(pSWNorm_wtp(1,:) <= p);
            fprintf('%d Group effects are significant p<=%.2f without correction for SWNorm_wtp\n',n_sig, p)
            
            n_sig = sum(pSWNorm_wtp(3,:) <= p);
            fprintf('%d Group*SES effects are significant p<=%.2f without correction for SWNorm_wtp\n',n_sig, p)
            
            save([savepath date '_resultsANOVAinter_SWNorm_wtp.mat'],'pSWNorm_wtp', 'resSWNorm_wtp');
            clear n_sig res
            clear resGENorm_wtp pGENorm_wtp resLLNorm_wtp pLLNorm_wtp resSWNorm_wtp pSWNorm_wtp
            
            %%%%%%%Sans covarier pour SES%%%
            x = 1;
            resGENorm_wtp = [];
            for idtr = 1:numel(itr)
                tr = itr(idtr);
                [pGENorm_wtp(:,x),res] = anovan(GENorm_wtp(:,idtr),{gr},'varnames',{'Group'});
                close hidden;
                resGENorm_wtp = [resGENorm_wtp, res];
                x = x + 1;
            end
            
            n_sig = sum(pGENorm_wtp(1,:) <= p);
            fprintf('%d Group effects are significant p<=%.2f without correction for GENorm_wtp\n',n_sig, p)
            
            save([savepath date '_resultsANOVApascov_GENorm_wtp.mat'],'pGENorm_wtp', 'resGENorm_wtp');
            clear n_sig res
            
            x = 1;
            resLLNorm_wtp = [];
            for idtr = 1:numel(itr)
                tr = itr(idtr);
                [pLLNorm_wtp(:,x),res] = anovan(LLNorm_wtp(:,idtr),{gr},'varnames',{'Group'});
                close hidden;
                resLLNorm_wtp = [resLLNorm_wtp, res];
                x = x + 1;
            end
            
            n_sig = sum(pLLNorm_wtp(1,:) <= p);
            fprintf('%d Group effects are significant p<=%.2f without correction for LLNorm_wtp\n',n_sig, p)
            
            save([savepath date '_resultsANOVApascov_LLNorm_wtp.mat'],'pLLNorm_wtp', 'resLLNorm_wtp');
            clear n_sig res
            
            x = 1;
            resSWNorm_wtp = [];
            for idtr = 1:numel(itr)-1
                tr = itr(idtr);
                [pSWNorm_wtp(:,x),res] = anovan(SWNorm_wtp(:,idtr),{gr},'varnames',{'Group'});
                close hidden;
                resSWNorm_wtp = [resSWNorm_wtp, res];
                x = x + 1;
            end
            
            n_sig = sum(pSWNorm_wtp(1,:) <= p);
            fprintf('%d Group effects are significant p<=%.2f without correction for SWNorm_wtp\n',n_sig, p)
            
            save([savepath date '_resultsANOVApascov_SWNorm_wtp.mat'],'pSWNorm_wtp', 'resSWNorm_wtp');
            clear n_sig res
            clear resGENorm_wtp pGENorm_wtp resLLNorm_wtp pLLNorm_wtp resSWNorm_wtp pSWNorm_wtp
        end
   end
   
   %% proportional threshold & binarized
   if proportionalthresh ==1 && binarizedmode ==1 %BTP
      %%%%%%%Effets principaux de Gr et SSE%%%
        x = 1;
        resGE_btp = [];
        for idtr = 1:numel(itr)
            tr = itr(idtr);
            [pGE_btp(:,x),res] = anovan(GE_btp(:,idtr),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
            close hidden;
            resGE_btp = [resGE_btp, res];
            x = x + 1;
        end

        n_sig = sum(pGE_btp(1,:) <= p);
        fprintf('%d Group effects are significant p<=%.2f without correction for GE_btp\n',n_sig, p)

        save([savepath date '_resultsANOVA_GE_btp.mat'],'pGE_btp', 'resGE_btp');
        clear n_sig res

        x = 1;
        resLL_btp = [];
        for idtr = 1:numel(itr)
            tr = itr(idtr);
            [pLL_btp(:,x),res] = anovan(LL_btp(:,idtr),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
            close hidden;
            resLL_btp = [resLL_btp, res];
            x = x + 1;
        end

        n_sig = sum(pLL_btp(1,:) <= p);
        fprintf('%d Group effects are significant p<=%.2f without correction for LL_btp\n',n_sig, p)

        save([savepath date '_resultsANOVA_LL_btp.mat'],'pLL_btp', 'resLL_btp');
        clear n_sig res

        x = 1;
        resSW_btp = [];
        for idtr = 1:numel(itr)-1
            tr = itr(idtr);
            [pSW_btp(:,x),res] = anovan(SW_btp(:,idtr),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
            close hidden;
            resSW_btp = [resSW_btp, res];
            x = x + 1;
        end

        n_sig = sum(pSW_btp(1,:) <= p);
        fprintf('%d Group effects are significant p<=%.2f without correction for SW_btp\n',n_sig, p)

        save([savepath date '_resultsANOVA_SW_btp.mat'],'pSW_btp', 'resSW_btp');
        clear n_sig res
        clear resGE_btp pGE_btp resLL_btp pLL_btp resSW_btp pSW_btp 

        %%%%%%%Ajout de l'interaction%%%
        x = 1;
        resGE_btp = [];
        for idtr = 1:numel(itr)
            tr = itr(idtr);
            [pGE_btp(:,x),res] = anovan(GE_btp(:,idtr),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
            close hidden;
            resGE_btp = [resGE_btp, res];
            x = x + 1;
        end

        n_sig = sum(pGE_btp(1,:) <= p);
        fprintf('%d Group effects are significant p<=%.2f without correction for GE_btp\n',n_sig, p)

        n_sig = sum(pGE_btp(3,:) <= p);
        fprintf('%d Group*SES effects are significant p<=%.2f without correction for GE_btp\n',n_sig, p)

        save([savepath date '_resultsANOVAinter_GE_btp.mat'],'pGE_btp', 'resGE_btp');
        clear n_sig res

        x = 1;
        resLL_btp = [];
        for idtr = 1:numel(itr)
            tr = itr(idtr);
            [pLL_btp(:,x),res] = anovan(LL_btp(:,idtr),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
            close hidden;
            resLL_btp = [resLL_btp, res];
            x = x + 1;
        end

        n_sig = sum(pLL_btp(1,:) <= p);
        fprintf('%d Group effects are significant p<=%.2f without correction for LL_btp\n',n_sig, p)

        n_sig = sum(pLL_btp(3,:) <= p);
        fprintf('%d Group*SES effects are significant p<=%.2f without correction for LL_btp\n',n_sig, p)

        save([savepath date '_resultsANOVAinter_LL_btp.mat'],'pLL_btp', 'resLL_btp');
        clear n_sig res

        x = 1;
        resSW_btp = [];
        for idtr = 1:numel(itr)-1
            tr = itr(idtr);
            [pSW_btp(:,x),res] = anovan(SW_btp(:,idtr),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
            close hidden;
            resSW_btp = [resSW_btp, res];
            x = x + 1;
        end

        n_sig = sum(pSW_btp(1,:) <= p);
        fprintf('%d Group effects are significant p<=%.2f without correction for SW_btp\n',n_sig, p)

        n_sig = sum(pSW_btp(3,:) <= p);
        fprintf('%d Group*SES effects are significant p<=%.2f without correction for SW_btp\n',n_sig, p)

        save([savepath date '_resultsANOVAinter_SW_btp.mat'],'pSW_btp', 'resSW_btp');
        clear n_sig res
        clear resGE_btp pGE_btp resLL_btp pLL_btp resSW_btp pSW_btp 

        %%%%%%%Sans covarier pour SES%%%
        x = 1;
        resGE_btp = [];
        for idtr = 1:numel(itr)
            tr = itr(idtr);
            [pGE_btp(:,x),res] = anovan(GE_btp(:,idtr),{gr},'varnames',{'Group'});
            close hidden;
            resGE_btp = [resGE_btp, res];
            x = x + 1;
        end

        n_sig = sum(pGE_btp(1,:) <= p);
        fprintf('%d Group effects are significant p<=%.2f without correction for GE_btp\n',n_sig, p)

        save([savepath date '_resultsANOVApascov_GE_btp.mat'],'pGE_btp', 'resGE_btp');
        clear n_sig res

        x = 1;
        resLL_btp = [];
        for idtr = 1:numel(itr)
            tr = itr(idtr);
            [pLL_btp(:,x),res] = anovan(LL_btp(:,idtr),{gr},'varnames',{'Group'});
            close hidden;
            resLL_btp = [resLL_btp, res];
            x = x + 1;
        end

        n_sig = sum(pLL_btp(1,:) <= p);
        fprintf('%d Group effects are significant p<=%.2f without correction for LL_btp\n',n_sig, p)

        save([savepath date '_resultsANOVApascov_LL_btp.mat'],'pLL_btp', 'resLL_btp');
        clear n_sig res
        
        x = 1;
        resSW_btp = [];
        for idtr = 1:numel(itr)-1
            tr = itr(idtr);
            [pSW_btp(:,x),res] = anovan(SW_btp(:,idtr),{gr},'varnames',{'Group'});
            close hidden;
            resSW_btp = [resSW_btp, res];
            x = x + 1;
        end
        
        n_sig = sum(pSW_btp(1,:) <= p);
        fprintf('%d Group effects are significant p<=%.2f without correction for SW_btp\n',n_sig, p)
        
        save([savepath date '_resultsANOVApascov_SW_btp.mat'],'pSW_btp', 'resSW_btp');
        clear n_sig res
        clear resGE_btp pGE_btp resLL_btp pLL_btp resSW_btp pSW_btp
        
        if normalizemode == 1
            %%%%%%%Effets principaux de Gr et SSE%%%
            x = 1;
            resGENorm_btp = [];
            for idtr = 1:numel(itr)
                tr = itr(idtr);
                [pGENorm_btp(:,x),res] = anovan(GENorm_btp(:,idtr),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
                close hidden;
                resGENorm_btp = [resGENorm_btp, res];
                x = x + 1;
            end
            
            n_sig = sum(pGENorm_btp(1,:) <= p);
            fprintf('%d Group effects are significant p<=%.2f without correction for GENorm_btp\n',n_sig, p)
            
            save([savepath date '_resultsANOVA_GENorm_btp.mat'],'pGENorm_btp', 'resGENorm_btp');
            clear n_sig res
            
            x = 1;
            resLLNorm_btp = [];
            for idtr = 1:numel(itr)
                tr = itr(idtr);
                [pLLNorm_btp(:,x),res] = anovan(LLNorm_btp(:,idtr),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
                close hidden;
                resLLNorm_btp = [resLLNorm_btp, res];
                x = x + 1;
            end
            
            n_sig = sum(pLLNorm_btp(1,:) <= p);
            fprintf('%d Group effects are significant p<=%.2f without correction for LLNorm_btp\n',n_sig, p)
            
            save([savepath date '_resultsANOVA_LLNorm_btp.mat'],'pLLNorm_btp', 'resLLNorm_btp');
            clear n_sig res
            
            x = 1;
            resSWNorm_btp = [];
            for idtr = 1:numel(itr)-1
                tr = itr(idtr);
                [pSWNorm_btp(:,x),res] = anovan(SWNorm_btp(:,idtr),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
                close hidden;
                resSWNorm_btp = [resSWNorm_btp, res];
                x = x + 1;
            end
            
            n_sig = sum(pSWNorm_btp(1,:) <= p);
            fprintf('%d Group effects are significant p<=%.2f without correction for SWNorm_btp\n',n_sig, p)
            
            save([savepath date '_resultsANOVA_SWNorm_btp.mat'],'pSWNorm_btp', 'resSWNorm_btp');
            clear n_sig res
            clear resGENorm_btp pGENorm_btp resLLNorm_btp pLLNorm_btp resSWNorm_btp pSWNorm_btp
            
            %%%%%%%Ajout de l'interaction%%%
            x = 1;
            resGENorm_btp = [];
            for idtr = 1:numel(itr)
                tr = itr(idtr);
                [pGENorm_btp(:,x),res] = anovan(GENorm_btp(:,idtr),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
                close hidden;
                resGENorm_btp = [resGENorm_btp, res];
                x = x + 1;
            end
            
            n_sig = sum(pGENorm_btp(1,:) <= p);
            fprintf('%d Group effects are significant p<=%.2f without correction for GENorm_btp\n',n_sig, p)
            
            n_sig = sum(pGENorm_btp(3,:) <= p);
            fprintf('%d Group*SES effects are significant p<=%.2f without correction for GENorm_btp\n',n_sig, p)
            
            save([savepath date '_resultsANOVAinter_GENorm_btp.mat'],'pGENorm_btp', 'resGENorm_btp');
            clear n_sig res
            
            x = 1;
            resLLNorm_btp = [];
            for idtr = 1:numel(itr)
                tr = itr(idtr);
                [pLLNorm_btp(:,x),res] = anovan(LLNorm_btp(:,idtr),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
                close hidden;
                resLLNorm_btp = [resLLNorm_btp, res];
                x = x + 1;
            end
            
            n_sig = sum(pLLNorm_btp(1,:) <= p);
            fprintf('%d Group effects are significant p<=%.2f without correction for LLNorm_btp\n',n_sig, p)
            
            n_sig = sum(pLLNorm_btp(3,:) <= p);
            fprintf('%d Group*SES effects are significant p<=%.2f without correction for LLNorm_btp\n',n_sig, p)
            
            save([savepath date '_resultsANOVAinter_LLNorm_btp.mat'],'pLLNorm_btp', 'resLLNorm_btp');
            clear n_sig res
            
            x = 1;
            resSWNorm_btp = [];
            for idtr = 1:numel(itr)-1
                tr = itr(idtr);
                [pSWNorm_btp(:,x),res] = anovan(SWNorm_btp(:,idtr),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
                close hidden;
                resSWNorm_btp = [resSWNorm_btp, res];
                x = x + 1;
            end
            
            n_sig = sum(pSWNorm_btp(1,:) <= p);
            fprintf('%d Group effects are significant p<=%.2f without correction for SWNorm_btp\n',n_sig, p)
            
            n_sig = sum(pSWNorm_btp(3,:) <= p);
            fprintf('%d Group*SES effects are significant p<=%.2f without correction for SWNorm_btp\n',n_sig, p)
            
            save([savepath date '_resultsANOVAinter_SWNorm_btp.mat'],'pSWNorm_btp', 'resSWNorm_btp');
            clear n_sig res
            clear resGENorm_btp pGENorm_btp resLLNorm_btp pLLNorm_btp resSWNorm_btp pSWNorm_btp
            
            %%%%%%%Sans covarier pour SES%%%
            x = 1;
            resGENorm_btp = [];
            for idtr = 1:numel(itr)
                tr = itr(idtr);
                [pGENorm_btp(:,x),res] = anovan(GENorm_btp(:,idtr),{gr},'varnames',{'Group'});
                close hidden;
                resGENorm_btp = [resGENorm_btp, res];
                x = x + 1;
            end
            
            n_sig = sum(pGENorm_btp(1,:) <= p);
            fprintf('%d Group effects are significant p<=%.2f without correction for GENorm_btp\n',n_sig, p)
            
            save([savepath date '_resultsANOVApascov_GENorm_btp.mat'],'pGENorm_btp', 'resGENorm_btp');
            clear n_sig res
            
            x = 1;
            resLLNorm_btp = [];
            for idtr = 1:numel(itr)
                tr = itr(idtr);
                [pLLNorm_btp(:,x),res] = anovan(LLNorm_btp(:,idtr),{gr},'varnames',{'Group'});
                close hidden;
                resLLNorm_btp = [resLLNorm_btp, res];
                x = x + 1;
            end
            
            n_sig = sum(pLLNorm_btp(1,:) <= p);
            fprintf('%d Group effects are significant p<=%.2f without correction for LLNorm_btp\n',n_sig, p)
            
            save([savepath date '_resultsANOVApascov_LLNorm_btp.mat'],'pLLNorm_btp', 'resLLNorm_btp');
            clear n_sig res
            
            x = 1;
            resSWNorm_btp = [];
            for idtr = 1:numel(itr)-1
                tr = itr(idtr);
                [pSWNorm_btp(:,x),res] = anovan(SWNorm_btp(:,idtr),{gr},'varnames',{'Group'});
                close hidden;
                resSWNorm_btp = [resSWNorm_btp, res];
                x = x + 1;
            end
            
            n_sig = sum(pSWNorm_btp(1,:) <= p);
            fprintf('%d Group effects are significant p<=%.2f without correction for SWNorm_btp\n',n_sig, p)
            
            save([savepath date '_resultsANOVApascov_SWNorm_btp.mat'],'pSWNorm_btp', 'resSWNorm_btp');
            clear n_sig res
            clear resGENorm_btp pGENorm_btp resLLNorm_btp pLLNorm_btp resSWNorm_btp pSWNorm_btp
        end
   end
end
%save([savepath 'workspace.mat'])
