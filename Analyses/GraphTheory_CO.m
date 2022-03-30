%%%%%%%%%%%%%%%%%%%%%%%%%% Graph Theory%%%%%%%%%%%%%%%%%%%%
tic
datapath = 'C:\data\Malnutrition\Resting\NIRS\Analyses\Stats\CORR0,01_0,08\PhysioSatRespEKG\CorrALL\nofisher\';

load ([datapath 'workspace.mat'])
load ([datapath 'workspacemat.mat'])

savepath = [savepath 'GraphTheory\'];
if ~isfolder(savepath)
    mkdir(savepath)
end

%% PARAMETERS%%%%%%%%%%%
calculatemode = 1;
graphmode = 1; %only if calculatemode = 1
analysemode = 1;
normalizemode = 1;

weightedmode = 1;
binarizedmode = 1;
absolutethresh = 1;
proportionalthresh = 1;

%BCT threshold
%Weigthed = garder les valeurs de corr, Binary = remplacer tout par 0 ou 1.%
itra = 0.01:0.01:0.89;% Seuil variable BCT soustraction
itrp = 0.01:0.01:0.89;
p = 0.05;
fixetr = 0.2; %pour les figures

if calculatemode == 1

    %%% Pour Kass, mettre en commentaires%%%%%
    MATall(:,:,54) = [];
    MATallG2(:,:,29) = [];
    idG2(29) = [];
    gr(54) = [];
    ses(54)= [];

    %% Absolute threshold, weighted
    if absolutethresh ==1 && weightedmode ==1 %(WTA)

        disp('Computing absolute threshold & weighted metrics, running')
        LE_wta = zeros(size(MATall,1),size(MATall,3), numel(itra));
        LERR_wta = zeros(size(MATall,1),size(MATall,3), numel(itra));
        GE_wta = zeros(size(MATall,3), numel(itra)); %%erreur dimension
        GERR_wta = zeros(size(MATall,3), numel(itra));
        LL_wta = zeros(size(MATall,3), numel(itra));
        LLRR_wta = zeros(size(MATall,3), numel(itra));
        CC_wta = zeros(size(MATall,1),size(MATall,3), numel(itra));
        CCRR_wta = zeros(size(MATall,1),size(MATall,3), numel(itra));
        RR_wta = zeros(size(MATall,1),size(MATall,2),100);
        for isubject = 42:size(MATall,3)
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

                A = squeeze(MATall(:,:,isubject)); %prendre la matrice de chaque participant%
                A_wta = threshold_absolute(A, itrv); %Thresholder la matrice (weigthed)%
                A_wta(isnan(A_wta)) = 0;  %% Mettre les NAN à 0
                for indrand = 1:100 %% For each subject, 100 random networks that preserved the degree distribution of the original network were estimated%
                    RR_wta (:,:, indrand) = randmio_und(A_wta,1); %% random matrix
                end
                RRmean_wta = mean(RR_wta,3); %% average of the 100 estimated random matrices
                RRmean_wta(isnan(RRmean_wta)) = 0; %Mettre les NAN à 0

                LE_wta(:,isubject, idtr) = efficiency_wei(A_wta,2); %calculer le LE pour chaque canal, de chaque matrice binarisée de chaque participant%
                LERR_wta(:,isubject, idtr) = efficiency_wei(RRmean_wta,2); %calculer le LE pour chaque canal, de chaque matrice binarisée de chaque participant%

                GE_wta(isubject, idtr) = efficiency_wei(A_wta,0); %calculer le GE pour chaque matrice binarisée de chaque participant%
                GERR_wta(isubject, idtr) = efficiency_wei(RRmean_wta,0); %calculer le GE pour chaque matrice binarisée de chaque participant%

                D = distance_wei(weight_conversion(A_wta,'lengths')); %calculer la matrice de distance%
                [lambda,~,~,~,~] = charpath(D,1,0);
                LL_wta(isubject, idtr) = lambda; %calculer le LL pour chaque matrice de distance binarisée de chaque participant%
                DRR = distance_wei(weight_conversion(RRmean_wta,'lengths'));
                [lambdaRR,~,~,~,~] = charpath(DRR,1,0);
                LLRR_wta(isubject, idtr) = lambdaRR; %calculer le LL pour chaque matrice de distance binarisée de chaque participant%

                CC_wta(:,isubject, idtr) = clustering_coef_wu(A_wta); %calculer le CC pour chaque canal, de chaque matrice binarisée de chaque participant%
                CCRR_wta(:,isubject, idtr) = clustering_coef_wu(RRmean_wta); %calculer le CC pour chaque canal, de chaque matrice binarisée de chaque participant%
            end
            fprintf(1,'\n')
        end
        disp('Computing absolute threshold & weighted metrics, done')

        %% Local Efficiency Input Distance matrix, une valeur par canal, par threshold, par participant%
        fprintf('\t Local Efficiency_wta, running\n')
        label = {'LocalEff'; 'wta'};
        
        LEmean_wta = squeeze(mean(LE_wta,1)); %faire la moyenne du LE de tous les canaux, pour chaque matrice de chaque participant%
        LEmeanRR_wta = squeeze(mean(LERR_wta,1)); %average the random LE of all channels, for each matrix of each participant%

        LENorm_wta = LE_wta./LERR_wta;
        LENormmean_wta = squeeze(nanmean(LENorm_wta,1));

        save([savepath date '_LE_wta.mat'],'LE_wta','LEmean_wta','LERR_wta','LEmeanRR_wta','LENorm_wta','LENormmean_wta');

        if graphmode == 1
            Graph_GraphTheory_CO(LEmean_wta, LEmeanRR_wta, LENormmean_wta, itra, fixetr, label, idG1, idG2, savepath)
        end

        fprintf('\t Local Efficiency_wta, done\n')

        %%Global efficiency, une valeur par threshold par participant
        fprintf('\t Global Efficiency_wta, running\n')
        label = {'GlobalEff'; 'wta'};

        GENorm_wta = GE_wta./GERR_wta;
        save([savepath date '_GE_wta.mat'],'GE_wta', 'GENorm_wta','GERR_wta');

        if graphmode == 1
            Graph_GraphTheory_CO(GE_wta, GERR_wta, GENorm_wta, itra, fixetr, label, idG1, idG2, savepath)
        end

        fprintf('\t Global Efficiency_wta, done\n')

        %%Characteristic pathlenght charpath(D)(BU, BD, WU, WD networks), une valeur par threshold par participant.
        fprintf('\t Characteristic path length_wta, running\n')
        label = {'CharPathLength'; 'wta'};

        LLNorm_wta = LL_wta./LLRR_wta;
        save([savepath date '_LL_wta.mat'],'LL_wta','LLRR_wta','LLNorm_wta');

        if graphmode == 1
            Graph_GraphTheory_CO(LL_wta, LLRR_wta, LLNorm_wta, itra, fixetr, label, idG1, idG2, savepath)
        end

        fprintf('\t Characteristic path length_wta, done\n')

        %%Clustering Coefficient clustering_coef_bu(A), une valeur par canal, par threshold par participant
        fprintf('\t Clustering coefficient_wta, running\n')
        label = {'ClustCoeff'; 'wta'};

        CCmean_wta = squeeze(mean(CC_wta,1)); %faire la moyenne du CC de tous les canaux, pour chaque matrice binarisée de chaque participant%
        CCmeanRR_wta = squeeze(mean(CCRR_wta,1)); %%average random CC of all channels, for each binarized matrix of each participant%
        
        CCNorm_wta = CC_wta./CCRR_wta;
        CCNormmean_wta = CCmean_wta./CCmeanRR_wta;
        save([savepath date '_CC_wta.mat'],'CC_wta','CCmean_wta','CCmeanRR_wta','CCNormmean_wta');

        if graphmode == 1
            Graph_GraphTheory_CO(CCmean_wta, CCmeanRR_wta, CCNormmean_wta, itra, fixetr, label, idG1, idG2, savepath)
        end

        fprintf('\t Clustering coefficient_wta, done\n')

        %%Small worldness index%
        fprintf('\t Small Worldness index, running\n')
        label = {'Small Worldness index'; 'wta'};

        SW_wta = CCmean_wta./LL_wta;
        SWNorm_wta = CCNormmean_wta./LLNorm_wta;
        save([savepath date '_SW_wta.mat'],'SW_wta','SWNorm_wta');

        if graphmode == 1
            Graph_GraphTheorySW_CO(SW_wta, SWNorm_wta, itra, fixetr, label, idG1, idG2, savepath)
        end

        fprintf('\t Small Worldness index, done\n')
    end

    %% Absolute threshold, binarized
    if absolutethresh == 1 && binarizedmode ==1 %BTA

        disp('Computing absolute threshold & binarized metrics, running')

        LE_bta = zeros(size(MATall,1),size(MATall,3), numel(itra));
        LERR_bta = zeros(size(MATall,1),size(MATall,3), numel(itra));
        GE_bta = zeros(size(MATall,3), numel(itra));
        GERR_bta = zeros(size(MATall,3), numel(itra));
        LL_bta = zeros(size(MATall,3), numel(itra));
        LLRR_bta = zeros(size(MATall,3), numel(itra));
        CC_bta = zeros(size(MATall,1),size(MATall,3), numel(itra));
        CCRR_bta = zeros(size(MATall,1),size(MATall,3), numel(itra));
        RR_bta = zeros(size(MATall,1),size(MATall,2),100);
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
                %fprintf('\t \t Computing threshold %.2f\n',itrv);
                A = squeeze(MATall(:,:,isubject)); %prendre la matrice de chaque participant%
                A_wta = threshold_absolute(A, itrv); %Thresholder la matrice (weigthed)%
                A_bta = weight_conversion(A_wta,'binarize'); %Binarizer la matrice (binarized)
                A_bta(isnan(A_bta)) = 0;  %% Mettre les NAN à 0
                for indrand = 1:100 %% For each subject, 100 random networks that preserved the degree distribution of the original network were estimated%
                    RR_bta (:,:, indrand) = randomizer_bin_und(A_bta,1); %% random matrix
                end
                RRmean_bta = mean(RR_bta,3); %% average of the 100 estimated random matrices
                RRmean_bta(isnan(RRmean_bta)) = 0; %Mettre les NAN à 0

                LE_bta(:,isubject, idtr) = efficiency_bin(A_bta,1); %calculer le LE pour chaque canal, de chaque matrice binarisée de chaque participant%
                LERR_bta(:,isubject, idtr) = efficiency_bin(RRmean_bta,1); %calculer le LE pour chaque canal, de chaque matrice binarisée de chaque participant%

                GE_bta(isubject, idtr) = efficiency_bin(A_bta); %calculer le GE pour chaque matrice de chaque participant%
                GERR_bta(isubject, idtr) = efficiency_bin(RRmean_bta); %calculer le LE pour chaque canal, de chaque matrice de chaque participant%

                D = distance_bin(A_bta); %calculer la matrice de distance%
                [lambda,~,~,~,~] = charpath(D,1,0);
                LL_bta(isubject, idtr) = lambda; %calculer le LL pour chaque matrice de distance binarisée de chaque participant%
                DRR = distance_bin(RRmean_bta);
                [lambdaRR,~,~,~,~] = charpath(DRR,1,0);
                LLRR_bta(isubject, idtr) = lambdaRR; %calculer le LL pour chaque matrice de distance binarisée de chaque participant%

                CC_bta(:,isubject, idtr) = clustering_coef_bu(A_bta); %calculer le CC pour chaque canal, de chaque matrice de chaque participant%
                CCRR_bta(:,isubject, idtr) = clustering_coef_bu(RRmean_bta); %calculer le LE pour chaque canal, de chaque matrice binarisée de chaque participant%
            end
            fprintf(1,'\n')
        end
        disp('Computing absolute threshold & binarized metrics, done')

        %% Local Efficiency Input Distance matrix, une valeur par canal, par threshold, par participant%
        fprintf(' \t Local Efficiency_bta, running\n');
        label = {'LocalEff'; 'bta'};
        
        %faire la moyenne du LE de tous les canaux, pour chaque matrice de chaque participant%
        LEmean_bta = squeeze(mean(LE_bta,1));
        LEmeanRR_bta = squeeze(mean(LERR_bta,1));

        LENorm_bta = LE_bta./LERR_bta;
        LENormmean_bta = squeeze(nanmean(LENorm_bta,1));

        save([savepath date '_LE_bta.mat'],'LE_bta','LEmean_bta','LERR_bta','LEmeanRR_bta','LENorm_bta','LENormmean_bta');

        if graphmode == 1
            Graph_GraphTheory_CO(LEmean_bta, LEmeanRR_bta, LENormmean_bta, itra, fixetr, label, idG1, idG2, savepath)
        end

        fprintf('\t Local Efficiency_bta, done\n');

        %%Global efficiency, une valeur par threshold par participant
        fprintf('\t Global Efficiency_bta, running\n');
        label = {'GlobalEff'; 'bta'};

        GENorm_bta = GE_bta./GERR_bta;
        save([savepath date '_GE_bta.mat'],'GE_bta', 'GENorm_bta','GERR_bta');

        if graphmode == 1
            Graph_GraphTheory_CO(GE_bta, GERR_bta, GENorm_bta, itra, fixetr, label, idG1, idG2, savepath)
        end

        fprintf('\t Global Efficiency_bta, done\n');

        %%Characteristic pathlenght charpath(D)(BU, BD, WU, WD networks), une valeur par threshold par participant.
        fprintf('\t Characteristic path length_bta, running\n');
        label = {'CharPathLength'; 'bta'};

        LLNorm_bta = LL_bta./LLRR_bta;
        save([savepath date '_LL_bta.mat'],'LL_bta','LLRR_bta','LLNorm_bta');

        if graphmode == 1
            Graph_GraphTheory_CO(LL_bta, LLRR_bta, LLNorm_bta, itra, fixetr, label, idG1, idG2, savepath)
        end

        fprintf('\t Characteristic path length_bta, done\n');

        %%Clustering Coefficient clustering_coef_bu(A), une valeur par canal, par threshold par participant
        fprintf('\t Clustering coefficient_bta, running\n');
        label = {'ClustCoeff'; 'bta'};

        CCmean_bta = squeeze(mean(CC_bta,1)); %faire la moyenne du CC de tous les canaux, pour chaque matrice binarisée de chaque participant%
        CCmeanRR_bta = squeeze(mean(CCRR_bta,1)); %%average random CC of all channels, for each binarized matrix of each participant%

        CCNorm_bta = CC_bta./CCRR_bta;
        CCNormmean_bta = CCmean_bta./CCmeanRR_bta;
        save([savepath date '_CC_bta.mat'],'CC_bta','CCmean_bta','CCmeanRR_bta','CCNormmean_bta');

        if graphmode == 1
            Graph_GraphTheory_CO(CCmean_bta, CCmeanRR_bta, CCNormmean_bta, itra, fixetr, label, idG1, idG2, savepath)
        end

        fprintf('\t Clustering coefficient_bta, done\n');

        %%Small worldness index%
        fprintf('\t Small Worldness index_bta, running\n');
        label = {'Small Worldness index'; 'bta'};

        SW_bta = CCmean_bta./LL_bta;
        SWNorm_bta = CCNormmean_bta./LLNorm_bta;
        save([savepath date '_SW_bta.mat'],'SW_bta','SWNorm_bta');

        if graphmode == 1
            Graph_GraphTheorySW_CO(SW_bta, SWNorm_bta, itra, fixetr, label, idG1, idG2, savepath)
        end

        fprintf('\t Small Worldness index_bta, done\n');
    end

    %% Proportional threshold & weighted
    if proportionalthresh ==1 && weightedmode ==1 %WTP

        disp('Computing proportional threshold & weighted metrics, running')

        LE_wtp = zeros(size(MATall,1),size(MATall,3), numel(itrp));
        LERR_wtp = zeros(size(MATall,1),size(MATall,3), numel(itrp));
        GE_wtp = zeros(size(MATall,3), numel(itrp)); %%erreur dimension
        GERR_wtp = zeros(size(MATall,3), numel(itrp));
        LL_wtp = zeros(size(MATall,3), numel(itrp));
        LLRR_wtp = zeros(size(MATall,3), numel(itrp));
        CC_wtp = zeros(size(MATall,1),size(MATall,3), numel(itrp));
        CCRR_wtp = zeros(size(MATall,1),size(MATall,3), numel(itrp));
        RR_wtp = zeros(size(MATall,1),size(MATall,2),100);
        for isubject = 40:size(MATall,3)
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
                %fprintf('\t \t Computing threshold %.2f\n',itrv);
                A = squeeze(MATall(:,:,isubject)); %prendre la matrice de chaque participant%
                A_wtp = threshold_proportional(A,itrv); %Thresholder la matrice (weigthed)%
                A_wtp(isnan(A_wtp)) = 0;  %% Mettre les NAN à 0
%                 if round(itrv,2) == 0.21 && isubject == 35 %%Quick fix pour le sujet 39/27 dont le threshold 0.17 ne fonctionne pas
%                     RRmean_wtp = nan(length(A_wtp));
%                     RRmean_wtp(isnan(RRmean_wtp)) = 0; %Mettre les NAN à 0
%                 elseif round(itrv,2) == 0.17 && (isubject == 27 || isubject == 39)
%                     RRmean_wtp = nan(length(A_wtp));
%                     RRmean_wtp(isnan(RRmean_wtp)) = 0; %Mettre les NAN à 0
%                 else
                    for indrand = 1:100 %% For each subject, 100 random networks that preserved the degree distribution of the original network were estimated%
                        RR_wtp (:,:, indrand) = randmio_und(A_wtp,1); %% random matrix
                    end
                    RRmean_wtp = mean(RR_wtp,3); %% average of the 100 estimated random matrices
                    RRmean_wtp(isnan(RRmean_wtp)) = 0; %Mettre les NAN à 0
%                 end

                LE_wtp(:,isubject, idtr) = efficiency_wei(A_wtp,2); %calculer le LE pour chaque canal, de chaque matrice binarisée de chaque participant%
                LERR_wtp(:,isubject, idtr) = efficiency_wei(RRmean_wtp,2); %calculer le LE pour chaque canal, de chaque matrice binarisée de chaque participant%

                GE_wtp(isubject, idtr) = efficiency_wei(A_wtp,0); %calculer le GE pour chaque matrice de chaque participant%
                GERR_wtp(isubject, idtr) = efficiency_wei(RRmean_wtp,0); %calculer le LE pour chaque canal, de chaque matrice de chaque participant%

                D = distance_wei(weight_conversion(A_wtp,'lengths')); %calculer la matrice de distance% %%Pas certaine du bon input pour la fonction de distance
                [lambda,~,~,~,~] = charpath(D,1,0);
                LL_wtp(isubject, idtr) = lambda; %calculer le LL pour chaque matrice de distance binarisée de chaque participant%
                DRR = distance_wei(weight_conversion(RRmean_wtp,'lengths'));
                [lambdaRR,~,~,~,~] = charpath(DRR,1,0);
                LLRR_wtp(isubject, idtr) = lambdaRR; %calculer le LL pour chaque matrice de distance binarisée de chaque participant%

                CC_wtp(:,isubject, idtr) = clustering_coef_wu(A_wtp); %calculer le CC pour chaque canal, de chaque matrice de chaque participant%
                CCRR_wtp(:,isubject, idtr) = clustering_coef_wu(RRmean_wtp); %calculer le LE pour chaque canal, de chaque matrice binarisée de chaque participant%
            end
            fprintf(1,'\n')
        end
        disp('Computing proportional threshold & weighted metrics, done')

        %% Local Efficiency Input Distance matrix, une valeur par canal, par threshold, par participant%
        fprintf('\t Local Efficiency_wtp, running\n')
        label = {'LocalEff'; 'wtp'};
        
        LEmean_wtp = squeeze(mean(LE_wtp,1)); %faire la moyenne du LE de tous les canaux, pour chaque matrice de chaque participant%
        LEmeanRR_wtp = squeeze(mean(LERR_wtp,1)); %average the random LE of all channels, for each matrix of each participant%.
        LENorm_wtp = LE_wtp./LERR_wtp;
        LENormmean_wtp = squeeze(nanmean(LENorm_wtp,1));
        save([savepath date '_LE_wtp.mat'],'LE_wtp','LEmean_wtp','LERR_wtp','LEmeanRR_wtp','LENorm_wtp','LENormmean_wtp');

        if graphmode == 1
            Graph_GraphTheory_CO(LEmean_wtp, LEmeanRR_wtp, LENormmean_wtp, itrp, fixetr, label, idG1, idG2, savepath)
        end

        fprintf('\t Local Efficiency_wtp, done\n')

        %%Global efficiency, une valeur par threshold par participant
        fprintf('\t Global Efficiency_wtp, running\n')
        label = {'GlobalEff'; 'wtp'};

        GENorm_wtp = GE_wtp./GERR_wtp;
        save([savepath date '_GE_wtp.mat'],'GE_wtp', 'GENorm_wtp','GERR_wtp');

        if graphmode == 1
            Graph_GraphTheory_CO(GE_wtp, GERR_wtp, GENorm_wtp, itrp, fixetr, label, idG1, idG2, savepath)
        end

        fprintf('\t Global Efficiency_wtp, done\n')

        %%Characteristic pathlenght charpath(D)(BU, BD, WU, WD networks), une valeur par threshold par participant.
        fprintf('\t Characteristic path length_wtp, running\n')
        label = {'CharPathLength'; 'wtp'};

        LLNorm_wtp = LL_wtp./LLRR_wtp;
        save([savepath date '_LL_wtp.mat'],'LL_wtp','LLRR_wtp','LLNorm_wtp');

        if graphmode == 1
            Graph_GraphTheory_CO(LL_wtp, LLRR_wtp, LLNorm_wtp, itrp, fixetr, label, idG1, idG2, savepath)
        end

        fprintf('\tCharacteristic path length_wtp, done\n')

        %%Clustering Coefficient clustering_coef_bu(A), une valeur par canal, par threshold par participant
        fprintf('\t Clustering coefficient_wtp, running\n')
        label = {'ClustCoeff'; 'wtp'};

        CCmean_wtp = squeeze(mean(CC_wtp,1)); %faire la moyenne du CC de tous les canaux, pour chaque matrice binarisée de chaque participant%
        CCmeanRR_wtp = squeeze(mean(CCRR_wtp,1)); %%average random CC of all channels, for each binarized matrix of each participant%
        CCNorm_wtp = CC_wtp./CCRR_wtp;
        CCNormmean_wtp = CCmean_wtp./CCmeanRR_wtp;
        save([savepath date '_CC_wtp.mat'],'CC_wtp','CCmean_wtp','CCmeanRR_wtp','CCNormmean_wtp');

        if graphmode == 1
            Graph_GraphTheory_CO(CCmean_wtp, CCmeanRR_wtp, CCNormmean_wtp, itrp, fixetr, label, idG1, idG2, savepath)
        end

        fprintf('\t Clustering coefficient_wtp, done\n')

        %%Small worldness index%
        fprintf('\t Small Worldness index_wtp, running\n')
        label = {'Small Worldness index'; 'wtp'};

        SW_wtp = CCmean_wtp./LL_wtp;
        SWNorm_wtp = CCNormmean_wtp./LLNorm_wtp;
        save([savepath date '_SW_wtp.mat'],'SW_wtp','SWNorm_wtp');

        if graphmode == 1
            Graph_GraphTheorySW_CO(SW_wtp, SWNorm_wtp, itrp, fixetr, label, idG1, idG2, savepath)
        end

        fprintf('\t Small Worldness index_wtp, done\n')

    end
    %% Proportional threshold & binarized
    if proportionalthresh ==1 && binarizedmode ==1 %BTP

        disp('Computing proportional threshold & binarized metrics, running')

        LE_btp = zeros(size(MATall,1),size(MATall,3), numel(itrp));
        LERR_btp = zeros(size(MATall,1),size(MATall,3), numel(itrp));
        GE_btp = zeros(size(MATall,3), numel(itrp)); %%erreur dimension
        GERR_btp = zeros(size(MATall,3), numel(itrp));
        LL_btp = zeros(size(MATall,3), numel(itrp));
        LLRR_btp = zeros(size(MATall,3), numel(itrp));
        CC_btp = zeros(size(MATall,1),size(MATall,3), numel(itrp));
        CCRR_btp = zeros(size(MATall,1),size(MATall,3), numel(itrp));
        RR_btp = zeros(size(MATall,1),size(MATall,2),100);
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
                %fprintf('\t \t Computing threshold %.2f\n',itrv);
                A = squeeze(MATall(:,:,isubject)); %prendre la matrice de chaque participant%
                A_wtp = threshold_proportional(A,itrv); %Thresholder la matrice (weigthed)%
                A_btp = weight_conversion(A_wtp,'binarize'); %Binarizer la matrice (binarized)
                A_btp(isnan(A_btp)) = 0;  %% Mettre les NAN à 0
                for indrand = 1:100 %% For each subject, 100 random networks that preserved the degree distribution of the original network were estimated%
                    RR_btp (:,:, indrand) = randomizer_bin_und(A_btp,1); %% random matrix
                end
                RRmean_btp = mean(RR_btp,3); %% average of the 100 estimated random matrices
                RRmean_btp(isnan(RRmean_btp)) = 0; %Mettre les NAN à 0

                LE_btp(:,isubject, idtr) = efficiency_bin(A_btp,1); %calculer le LE pour chaque canal, de chaque matrice binarisée de chaque participant%
                LERR_btp(:,isubject, idtr) = efficiency_bin(RRmean_btp,1); %calculer le LE pour chaque canal, de chaque matrice binarisée de chaque participant%

                GE_btp(isubject, idtr) = efficiency_bin(A_btp); %calculer le GE pour chaque matrice de chaque participant%
                GERR_btp(isubject, idtr) = efficiency_bin(RRmean_btp); %calculer le LE pour chaque canal, de chaque matrice de chaque participant%

                D = distance_bin(A_btp); %calculer la matrice de distance%
                [lambda,~,~,~,~] = charpath(D,1,0);
                LL_btp(isubject, idtr) = lambda; %calculer le LL pour chaque matrice de distance binarisée de chaque participant%
                DRR = distance_bin(RRmean_btp);
                [lambdaRR,~,~,~,~] = charpath(DRR,1,0);
                LLRR_btp(isubject, idtr) = lambdaRR; %calculer le LL pour chaque matrice de distance binarisée de chaque participant%

                CC_btp(:,isubject, idtr) = clustering_coef_bu(A_btp); %calculer le CC pour chaque canal, de chaque matrice de chaque participant%
                CCRR_btp(:,isubject, idtr) = clustering_coef_bu(RRmean_btp); %calculer le LE pour chaque canal, de chaque matrice binarisée de chaque participant%
            end
            fprintf(1,'\n')
        end
        disp('Computing proportional threshold & binarized metrics, done')

        %% Local Efficiency Input Distance matrix, une valeur par canal, par threshold, par participant%
        fprintf('\t Local Efficiency_btp, running\n');
        label = {'LocalEff'; 'btp'};

        LEmean_btp = squeeze(mean(LE_btp,1)); %faire la moyenne du LE de tous les canaux, pour chaque matrice de chaque participant%
        LEmeanRR_btp = squeeze(mean(LERR_btp,1)); %average the random LE of all channels, for each matrix of each participant%.
        LENorm_btp = LE_btp./LERR_btp;
        LENormmean_btp = squeeze(nanmean(LENorm_btp,1));
        save([savepath date '_LE_btp.mat'],'LE_btp','LEmean_btp','LERR_btp','LEmeanRR_btp','LENorm_btp','LENormmean_btp');

        if graphmode == 1
            Graph_GraphTheory_CO(LEmean_btp, LEmeanRR_btp, LENormmean_btp, itrp, fixetr, label, idG1, idG2, savepath)
        end

        fprintf('\t Local Efficiency_btp, done\n');

        %%Global efficiency, une valeur par threshold par participant
        fprintf('\t Global Efficiency_btp, running\n');
        label = {'GlobalEff'; 'btp'};

        GENorm_btp = GE_btp./GERR_btp;
        save([savepath date '_GE_btp.mat'],'GE_btp', 'GENorm_btp','GERR_btp');

        if graphmode == 1
            Graph_GraphTheory_CO(GE_btp, GERR_btp, GENorm_btp, itrp, fixetr, label, idG1, idG2, savepath)
        end

        fprintf('\t Global Efficiency_btp, done\n');

        %%Characteristic pathlenght charpath(D)(BU, BD, WU, WD networks), une valeur par threshold par participant.
        fprintf('\t Characteristic path length_btp, running\n');
        label = {'CharPathLength'; 'btp'};

        LLNorm_btp = LL_btp./LLRR_btp;
        save([savepath date '_LL_btp.mat'],'LL_btp','LLRR_btp','LLNorm_btp');

        if graphmode == 1
            Graph_GraphTheory_CO(LL_btp, LLRR_btp, LLNorm_btp, itrp, fixetr, label, idG1, idG2, savepath)
        end

        fprintf('\t Characteristic path length_btp, done\n');

        %%Clustering Coefficient clustering_coef_bu(A), une valeur par canal, par threshold par participant
        fprintf('\t Clustering coefficient_btp, running\n');
        label = {'ClustCoeff'; 'btp'};

        CCmean_btp = squeeze(mean(CC_btp,1)); %faire la moyenne du CC de tous les canaux, pour chaque matrice binarisée de chaque participant%
        CCmeanRR_btp = squeeze(mean(CCRR_btp,1)); %%average random CC of all channels, for each binarized matrix of each participant%
        CCNorm_btp = CC_btp./CCRR_btp;
        CCNormmean_btp = CCmean_btp./CCmeanRR_btp;
        save([savepath date '_CC_btp.mat'],'CC_btp','CCmean_btp','CCmeanRR_btp','CCNormmean_btp');

        if graphmode == 1
            Graph_GraphTheory_CO(CCmean_btp, CCmeanRR_btp, CCNormmean_btp, itrp, fixetr, label, idG1, idG2, savepath)
        end

        fprintf('\t Clustering coefficient_btp, done\n');

        %%Small worldness index%
        fprintf('\t Small Worldness index_btp, running\n');
        label = {'Small Worldness index'; 'btp'};

        SW_btp = CCmean_btp./LL_btp;
        SWNorm_btp = CCNormmean_btp./LLNorm_btp;
        save([savepath date '_SW_btp.mat'],'SW_btp','SWNorm_btp');

        if graphmode == 1
            Graph_GraphTheorySW_CO(SW_btp, SWNorm_btp, itrp, fixetr, label, idG1, idG2, savepath)
        end

        fprintf('\t Small Worldness index_btp, done\n');

    end

    clear A A_bta A_btp A_wta A_wtp b D d1 d1R d2 d2R DRR f id idtr indrand isubject itrv lambda lambdaRR meand1 meand2 tr x y
end

%% ANOVA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if analysemode == 1

    disp(['Computing analyses on ' savepath])

    %%absolute threshold & weighted
    if absolutethresh ==1 && weightedmode ==1 %WTA
        fprintf('Computing analyses for absolute threshold & weighted metrics, running\n');

        label = ('wta');
        ANCOVA_GraphTheory_CO(itra, LEmean_wta, GE_wta, LL_wta, CCmean_wta, SW_wta,...
            LENormmean_wta, GENorm_wta, LLNorm_wta, CCNormmean_wta, SWNorm_wta, p, gr, ses, label, normalizemode, savepath)

        fprintf('Computing analyses for absolute threshold & weighted metrics, done\n\n');

    end

    %%absolute threshold & binarized
    if absolutethresh ==1 && binarizedmode ==1 %BTA
        fprintf('Computing analyses for absolute threshold & binarized metrics, running\n');

        label = ('bta');
        ANCOVA_GraphTheory_CO(itra, LEmean_bta, GE_bta, LL_bta, CCmean_bta, SW_bta,...
            LENormmean_bta, GENorm_bta, LLNorm_bta, CCNormmean_bta, SWNorm_bta, p, gr, ses, label, normalizemode, savepath)

        fprintf('Computing analyses for absolute threshold & binarized metrics, done\n\n');
    end

    %%proportional threshold & weighted
    if proportionalthresh ==1 && weightedmode ==1 %WTP
        fprintf('Computing analyses for proportional threshold & weighted metrics, running\n');

        label = ('wtp');
        ANCOVA_GraphTheory_CO(itrp, LEmean_wtp, GE_wtp, LL_wtp, CCmean_wtp, SW_wtp,...
            LENormmean_wtp, GENorm_wtp, LLNorm_wtp, CCNormmean_wtp, SWNorm_wtp, p, gr, ses, label, normalizemode, savepath)

        fprintf('Computing analyses for proportional threshold & weighted metrics, done\n\n');
    end

    %%proportional threshold & binarized
    if proportionalthresh ==1 && binarizedmode ==1 %BTP
        fprintf('Computing analyses for proportional threshold & binarized metrics, running\n');

        label = ('btp');
        ANCOVA_GraphTheory_CO(itrp, LEmean_btp, GE_btp, LL_btp, CCmean_btp, SW_btp,...
            LENormmean_btp, GENorm_btp, LLNorm_btp, CCNormmean_btp, SWNorm_btp, p, gr, ses, label, normalizemode, savepath)

        fprintf('Computing analyses for proportional threshold & binarized metrics, done\n');
    end
end
disp('Computing analyses, done')
toc
save([savepath date '_workspace.mat'])

%% GRAPHS%%%%%%%%
function Graph_GraphTheory_CO(metric, metricRR, metricNorm, itr, fixetr, label, idG1, idG2, savepath)
try
    figure
    subplot(1,2,1)
    hold on
    plot(itr,mean(metric(idG1,:),1),'r','displayname','Mal')
    plot(itr,mean(metricRR(idG1,:),1),'r--','displayname','MalRndm')
    plot(itr,mean(metric(idG2,:),1),'b','displayname','Ctl')
    plot(itr,mean(metricRR(idG2,:),1),'b--','displayname','CtlRndm')
    xlabel('Threshold','fontsize',12)
    ylabel([label{1} ' ' label{2}],'fontsize',12)
    legend
    box on

    id = sum(itr< fixetr);
    clear d1 d2 d1R d2R
    d1 = metric(idG1,id);
    d2 = metric(idG2,id);
    d1R = metricRR(idG1,id);
    d2R = metricRR(idG2,id);

    subplot(1,2,2)
    hold on
    %b = bar([mean(d1), mean(d1R); mean(d2), mean(d2R)]); %,'FaceColor','flat'
    %b(1).FaceColor = 'b';
    %b(2).FaceColor = 'r';
    b1 = bar(1,[mean(d1), mean(d1R)],'grouped','red');
    b1(2).FaceAlpha = 0.5;
    b2 = bar(2,[mean(d2), mean(d2R)],'grouped','blue');
    b2(2).FaceAlpha = 0.5;
    plot(0.85,d1,'x')
    plot(1.15,d1R,'x')
    plot(1.85,d2,'x')
    plot(2.15,d2R,'x')
    title(['Threshold <', num2str(fixetr)])
    xlabel('Groups','fontsize',12)
    ylabel([label{1} ' ' label{2}],'fontsize',12)
    xlim([0, 3])
    hold off
    legend([b1 b2],{'Mal','MalRndm','Ctl','CtlRndm'})
    box on

    savefig([savepath date '_' label{1} '_' label{2} '.fig'])
    f = gcf;
    exportgraphics(f,[savepath date '_' label{1} '_' label{2} '.png'])
    close

    figure
    subplot(1,2,1)
    hold on
    plot(itr,mean(metricNorm(idG1,:),1),'r','displayname','Mal')
    plot(itr,mean(metricNorm(idG2,:),1),'b','displayname','Ctl')
    xlabel('Threshold','fontsize',12)
    ylabel(['Norm' label{1} ' ' label{2}],'fontsize',12)
    xlim([0, 1])
    legend

    id = sum(itr< fixetr);
    clear d1 d2
    d1 = metricNorm(idG1,id);
    d2 = metricNorm(idG2,id);

    subplot(1,2,2)
    hold on
    bar(1,nanmean(d1),'r')
    bar(2,nanmean(d2),'b')
    plot(1,d1,'x')
    plot(2,d2,'x')
    title(['Threshold <', num2str(fixetr)])
    xlabel('Groups','fontsize',12)
    ylabel(['Norm' label{1} ' ' label{2}],'fontsize',12)
    xlim([0 , 3])
    legend Mal Ctl

    savefig([savepath date '_' label{1} 'Norm_' label{2} '.fig'])
    f=gcf;
    exportgraphics(f,[savepath date '_' label{1} 'Norm_' label{2} '.png'])
    close

catch
    disp('Error during graphic creation')
end
end

function Graph_GraphTheorySW_CO(metric, metricNorm, itr, fixetr, label, idG1, idG2, savepath)
try
    id = sum(itr < fixetr);
    figure
    subplot(1,2,1)
    hold on
    plot(itr,squeeze(nanmean(metric(idG1,:),1)),'r','displayname','Mal')
    plot(itr,squeeze(nanmean(metric(idG2,:),1)),'b','displayname','Ctl')
    xlabel('Threshold','fontsize',12)
    ylabel('Clustering/Path length','fontsize',12)
    title(['Small World Index ' label{2}])
    legend

    subplot(1,2,2)
    hold on
    meand1 = mean(metric(idG1,id));
    meand2 = mean(metric(idG2,id));
    b1 = bar(1,real(meand1),'red');
    b2 = bar(2,real(meand2),'blue');
    xlabel('Groups','fontsize',12)
    ylabel(['Small World Index ' label{2}],'fontsize',12)
    xlim([0 , 3])
    title(['Threshold <', num2str(fixetr)])
    hold off
    box on
    legend Mal Ctl

    savefig([savepath date,'_SmallWorld_' label{2} '.fig'])
    f=gcf;
    exportgraphics(f,[savepath date '_SmallWorld_' label{2} '.png'])
    close

    figure
    subplot(1,2,1)
    hold on
    plot(itr,squeeze(nanmean(metricNorm(idG1,:),1)),'r','displayname','Mal')
    plot(itr,squeeze(nanmean(metricNorm(idG2,:),1)),'b','displayname','Ctl')
    xlabel('Threshold','fontsize',12)
    ylabel('Norm Clustering/Norm Path length','fontsize',12)
    xlim([0, 1])
    title(['Small World Index Norm ' label{2}])
    legend

    subplot(1,2,2)
    hold on
    meand1 = mean(metricNorm(idG1,id));
    meand2 = mean(metricNorm(idG2,id));
    b1 = bar(1,real(meand1),'r');
    b2 = bar(2,real(meand2),'b');
    xlabel('Groups','fontsize',12)
    ylabel(['Small World Index Norm ' label{2}],'fontsize',12)
    xlim([0 , 3])
    title(['Threshold <', num2str(fixetr)])
    hold off
    box on
    legend Mal Ctl

    savefig([savepath date,'_SmallWorldNorm_' label{2} '.fig'])
    f=gcf;
    exportgraphics(f,[savepath date '_SmallWorldNorm_' label{2} '.png'])
    close

catch
    disp('Error during graphic creation')
end
end

%% ANCOVA%%%%%%%%%%
function ANCOVA_GraphTheory_CO(itr, LEmean, GE, LL, CCmean, SW, LENormmean, GENorm, LLNorm, CCNormmean, SWNorm, p, gr, ses, label, normalizemode, savepath)

%%%%%%%Effets principaux de Gr et SSE%%%
fprintf('\t Computing ANOVA with Group & SES\n')
x = 1;
resLE = [];
for idtr = 1:numel(itr)
    tr = itr(idtr);
    [pLE(:,x),res] = anovan(real(LEmean(:,idtr)),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
    close hidden;
    resLE = [resLE, res];
    x = x + 1;
end

n_sig = sum(pLE(1,:) <= p);
fprintf('\t \t %d Group effects are significant p<=%.2f without correction for LEmean_%s\n',n_sig, p, label)

n_sig = sum(pLE(2,:) <= p);
fprintf('\t \t %d SES effects are significant p<=%.2f without correction for LEmean_%s\n',n_sig, p, label)
clear n_sig res

x = 1;
resGE = [];
for idtr = 1:numel(itr)
    tr = itr(idtr);
    [pGE(:,x),res] = anovan(real(GE(:,idtr)),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
    close hidden;
    resGE = [resGE, res];
    x = x + 1;
end

n_sig = sum(pGE(1,:) <= p);
fprintf('\t \t %d Group effects are significant p<=%.2f without correction for GE_%s\n',n_sig, p, label)

n_sig = sum(pGE(2,:) <= p);
fprintf('\t \t %d SES effects are significant p<=%.2f without correction for GE_%s\n',n_sig, p, label)
clear n_sig res

x = 1;
resLL = [];
for idtr = 1:numel(itr)
    tr = itr(idtr);
    [pLL(:,x),res] = anovan(real(LL(:,idtr)),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
    close hidden;
    resLL = [resLL, res];
    x = x + 1;
end

n_sig = sum(pLL(1,:) <= p);
fprintf('\t \t %d Group effects are significant p<=%.2f without correction for LL_%s\n',n_sig, p, label)

n_sig = sum(pLL(2,:) <= p);
fprintf('\t \t %d SES effects are significant p<=%.2f without correction for LL_%s\n',n_sig, p, label)
clear n_sig res

x = 1;
resCC = [];
for idtr = 1:numel(itr)
    tr = itr(idtr);
    [pCC(:,x),res] = anovan(real(CCmean(:,idtr)),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
    close hidden;
    resCC = [resCC, res];
    x = x + 1;
end

n_sig = sum(pCC(1,:) <= p);
fprintf('\t \t %d Group effects are significant p<=%.2f without correction for CCmean_%s\n',n_sig, p, label)

n_sig = sum(pCC(2,:) <= p);
fprintf('\t \t %d SES effects are significant p<=%.2f without correction for CCmean_%s\n',n_sig, p, label)
clear n_sig res

x = 1;
resSW = [];
for idtr = 1:numel(itr)-1
    tr = itr(idtr);
    [pSW(:,x),res] = anovan(real(SW(:,idtr)),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
    close hidden;
    resSW = [resSW, res];
    x = x + 1;
end

n_sig = sum(pSW(1,:) <= p);
fprintf('\t \t %d Group effects are significant p<=%.2f without correction for SW_%s\n',n_sig, p, label)

n_sig = sum(pSW(2,:) <= p);
fprintf('\t \t %d SES effects are significant p<=%.2f without correction for SW_%s\n',n_sig, p, label)
clear n_sig res

save([savepath date '_resultsANOVA_' label '.mat'],'pLE', 'resLE','pGE', 'resGE','pLL', 'resLL','pCC', 'resCC','pSW', 'resSW')
clear resLE pLE resGE pGE resLL pLL resCC pCC resSW pSW

%%%%%%%Ajout de l'interaction%%%
fprintf('\t Computing ANOVA with Group, SES & Group*SES\n')
x = 1;
resLE = [];
for idtr = 1:numel(itr)
    tr = itr(idtr);
    [pLE(:,x),res] = anovan(real(LEmean(:,idtr)),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
    close hidden;
    resLE = [resLE, res];
    x = x + 1;
end

n_sig = sum(pLE(1,:) <= p);
fprintf('\t \t %d Group effects are significant p<=%.2f without correction for LEmean_%s\n',n_sig, p, label)

n_sig = sum(pLE(2,:) <= p);
fprintf('\t \t %d SES effects are significant p<=%.2f without correction for LEmean_%s\n',n_sig, p, label)

n_sig = sum(pLE(3,:) <= p);
fprintf('\t \t %d Group*SES effects are significant p<=%.2f without correction for LEmean_%s\n',n_sig, p, label)
clear n_sig res

x = 1;
resGE = [];
for idtr = 1:numel(itr)
    tr = itr(idtr);
    [pGE(:,x),res] = anovan(real(GE(:,idtr)),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
    close hidden;
    resGE = [resGE, res];
    x = x + 1;
end

n_sig = sum(pGE(1,:) <= p);
fprintf('\t \t %d Group effects are significant p<=%.2f without correction for GE_%s\n',n_sig, p, label)

n_sig = sum(pGE(2,:) <= p);
fprintf('\t \t %d SES effects are significant p<=%.2f without correction for GE_%s\n',n_sig, p, label)

n_sig = sum(pGE(3,:) <= p);
fprintf('\t \t %d Group*SES effects are significant p<=%.2f without correction for GE_%s\n',n_sig, p, label)
clear n_sig res

x = 1;
resLL = [];
for idtr = 1:numel(itr)
    tr = itr(idtr);
    [pLL(:,x),res] = anovan(real(LL(:,idtr)),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
    close hidden;
    resLL = [resLL, res];
    x = x + 1;
end

n_sig = sum(pLL(1,:) <= p);
fprintf('\t \t %d Group effects are significant p<=%.2f without correction for LL_%s\n',n_sig, p, label)

n_sig = sum(pLL(2,:) <= p);
fprintf('\t \t %d SES effects are significant p<=%.2f without correction for LL_%s\n',n_sig, p, label)

n_sig = sum(pLL(3,:) <= p);
fprintf('\t \t %d Group*SES effects are significant p<=%.2f without correction for LL_%s\n',n_sig, p, label)
clear n_sig res

x = 1;
resCC = [];
for idtr = 1:numel(itr)
    tr = itr(idtr);
    [pCC(:,x),res] = anovan(real(CCmean(:,idtr)),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
    close hidden;
    resCC = [resCC, res];
    x = x + 1;
end

n_sig = sum(pCC(1,:) <= p);
fprintf('\t \t %d Group effects are significant p<=%.2f without correction for CCmean_%s\n',n_sig, p, label)

n_sig = sum(pCC(2,:) <= p);
fprintf('\t \t %d SES effects are significant p<=%.2f without correction for CCmean_%s\n',n_sig, p, label)

n_sig = sum(pCC(3,:) <= p);
fprintf('\t \t %d Group*SES effects are significant p<=%.2f without correction for CCmean_%s\n',n_sig, p, label)
clear n_sig res

x = 1;
resSW = [];
for idtr = 1:numel(itr)-1
    tr = itr(idtr);
    [pSW(:,x),res] = anovan(real(SW(:,idtr)),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
    close hidden;
    resSW = [resSW, res];
    x = x + 1;
end

n_sig = sum(pSW(1,:) <= p);
fprintf('\t \t %d Group effects are significant p<=%.2f without correction for SW_%s\n',n_sig, p, label)

n_sig = sum(pSW(2,:) <= p);
fprintf('\t \t %d SES effects are significant p<=%.2f without correction for SW_%s\n',n_sig, p, label)

n_sig = sum(pSW(3,:) <= p);
fprintf('\t \t %d Group*SES effects are significant p<=%.2f without correction for SW_%s\n',n_sig, p, label)
clear n_sig res

save([savepath date '_resultsANOVAinter_' label '.mat'],'pLE', 'resLE','pGE', 'resGE','pLL', 'resLL','pCC', 'resCC','pSW', 'resSW')
clear resLE pLE resGE pGE resLL pLL resCC pCC resSW pSW

%%%%%%%Sans covarier pour SES%%%
fprintf('\t Computing ANOVA with Group ONLY\n')
x = 1;
resLE = [];
for idtr = 1:numel(itr)
    tr = itr(idtr);
    [pLE(:,x),res] = anovan(real(LEmean(:,idtr)),{gr},'varnames',{'Group'});
    close hidden;
    resLE = [resLE, res];
    x = x + 1;
end

n_sig = sum(pLE(1,:) <= p);
fprintf('\t \t %d Group effects are significant p<=%.2f without correction for LEmean_%s\n',n_sig, p, label)
clear n_sig res

x = 1;
resGE = [];
for idtr = 1:numel(itr)
    tr = itr(idtr);
    [pGE(:,x),res] = anovan(real(GE(:,idtr)),{gr},'varnames',{'Group'});
    close hidden;
    resGE = [resGE, res];
    x = x + 1;
end

n_sig = sum(pGE(1,:) <= p);
fprintf('\t \t %d Group effects are significant p<=%.2f without correction for GE_%s\n',n_sig, p, label)
clear n_sig res

x = 1;
resLL = [];
for idtr = 1:numel(itr)
    tr = itr(idtr);
    [pLL(:,x),res] = anovan(real(LL(:,idtr)),{gr},'varnames',{'Group'});
    close hidden;
    resLL = [resLL, res];
    x = x + 1;
end

n_sig = sum(pLL(1,:) <= p);
fprintf('\t \t %d Group effects are significant p<=%.2f without correction for LL_%s\n',n_sig, p, label)
clear n_sig res

x = 1;
resCC = [];
for idtr = 1:numel(itr)
    tr = itr(idtr);
    [pCC(:,x),res] = anovan(real(CCmean(:,idtr)),{gr},'varnames',{'Group'});
    close hidden;
    resCC = [resCC, res];
    x = x + 1;
end

n_sig = sum(pCC(1,:) <= p);
fprintf('\t \t %d Group effects are significant p<=%.2f without correction for CCmean_%s\n',n_sig, p, label)
clear n_sig res

x = 1;
resSW = [];
for idtr = 1:numel(itr)-1
    tr = itr(idtr);
    [pSW(:,x),res] = anovan(real(SW(:,idtr)),{gr},'varnames',{'Group'});
    close hidden;
    resSW = [resSW, res];
    x = x + 1;
end

n_sig = sum(pSW(1,:) <= p);
fprintf('\t \t %d Group effects are significant p<=%.2f without correction for SW_%s\n\n',n_sig, p, label)
clear n_sig res

save([savepath date '_resultsANOVApascov_' label '.mat'],'pLE', 'resLE','pGE', 'resGE','pLL', 'resLL','pCC', 'resCC','pSW', 'resSW')
clear resLE pLE resGE pGE resLL pLL resCC pCC resSW pSW

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if normalizemode == 1
    %%%%%%%Effets principaux de Gr et SSE%%%
    fprintf('\t Computing ANOVA with Group & SES\n')
    x = 1;
    resLENorm = [];
    for idtr = 1:numel(itr)
        tr = itr(idtr);
        [pLENorm(:,x),res] = anovan(real(LENormmean(:,idtr)),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
        close hidden;
        resLENorm = [resLENorm, res];
        x = x + 1;
    end

    n_sig = sum(pLENorm(1,:) <= p);
    fprintf('\t \t %d Group effects are significant p<=%.2f without correction for LENormmean_%s\n',n_sig, p, label)

    n_sig = sum(pLENorm(2,:) <= p);
    fprintf('\t \t %d SES effects are significant p<=%.2f without correction for LENormmean_%s\n',n_sig, p, label)
    clear n_sig res

    x = 1;
    resGENorm = [];
    for idtr = 1:numel(itr)
        tr = itr(idtr);
        [pGENorm(:,x),res] = anovan(real(GENorm(:,idtr)),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
        close hidden;
        resGENorm = [resGENorm, res];
        x = x + 1;
    end

    n_sig = sum(pGENorm(1,:) <= p);
    fprintf('\t \t %d Group effects are significant p<=%.2f without correction for GENorm_%s\n',n_sig, p, label)

    n_sig = sum(pGENorm(2,:) <= p);
    fprintf('\t \t %d SES effects are significant p<=%.2f without correction for GENorm_%s\n',n_sig, p, label)
    clear n_sig res

    x = 1;
    resLLNorm = [];
    for idtr = 1:numel(itr)
        tr = itr(idtr);
        [pLLNorm(:,x),res] = anovan(real(LLNorm(:,idtr)),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
        close hidden;
        resLLNorm = [resLLNorm, res];
        x = x + 1;
    end

    n_sig = sum(pLLNorm(1,:) <= p);
    fprintf('\t \t %d Group effects are significant p<=%.2f without correction for LLNorm_%s\n',n_sig, p, label)

    n_sig = sum(pLLNorm(2,:) <= p);
    fprintf('\t \t %d SES effects are significant p<=%.2f without correction for LLNorm_%s\n',n_sig, p, label)
    clear n_sig res

    x = 1;
    resCCNorm = [];
    for idtr = 1:numel(itr)
        tr = itr(idtr);
        [pCCNorm(:,x),res] = anovan(real(CCNormmean(:,idtr)),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
        close hidden;
        resCCNorm = [resCCNorm, res];
        x = x + 1;
    end

    n_sig = sum(pCCNorm(1,:) <= p);
    fprintf('\t \t %d Group effects are significant p<=%.2f without correction for CCNormmean_%s\n',n_sig, p, label)

    n_sig = sum(pCCNorm(2,:) <= p);
    fprintf('\t \t %d SES effects are significant p<=%.2f without correction for CCNormmean_%s\n',n_sig, p, label)
    clear n_sig res

    x = 1;
    resSWNorm = [];
    for idtr = 1:numel(itr)-1
        tr = itr(idtr);
        [pSWNorm(:,x),res] = anovan(real(SWNorm(:,idtr)),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
        close hidden;
        resSWNorm = [resSWNorm, res];
        x = x + 1;
    end

    n_sig = sum(pSWNorm(1,:) <= p);
    fprintf('\t \t %d Group effects are significant p<=%.2f without correction for SWNorm_%s\n',n_sig, p, label)

    n_sig = sum(pSWNorm(2,:) <= p);
    fprintf('\t \t %d SES effects are significant p<=%.2f without correction for SWNorm_%s\n',n_sig, p, label)
    clear n_sig res

    save([savepath date '_resultsANOVANORM_' label '.mat'],'pLENorm', 'resLENorm','pGENorm', 'resGENorm','pLLNorm', 'resLLNorm','pCCNorm', 'resCCNorm','pSWNorm', 'resSWNorm')
    clear resLENorm pLENorm resGENorm pGENorm resLLNorm pLLNorm resCCNorm pCCNorm resSWNorm pSWNorm

    %%%%%%%Ajout de l'interaction%%%
    fprintf('\t Computing ANOVA with Group, SES & Group*SES\n')
    x = 1;
    resLENorm = [];
    for idtr = 1:numel(itr)
        tr = itr(idtr);
        [pLENorm(:,x),res] = anovan(real(LENormmean(:,idtr)),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
        close hidden;
        resLENorm = [resLENorm, res];
        x = x + 1;
    end

    n_sig = sum(pLENorm(1,:) <= p);
    fprintf('\t \t %d Group effects are significant p<=%.2f without correction for LENormmean_%s\n',n_sig, p, label)

    n_sig = sum(pLENorm(2,:) <= p);
    fprintf('\t \t %d SES effects are significant p<=%.2f without correction for LENormmean_%s\n',n_sig, p, label)

    n_sig = sum(pLENorm(3,:) <= p);
    fprintf('\t \t %d Group*SES effects are significant p<=%.2f without correction for LENormmean_%s\n',n_sig, p, label)
    clear n_sig res

    x = 1;
    resGENorm = [];
    for idtr = 1:numel(itr)
        tr = itr(idtr);
        [pGENorm(:,x),res] = anovan(real(GENorm(:,idtr)),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
        close hidden;
        resGENorm = [resGENorm, res];
        x = x + 1;
    end

    n_sig = sum(pGENorm(1,:) <= p);
    fprintf('\t \t %d Group effects are significant p<=%.2f without correction for GENorm_%s\n',n_sig, p, label)

    n_sig = sum(pGENorm(2,:) <= p);
    fprintf('\t \t %d SES effects are significant p<=%.2f without correction for GENorm_%s\n',n_sig, p, label)

    n_sig = sum(pGENorm(3,:) <= p);
    fprintf('\t \t %d Group*SES effects are significant p<=%.2f without correction for GENorm_%s\n',n_sig, p, label)
    clear n_sig res

    x = 1;
    resLLNorm = [];
    for idtr = 1:numel(itr)
        tr = itr(idtr);
        [pLLNorm(:,x),res] = anovan(real(LLNorm(:,idtr)),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
        close hidden;
        resLLNorm = [resLLNorm, res];
        x = x + 1;
    end

    n_sig = sum(pLLNorm(1,:) <= p);
    fprintf('\t \t %d Group effects are significant p<=%.2f without correction for LLNorm_%s\n',n_sig, p, label)

    n_sig = sum(pLLNorm(2,:) <= p);
    fprintf('\t \t %d SES effects are significant p<=%.2f without correction for LLNorm_%s\n',n_sig, p, label)

    n_sig = sum(pLLNorm(3,:) <= p);
    fprintf('\t \t %d Group*SES effects are significant p<=%.2f without correction for LLNorm_%s\n',n_sig, p, label)
    clear n_sig res

    x = 1;
    resCCNorm = [];
    for idtr = 1:numel(itr)
        tr = itr(idtr);
        [pCCNorm(:,x),res] = anovan(real(CCNormmean(:,idtr)),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
        close hidden;
        resCCNorm = [resCCNorm, res];
        x = x + 1;
    end

    n_sig = sum(pCCNorm(1,:) <= p);
    fprintf('\t \t %d Group effects are significant p<=%.2f without correction for CCNormmean_%s\n',n_sig, p, label)

    n_sig = sum(pCCNorm(2,:) <= p);
    fprintf('\t \t %d SES effects are significant p<=%.2f without correction for CCNormmean_%s\n',n_sig, p, label)

    n_sig = sum(pCCNorm(3,:) <= p);
    fprintf('\t \t %d Group*SES effects are significant p<=%.2f without correction for CCNormmean_%s\n',n_sig, p, label)
    clear n_sig res

    x = 1;
    resSWNorm = [];
    for idtr = 1:numel(itr)-1
        tr = itr(idtr);
        [pSWNorm(:,x),res] = anovan(real(SWNorm(:,idtr)),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
        close hidden;
        resSWNorm = [resSWNorm, res];
        x = x + 1;
    end

    n_sig = sum(pSWNorm(1,:) <= p);
    fprintf('\t \t %d Group effects are significant p<=%.2f without correction for SWNorm_%s\n',n_sig, p, label)

    n_sig = sum(pSWNorm(2,:) <= p);
    fprintf('\t \t %d SES effects are significant p<=%.2f without correction for SWNorm_%s\n',n_sig, p, label)

    n_sig = sum(pSWNorm(3,:) <= p);
    fprintf('\t \t %d Group*SES effects are significant p<=%.2f without correction for SWNorm_%s\n',n_sig, p, label)
    clear n_sig res

    save([savepath date '_resultsANOVAinterNORM_' label '.mat'],'pLENorm', 'resLENorm','pGENorm', 'resGENorm','pLLNorm', 'resLLNorm','pCCNorm', 'resCCNorm','pSWNorm', 'resSWNorm')
    clear resLENorm pLENorm resGENorm pGENorm resLLNorm pLLNorm resCCNorm pCCNorm resSWNorm pSWNorm

    %%%%%%%Sans covarier pour SES%%%
    fprintf('\t Computing ANOVA with Group ONLY\n')
    x = 1;
    resLENorm = [];
    for idtr = 1:numel(itr)
        tr = itr(idtr);
        [pLENorm(:,x),res] = anovan(real(LENormmean(:,idtr)),{gr},'varnames',{'Group'});
        close hidden;
        resLENorm = [resLENorm, res];
        x = x + 1;
    end

    n_sig = sum(pLENorm(1,:) <= p);
    fprintf('\t \t %d Group effects are significant p<=%.2f without correction for LENormmean_%s\n',n_sig, p, label)
    clear n_sig res

    x = 1;
    resGENorm = [];
    for idtr = 1:numel(itr)
        tr = itr(idtr);
        [pGENorm(:,x),res] = anovan(real(GENorm(:,idtr)),{gr},'varnames',{'Group'});
        close hidden;
        resGENorm = [resGENorm, res];
        x = x + 1;
    end

    n_sig = sum(pGENorm(1,:) <= p);
    fprintf('\t \t %d Group effects are significant p<=%.2f without correction for GENorm_%s\n',n_sig, p, label)
    clear n_sig res

    x = 1;
    resLLNorm = [];
    for idtr = 1:numel(itr)
        tr = itr(idtr);
        [pLLNorm(:,x),res] = anovan(real(LLNorm(:,idtr)),{gr},'varnames',{'Group'});
        close hidden;
        resLLNorm = [resLLNorm, res];
        x = x + 1;
    end

    n_sig = sum(pLLNorm(1,:) <= p);
    fprintf('\t \t %d Group effects are significant p<=%.2f without correction for LLNorm_%s\n',n_sig, p, label)
    clear n_sig res

    x = 1;
    resCCNorm = [];
    for idtr = 1:numel(itr)
        tr = itr(idtr);
        [pCCNorm(:,x),res] = anovan(real(CCNormmean(:,idtr)),{gr},'varnames',{'Group'});
        close hidden;
        resCCNorm = [resCCNorm, res];
        x = x + 1;
    end

    n_sig = sum(pCCNorm(1,:) <= p);
    fprintf('\t \t %d Group effects are significant p<=%.2f without correction for CCNormmean_%s\n',n_sig, p, label)
    clear n_sig res

    x = 1;
    resSWNorm = [];
    for idtr = 1:numel(itr)-1
        tr = itr(idtr);
        [pSWNorm(:,x),res] = anovan(real(SWNorm(:,idtr)),{gr},'varnames',{'Group'});
        close hidden;
        resSWNorm = [resSWNorm, res];
        x = x + 1;
    end

    n_sig = sum(pSWNorm(1,:) <= p);
    fprintf('\t \t %d Group effects are significant p<=%.2f without correction for SWNorm_%s\n',n_sig, p, label)
    clear n_sig res

    save([savepath date '_resultsANOVApascovNORM_' label '.mat'],'pLENorm', 'resLENorm','pGENorm', 'resGENorm','pLLNorm', 'resLLNorm','pCCNorm', 'resCCNorm','pSWNorm', 'resSWNorm')
    clear resLENorm pLENorm resGENorm pGENorm resLLNorm pLLNorm resCCNorm pCCNorm resSWNorm pSWNorm
end
end
