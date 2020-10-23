%%%%%%%%PCA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Computing PCA')
tic

datapath = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Stats\CORR\0,01_0,08\';
load ([datapath 'workspace.mat'])
% load ([data 'workspacemat.mat'])
savepath='C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Stats\CORR\0,01_0,08\PCA\Round3\';
if ~isdir(savepath)
    mkdir(savepath)
end
    
channelmode = 1;
roimode = 1;
initialmode = 0;
imputemode = 1;
rotationmode = 1;
norotationmode = 1;
tablegraphmode = 1;
recalculatemode = 0;

p = 0.07;

%%
%%%%%%%%%%%%%%%%%Canaux%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if channelmode == 1
    datach = table2array(tblch(1:end-1,6:end));
    crdatach = table2array(tblch(1:end-1,6:end)) - repmat(nanmean(table2array(tblch(1:end-1,6:end))),size(tblch(1:end-1,6:end),1),1);
    
    if initialmode == 1
        if norotationmode == 1
            tic
            disp('Computing PCA on Initial Channels Data, Unrotated Components')
            
            [coeff_pairch,score_pairch,latent_pairch,tsquared_pairch,explained_pairch,mu_pairch] = pca(datach,'Rows','pairwise');
            opt = statset('pca'); opt.MaxIter = 50000;
            [coeff_alsch,score_alsch,latent_alsch,tsquared_alsch,explained_alsch,mu_alsch] = pca(datach,'Algorithm','als','Options',opt); %trèslong!!
            
            save([savepath date '_resultsPCA_PairCh.mat'],'coeff_pairch', 'score_pairch', 'latent_pairch', 'tsquared_pairch', 'explained_pairch', 'mu_pairch');
            save([savepath date '_resultsPCA_AlsCh.mat'],'coeff_alsch', 'score_alsch', 'latent_alsch', 'tsquared_alsch', 'explained_alsch', 'mu_alsch');
            
            %plot explained variance of each component
            figure('Name','Explained Variance Pairwise Channels')
            hold on
            pareto(explained_pairch)
            yline(70)
            yline(80)
            xlabel('Principal Component')
            ylabel('Variance Explained (%)')
            title('Variance of each component Pairwise Channels')
            savefig([savepath date 'ExpVarPairCh']);
            
            figure('Name','Explained Variance ALS Channels')
            hold on
            pareto(explained_alsch)
            yline(70)
            yline(80)
            xlabel('Principal Component')
            ylabel('Variance Explained (%)')
            title('Variance of each component ALS Channels')
            savefig([savepath date 'ExpVarALSCh']);
            
            %Pair: 70% variance expliquée = 8 composantes, 80% = 12 composantes
            %ALS: 70% variance expliquée = 5 composantes, 80% = 8 composantes
            
            if tablegraphmode == 1
                for c = 1:12
                    name = sprintf('C%.0f',c);
                    figure('Name','Score PairCh')
                    plot(score_pairch(idG1(1:end-1,1),c),1,'r+',score_pairch(idG2(1:end-1,1),c),1,'bo')
                    label = [name ' Principal Component'];
                    xlabel(label);
                    titre = ['Score on ' name ' PairCh'];
                    title(titre)
                    savefig([savepath date 'PairCh ' name]);
                    exportgraphics(gcf,[savepath date 'PairCh ' name '.png'])
                end
                
                for c = 1:8
                    name = sprintf('C%.0f',c);
                    figure('Name','Score ALSCh')
                    plot(score_alsch(idG1(1:end-1,1),c),1,'r+',score_alsch(idG2(1:end-1,1),c),1,'bo')
                    label = [name ' Principal Component'];
                    xlabel(label);
                    titre = ['Score on ' name ' ALSCh'];
                    title(titre)
                    savefig([savepath date 'ALSCh ' name]);
                    exportgraphics(gcf,[savepath date 'ALSCh ' name '.png'])
                end
            end
            
            % estimcrdatach_pairch = score_pairch*coeff_pairch';
            % estimcrdatach_alsch = score_alsch*coeff_alsch';
            
            %recalculate pca with appropriate number of components
            if recalculatemode == 1;
                [coeff_pair8ch,score_pair8ch,latent_pair8ch,tsquared_pair8ch,explained_pair8ch,mu_pair8ch] = pca(datach,'Rows','pairwise','NumComponents',8);
                [coeff_als5ch,score_als5ch,latent_als5ch,tsquared_als5ch,explained_als5ch,mu_als5ch] = pca(datach,'Algorithm','als','NumComponents',5); %trèslong!!
                
                save([savepath date '_resultsPCA_Pair8Ch.mat'],'coeff_pair8ch', 'score_pair8ch', 'latent_pair8ch', 'tsquared_pair8ch', 'explained_pair8ch', 'mu_pair8ch');
                save([savepath date '_resultsPCA_Als5Ch.mat'],'coeff_als5ch', 'score_als5ch', 'latent_als5ch', 'tsquared_als5ch', 'explained_als5ch', 'mu_als5ch');
                
                [coeff_pair12ch,score_pair12ch,latent_pair12ch,tsquared_pair12ch,explained_pair12ch,mu_pair12ch] = pca(datach,'Rows','pairwise','NumComponents',12);
                [coeff_als8ch,score_als8ch,latent_als8ch,tsquared_als8ch,explained_als8ch,mu_als8ch] = pca(datach,'Algorithm','als','NumComponents',8); %trèslong!!
                
                save([savepath date '_resultsPCA_Pair12Ch.mat'],'coeff_pair12ch', 'score_pair12ch', 'latent_pair12ch', 'tsquared_pair12ch', 'explained_pair12ch', 'mu_pair12ch');
                save([savepath date '_resultsPCA_Als8Ch.mat'],'coeff_als8ch', 'score_als8ch', 'latent_als8ch', 'tsquared_als8ch', 'explained_als8ch', 'mu_als8ch');
                
                if tablegraphmode == 1
                    for c = 1:5
                        name = sprintf('C%.0f',c);
                        figure('Name','Score ALS5Ch')
                        plot(score_als5ch(idG1(1:end-1,1),c),1,'r+',score_als5ch(idG2(1:end-1,1),c),1,'bo')
                        label = [name ' Principal Component'];
                        xlabel(label);
                        titre = ['Score on ' name ' ALS5Ch'];
                        title(titre)
                        %savefig([savepath date 'ALS5Ch ' name]);
                        %exportgraphics(gcf,[savepath date 'ALS5Ch ' name '.png'])
                    end
                    
                    for c = 1:8
                        name = sprintf('C%.0f',c);
                        figure('Name','Score ALS8Ch')
                        plot(score_als8ch(idG1(1:end-1,1),c),1,'r+',score_als8ch(idG2(1:end-1,1),c),1,'bo')
                        label = [name ' Principal Component'];
                        xlabel(label);
                        titre = ['Score on ' name ' ALS8Ch'];
                        title(titre)
                        %savefig([savepath date 'ALS8Ch ' name]);
                        %exportgraphics(gcf,[savepath date 'ALS8Ch ' name '.png'])
                    end
                end
            end
            toc
            
%%%%%%%%%%%%%%%%%%ANCOVA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp('Computing ANCOVA on Initial Channels Data, Unrotated Components')
            tic
            
            %%% Pairwise non rotated%%%%%%%%%%%%
            % Trop de NAN dans les scores, seulement une vingaine de participant sans aucun nan%
            % x = 1;
            % ppcapairch = zeros (2,12);
            % respcapairch = [];
            % for c = 1:12
            %     [ppcapairch(:,x),res] = anovan(score_pairch(:,c),{gr(1:end-1,1),ses(1:end-1,1)},'Continuous',2,'varnames',{'Group','SES'});
            %     close hidden;
            %     respcapairch = [respcapairch, res];
            %     x = x + 1;
            % end
            %
            % n_sig = sum(ppcapair(1,:) <= p);
            % fprintf('%d tests are significant at p<=%.2f for ANCOVA on Pairwise non rotated Initial Channels components\n',n_sig, p)
            %
            % clear r x res n_sig
            
            %%% ALS non rotated%%%%%%%%%%%%
            x = 1;
            p_pcaalsch = zeros (2,8);
            res_pcaalsch = [];
            for c = 1:8
                [p_pcaalsch(:,x),res] = anovan(score_alsch(:,c),{gr(1:end-1,1),ses(1:end-1,1)},'Continuous',2,'varnames',{'Group','SES'});
                close hidden;
                res_pcaalsch = [res_pcaalsch, res];
                x = x + 1;
            end
            
            n_sig = sum(p_pcaalsch(1,:) <= p);
            fprintf('%d tests are significant p<=%.2f for ANCOVA on ALS non rotated Initial Channels components\n',n_sig, p)
            
            clear r x res n_sig
            
            if recalculatemode == 1
                %%Seulement ALS, car Pair a trop de NAN
                %%%ALS non rotated recalculé avec 8 composantes%%%%%%%%
                x = 1;
                ppcaals8 = zeros (2,8);
                respcaals8 = [];
                for c = 1:8
                    [ppcaals8(:,x),res] = anovan(score_als8ch(:,c),{gr(1:end-1,1),ses(1:end-1,1)},'Continuous',2,'varnames',{'Group','SES'});
                    close hidden;
                    respcaals8 = [respcaals8, res];
                    x = x + 1;
                end
                
                n_sig = sum(ppcaals8(1,:) <= p);
                fprintf('%d tests are significant p<=%.2f for ANCOVA on ALS8 non rotated Recalculated Channels components\n',n_sig, p)
                
                clear r x res n_sig
                
                %%%ALS non rotated recalculé avec 5 composantes%%%%%%%%
                x = 1;
                ppcaals5 = zeros (2,5);
                respcaals5 = [];
                for c = 1:5
                    [ppcaals5(:,x),res] = anovan(score_als5ch(:,c),{gr(1:end-1,1),ses(1:end-1,1)},'Continuous',2,'varnames',{'Group','SES'});
                    close hidden;
                    respcaals5 = [respcaals5, res];
                    x = x + 1;
                end
                
                n_sig = sum(ppcaals5(1,:) <= p);
                fprintf('%d tests are significant p<=%.2f for ANCOVA on ALS5 non rotated Recalculated Channels components\n',n_sig, p)
                
                clear r x res n_sig
            end
        end
            toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%Rotation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if rotationmode == 1
            disp('Computing PCA on Initial Channels Data, Rotated Components')
            tic
            [coeff_pairchrot, matrix_pairchrot] = rotatefactors(coeff_pairch, 'Method', 'promax','Maxit',1000);
            [coeff_alschrot, matrix_alschrot] = rotatefactors(coeff_alsch, 'Method', 'promax','Maxit',1000);

            score_pairchrot = crdatach/coeff_pairchrot';
            score_alschrot = crdatach/coeff_alschrot';

            latent_pairchrot = sort(nanvar(score_pairchrot),'descend')';
            latent_alschrot = sort(nanvar(score_alschrot),'descend')';

            explained_pairchrot = (latent_pairchrot*100/sum(latent_pairchrot));
            explained_alschrot = (latent_alschrot*100/sum(latent_alschrot));

            %Problème avec ALS, comme je recalcule avec les données originales = NAN%%

            save([savepath date '_resultsPCA_PairChRot.mat'],'coeff_pairchrot', 'score_pairchrot', 'latent_pairchrot', 'explained_pairchrot');
            save([savepath date '_resultsPCA_AlsChRot.mat'],'coeff_alschrot', 'score_alschrot', 'latent_alschrot', 'explained_alschrot');

            if tablegraphmode == 1
                for c = 1:6
                    name = sprintf('C%.0f',c);
                    figure('Name','Score Pairwise ChRot')
                    plot(score_pairchrot(idG1(1:end-1,1),c),1,'r+',score_pairchrot(idG2(1:end-1,1),c),1,'bo')
                    label = [name ' Principal Component'];
                    xlabel(label);
                    titre = ['Score on ' name ' Pairwise ChRot'];
                    title(titre)
                    savefig([savepath date 'PairChRot ' name]);
                    exportgraphics(gcf,[savepath date 'PairChRot ' name '.png'])
                end

                clear label titre c name

                for c = 1:6
                    name = sprintf('C%.0f',c);
                    figure('Name','Score ALS ChRot')
                    plot(score_alschrot(idG1(1:end-1,1),c),1,'r+',score_alschrot(idG2(1:end-1,1),c),1,'bo')
                    label = [name ' Principal Component'];
                    xlabel(label);
                    titre = ['Score on ' name ' ALS ChRot'];
                    title(titre)
                    savefig([savepath date 'ALSChRot ' name]);
                    exportgraphics(gcf,[savepath date 'ALSChRot ' name '.png'])
                end

                clear label titre c name

                %sort tsquared and extract extreme values
                %[st2_pairch,index_pairch] = sort(tsquared_pairch,'descend'); % sort in descending order
                %extreme_pairch = index_pairch(1);
                %labelch(:,extreme_pairch)

                %[st2_alsch,index_alsch] = sort(tsquared_alsch,'descend'); % sort in descending order
                %extreme_alsch = index_alsch(1);
                %labelch(:,extreme_alsch)
            end
            toc

%%%%%%%%%%%%%%%%%%ANCOVA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp('Computing ANCOVA on Initial Channels Data, Rotated Components - Too many NANs')
            tic
            p = 0.07;
            %%Pas la peine de faire des ANOVA sur les données rotatées, car trop de NAN avec les deux méthodes%%%

            %%%Pairwise rotated
        %     x = 1;
        %     p_pcapairchrot = zeros (2,12);
        %     res_pcapairchrot = [];
        %     for c = 1:12
        %         [p_pcapairchrot(:,x),res] = anovan(score_pairchrot(:,c),{gr(1:end-1,1),ses(1:end-1,1)},'Continuous',2,'varnames',{'Group','SES'});
        %         close hidden;
        %         res_pcapairchrot = [res_pcapairchrot, res];
        %         x = x + 1;
        %     end
        %     
        %     n_sig = sum(p_pcapairchrot(1,:) <= p);
        %     fprintf('%d tests are significant p<=%.2f for ANCOVA on Pairwise Rotated Initial Channels components\n',n_sig, p)
        %     
        %     clear r x res n_sig
        %     
        %     %%% ALS rotated%%%%%%%%%%%%
        %     x = 1;
        %     p_pcaalschrot = zeros (2,8);
        %     res_pcaalschrot = [];
        %     for c = 1:8
        %         [p_pcaalschrot(:,x),res] = anovan(score_alschrot(:,c),{gr(1:end-1,1),ses(1:end-1,1)},'Continuous',2,'varnames',{'Group','SES'});
        %         close hidden;
        %         res_pcaalschrot = [res_pcaalschrot, res];
        %         x = x + 1;
        %     end
        %     
        %     n_sig = sum(p_pcaalschrot(1,:) <= p);
        %     fprintf('%d tests are significant p<=%.2f for ANCOVA on ALS Rotated Initial Channels components\n',n_sig, p)
        %     
        %     clear r x res n_sig
        end
    end
    toc
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%PCA sur données imputées%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if imputemode == 1
        if norotationmode == 1
            disp('Computing PCA on Imputed Channels Data, Unrotated Components')
            tic
            
            %%Imputations%%%%%%%%
            
            %load ('MDItoolbox_resultsch.mat')
            
            %varargout = NumberComponents(varargin);
            %expvarch = NumberComponents(datach);
            
            %%On ne voit que les 10 premières composantes, devrait être de 8/12 comme établi plus haut?%%%%
            % 70% de variance = 11 composantes; 80% de variance = 19
            
            %[impdata,impmean,impcov,recdata] = pcambda(data,Nchains,Lchains,Ncomp);
            %[impdatach,impmeanch,impcovch,recdatach] = pcambda(datach,10,100,10); %matrice not positive definite%
            
            %[impdata,impmean,impcov,It,diff,recdata]=pcambia(data,Ncomp,maxiter,tolerance)
            [dataimpch,meanimpch,covimpch,Itch,diffch,recdatach] = pcambia(datach,11,10000,10E-10);
            
            %knnimpute
            
            crdataimpch = dataimpch - repmat(nanmean(dataimpch),size(dataimpch,1),1);
            
            %%PCA%%%%%%%
            [coeff_pairimpch,score_pairimpch,latent_pairimpch,tsquared_pairimpch,explained_pairimpch,mu_pairimpch] = pca(dataimpch,'Rows','pairwise');
            opt = statset('pca'); opt.MaxIter = 10000;
            [coeff_alsimpch,score_alsimpch,latent_alsimpch,tsquared_alsimpch,explained_alsimpch,mu_alsimpch] = pca(dataimpch,'Algorithm','als','Options',opt);
            
            save([savepath date '_resultsPCA_PairImpCh.mat'],'coeff_pairimpch', 'score_pairimpch', 'latent_pairimpch', 'tsquared_pairimpch', 'explained_pairimpch', 'mu_pairimpch');
            save([savepath date '_resultsPCA_AlsImpCh.mat'],'coeff_alsimpch', 'score_alsimpch', 'latent_alsimpch', 'tsquared_alsimpch', 'explained_alsimpch', 'mu_alsimpch');
            
            %plot explained variance of each component
            figure('Name','Explained Variance Pairwise ImpChannels')
            hold on
            pareto(explained_pairimpch)
            yline(70)
            yline(80)
            xlabel('Principal Component')
            ylabel('Variance Explained (%)')
            title('Variance of each component Pairwise ImpCh')
            savefig([savepath date 'ExpVarPairImpCh']);
            
            figure('Name','Explained Variance ALS ImpChannels')
            hold on
            pareto(explained_alsimpch)
            yline(70)
            yline(80)
            xlabel('Principal Component')
            ylabel('Variance Explained (%)')
            title('Variance of each component ALS ImpCh')
            savefig([savepath date 'ExpVarALSImpCh']);
            
            %Pair: 70% variance expliquée = 4 composantes, 80% = 6 composantes
            %ALS: 70% variance expliquée = 4 composantes, 80% = 6 composantes
            
            % estimcrdatach_pair = score_pair*coeff_pair';
            % estimcrdatach_als = score_als*coeff_als';
            
            if tablegraphmode == 1
                
                for c = 1:6
                    name = sprintf('C%.0f',c);
                    figure('Name','Score Pairwise ImpCh')
                    plot(score_pairimpch(idG1(1:end-1,1),c),1,'r+',score_pairimpch(idG2(1:end-1,1),c),1,'bo')
                    label = [name ' Principal Component'];
                    xlabel(label);
                    titre = ['Score on ' name ' Pairwise ImpCh'];
                    title(titre)
                    savefig([savepath date 'PairImpCh ' name]);
                    exportgraphics(gcf,[savepath date 'PairImpCh ' name '.png'])
                end
                
                clear label titre c name
                
                for c = 1:6
                    name = sprintf('C%.0f',c);
                    figure('Name','Score ALS ImpCh')
                    plot(score_alsimpch(idG1(1:end-1,1),c),1,'r+',score_alsimpch(idG2(1:end-1,1),c),1,'bo')
                    label = [name ' Principal Component'];
                    xlabel(label);
                    titre = ['Score on ' name ' ALS ImpCh'];
                    title(titre)
                    savefig([savepath date 'ALSImpCh ' name]);
                    exportgraphics(gcf,[savepath date 'ALSImpCh ' name '.png'])
                end
                
                clear label titre c name
                
                %sort tsquared and extract extreme values
                %[st2_pairimpch,index_pairimpch] = sort(tsquared_pairimpch,'descend'); % sort in descending order
                %extreme_pairimpch = index_pairimpch(1);
                %labelch(:,extreme_pairimpch)
                
                %[st2_alsimpch,index_alsimpch] = sort(tsquared_alsimpch,'descend'); % sort in descending order
                %extreme_alsimpch = index_alsimpch(1);
                %labelch(:,extreme_alsimpch)
            end
            toc
%%%%%%%%%%%%%%%%%%ANCOVA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic
        disp('Computing ANCOVA on Imputed Channels Data, Unrotated Components')
        p = 0.07;

        %Pairwise Not Rotated%%%%%
        x = 1;
        p_pcapairimpch = zeros (2,6);
        res_pcapairimpch = [];
        for c = 1:6
            [p_pcapairimpch(:,x),res] = anovan(score_pairimpch(:,c),{gr(1:end-1,1),ses(1:end-1,1)},'Continuous',2,'varnames',{'Group','SES'});
            close hidden;
            res_pcapairimpch = [res_pcapairimpch, res];
            x = x + 1;
        end

        n_sig = sum(p_pcapairimpch(1,:) <= p);
        fprintf('%d tests are significant p<=%.2f for ANCOVA on Pairwise non rotated Imputed Channels components\n',n_sig, p)

        clear r x res n_sig

        %ALS Not Rotated%%%%%
        x = 1;
        p_pcaalsimpch = zeros (2,6);
        res_pcaalsimpch = [];
        for c = 1:6
            [p_pcaalsimpch(:,x),res] = anovan(score_alsimpch(:,c),{gr(1:end-1,1),ses(1:end-1,1)},'Continuous',2,'varnames',{'Group','SES'});
            close hidden;
            res_pcaalsimpch = [res_pcaalsimpch, res];
            x = x + 1;
        end

        n_sig = sum(p_pcaalsimpch(1,:) <= p);
        fprintf('%d tests are significant p<=%.2f for ANCOVA on ALS non rotated Imputed Channels components\n',n_sig, p)

        clear r x res n_sig
        
        end
        toc 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%Rotation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if rotationmode == 1
            tic
            disp('Computing PCA on Imputed Channels Data, Rotated Components')
            
            [coeff_pairimpchrot, matrix_pairimpchrot] = rotatefactors(coeff_pairimpch, 'Method', 'promax','Maxit',1000);
            [coeff_alsimpchrot, matrix_alsimpchrot] = rotatefactors(coeff_alsimpch, 'Method', 'promax','Maxit',1000);
            %[coeff_mdiimpchrot, matrix_mdiimpchrot] = rotatefactors(MDItoolbox_results.Loadings, 'Method', 'promax','Maxit',1000);            
            
            score_pairimpchrot = crdataimpch/coeff_pairimpchrot';
            score_alsimpchrot = crdataimpch/coeff_alsimpchrot';
            %score_mdiimpchrot = crdataimpch/coeff_mdiimpchrot';
            
            latent_pairimpchrot = sort(nanvar(score_pairimpchrot),'descend')';
            latent_alsimpchrot = sort(nanvar(score_alsimpchrot),'descend')';
            %latent_mdiimpchrot = sort(nanvar(score_mdiimpchrot),'descend')';
            
            explained_pairimpchrot = (latent_pairimpchrot*100/sum(latent_pairimpchrot));
            explained_alsimpchrot = (latent_alsimpchrot*100/sum(latent_alsimpchrot));
            %explained_mdiimpchrot = (latent_mdiimpchrot*100/sum(latent_mdiimpchrot));
            
            save([savepath date '_resultsPCA_PairImpChRot.mat'],'coeff_pairimpchrot', 'score_pairimpchrot', 'latent_pairimpchrot', 'explained_pairimpchrot');
            save([savepath date '_resultsPCA_AlsImpChRot.mat'],'coeff_alsimpchrot', 'score_alsimpchrot', 'latent_alsimpchrot', 'explained_alsimpchrot');
            
            
            if tablegraphmode == 1
                
                for c = 1:6
                    name = sprintf('C%.0f',c);
                    figure('Name','Score Pairwise ImpChRot')
                    plot(score_pairimpchrot(idG1(1:end-1,1),c),1,'r+',score_pairimpchrot(idG2(1:end-1,1),c),1,'bo')
                    label = [name ' Principal Component'];
                    xlabel(label);
                    titre = ['Score on ' name ' Pairwise ImpChRot'];
                    title(titre)
                    savefig([savepath date 'PairImpChRot ' name]);
                    exportgraphics(gcf,[savepath date 'PairImpChRot ' name '.png'])
                end
                
                clear label titre c name
                
                for c = 1:6
                    name = sprintf('C%.0f',c);
                    figure('Name','Score ALS ImpChRot')
                    plot(score_alsimpchrot(idG1(1:end-1,1),c),1,'r+',score_alsimpchrot(idG2(1:end-1,1),c),1,'bo')
                    label = [name ' Principal Component'];
                    xlabel(label);
                    titre = ['Score on ' name ' ALS ImpChRot'];
                    title(titre)
                    savefig([savepath date 'ALSImpChRot ' name]);
                    exportgraphics(gcf,[savepath date 'ALSImpChRot ' name '.png'])
                end
                
                clear label titre c name
            end
            toc


%%%%%%%%%%%%%%%%%%ANCOVA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic
        disp('Computing ANCOVA on Imputed Channels Data, Rotated Components')
        p = 0.07;

        %Pairwise RotCh%%%%%
        x = 1;
        p_pcapairimpchrot = zeros (2,6);
        res_pcapairimpchrot = [];
        for c = 1:6
            [p_pcapairimpchrot(:,x),res] = anovan(score_pairimpchrot(:,c),{gr(1:end-1,1),ses(1:end-1,1)},'Continuous',2,'varnames',{'Group','SES'});
            close hidden;
            res_pcapairimpchrot = [res_pcapairimpchrot, res];
            x = x + 1;
        end

        n_sig = sum(p_pcapairimpchrot(1,:) <= p);
        fprintf('%d tests are significant p<=%.2f for ANCOVA on Pairwise Rotated Imputed Channels components\n',n_sig, p)

        clear r x res n_sig

        %ALS RotCh%%%%%
        x = 1;
        p_pcaalsimpchrot = zeros (2,6);
        res_pcaalsimpchrot = [];
        for c = 1:6
            [p_pcaalsimpchrot(:,x),res] = anovan(score_alsimpchrot(:,c),{gr(1:end-1,1),ses(1:end-1,1)},'Continuous',2,'varnames',{'Group','SES'});
            close hidden;
            res_pcaalsimpchrot = [res_pcaalsimpchrot, res];
            x = x + 1;
        end

        n_sig = sum(p_pcaalsimpchrot(1,:) <= p);
        fprintf('%d tests are significant p<=%.2f for ANCOVA on ALS Rotated Imputed Channels components\n',n_sig, p)

        clear r x res n_sig
        
        end
    end
end

toc

%%
%%%%%%%%%%%%%%%%%%%%%%RROI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if roimode == 1
    dataRroi = table2array(tblmeanRroi(1:end-1, 6:end));
    crdataRroi = table2array(tblmeanRroi(1:end-1, 6:end)) - repmat(nanmean(table2array(tblmeanRroi(1:end-1, 6:end))),size(tblmeanRroi(1:end-1, 6:end),1),1);
    if initialmode == 1
        tic
        disp('Computing PCA on Initial Rroi Data, Unrotated components')
        R = 2;
        
        [coeff_pairRroi,score_pairRroi,latent_pairRroi,tsquared_pairRroi,explained_pairRroi,mu_pairRroi] = pca(dataRroi,'Rows','pairwise');
        opt = statset('pca'); opt.MaxIter = 10000;
        [coeff_alsRroi,score_alsRroi,latent_alsRroi,tsquared_alsRroi,explained_alsRroi,mu_alsRroi] = pca(dataRroi,'Algorithm','als','Options',opt);

        save([savepath date '_resultsPCA_PairRroi.mat'],'coeff_pairRroi', 'score_pairRroi', 'latent_pairRroi', 'tsquared_pairRroi', 'explained_pairRroi', 'mu_pairRroi');
        save([savepath date '_resultsPCA_AlsRroi.mat'],'coeff_alsRroi', 'score_alsRroi', 'latent_alsRroi', 'tsquared_alsRroi', 'explained_alsRroi', 'mu_alsRroi');
        
        
        if norotationmode == 1
        %plot explained variance of each component
        figure('Name','Explained Variance Pairwise Rroi')
        hold on
        pareto(explained_pairRroi)
        yline(70)
        yline(80)
        xlabel('Principal Component')
        ylabel('Variance Explained (%)')
        title('Variance of each component Pairwise Rroi')
        savefig([savepath date 'ExpVarPairRroi']);

        figure('Name','Explained Variance ALS Rroi')
        hold on
        pareto(explained_alsRroi)
        yline(70)
        yline(80)
        xlabel('Principal Component')
        ylabel('Variance Explained (%)')
        title('Variance of each component ALS Rroi')
        savefig([savepath date 'ExpVarALSRroi']);

        %Pair: 70% variance expliquée = 6 composantes, 80% = 10 composantes
        %ALS: 70% variance expliquée = 4 composantes, 80% = 6 composantes

        if tablegraphmode == 1
            for c = 1:10
                name = sprintf('C%.0f',c);
                figure('Name','Score PairRroi')
                plot(score_pairRroi(idG1(1:end-1,1),c),1,'r+',score_pairRroi(idG2(1:end-1,1),c),1,'bo')
                label = [name ' Principal Component'];
                xlabel(label);
                titre = ['Score on ' name ' PairRroi'];
                title(titre)
                savefig([savepath date 'PairRroi ' name]);
                exportgraphics(gcf,[savepath date 'PairRroi ' name '.png'])
            end

            for c = 1:6
                name = sprintf('C%.0f',c);
                figure('Name','Score ALSRroi')
                plot(score_alsRroi(idG1(1:end-1,1),c),1,'r+',score_alsRroi(idG2(1:end-1,1),c),1,'bo')
                label = [name ' Principal Component'];
                xlabel(label);
                titre = ['Score on ' name ' ALSRroi'];
                title(titre)
                savefig([savepath date 'ALSRroi ' name]);
                exportgraphics(gcf,[savepath date 'ALSRroi ' name '.png'])
            end
        end
        toc
        % estimcrdata_pairRroi = score_pairRroi*coeff_pairRroi';
        % estimcrdata_alsRroi = score_alsRroi*coeff_alsRroi';

%%%%%%%%%%%%%%%%%%ANCOVA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         tic
         disp('Computing ANCOVA on Initial Rroi Data Components, Unrotated Components')
         p = 0.07;

         %%% Pairwise non rotated%%%%%%%%%%%%
         % Trop de NAN dans les scores, seulement une vingaine de participant sans aucun nan%
         x = 1;
         p_pcapairRroi = zeros (2,10);
         res_pcapairRroi = [];
         for c = 1:10
             [p_pcapairRroi(:,x),res] = anovan(score_pairRroi(:,c),{gr(1:end-1,1),ses(1:end-1,1)},'Continuous',2,'varnames',{'Group','SES'});
             close hidden;
             res_pcapairRroi = [res_pcapairRroi, res];
             x = x + 1;
         end

         n_sig = sum(p_pcapairRroi(1,:) <= p);
         fprintf('%d tests are significant p<=%.2f for ANCOVA on Pairwise NonRotated Initial Rroi components\n',n_sig, p)

         clear r x res n_sig

         %%% ALS non rotated%%%%%%%%%%%%
         x = 1;
         p_pcaalsRroi = zeros (2,6);
         res_pcaalsRroi = [];
         for c = 1:6
             [p_pcaalsRroi(:,x),res] = anovan(score_alsRroi(:,c),{gr(1:end-1,1),ses(1:end-1,1)},'Continuous',2,'varnames',{'Group','SES'});
             close hidden;
             res_pcaalsRroi = [res_pcaalsRroi, res];
             x = x + 1;
         end

         n_sig = sum(p_pcaalsRroi(1,:) <= p);
         fprintf('%d tests are significant p<=%.2f for ANCOVA on ALS NonRotated Initial Rroi components\n',n_sig, p)

         clear r x res n_sig
        end
        toc
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%Rotation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic
        if rotationmode == 1
        disp('Computing PCA on Initial Rroi Data, Rotated components')
        
        [coeff_pairRroirot, matrix_pairRroirot] = rotatefactors(coeff_pairRroi, 'Method', 'promax','Maxit',1000);
        [coeff_alsRroirot, matrix_alsRroirot] = rotatefactors(coeff_alsRroi, 'Method', 'promax','Maxit',1000);

        score_pairRroirot = crdataRroi/coeff_pairRroirot';
        score_alsRroirot = crdataRroi/coeff_alsRroirot';

        latent_pairRroirot = sort(nanvar(score_pairRroirot),'descend')';
        latent_alsRroirot = sort(nanvar(score_alsRroirot),'descend')';

        explained_pairRroirot = (latent_pairRroirot*100/sum(latent_pairRroirot));
        explained_alsRroirot = (latent_alsRroirot*100/sum(latent_alsRroirot));

        %Problème avec ALS, comme je recalcule avec les données originales = NAN%%

        save([savepath date '_resultsPCA_PairRroiRot.mat'],'coeff_pairRroirot', 'score_pairRroirot', 'latent_pairRroirot', 'explained_pairRroirot');
        save([savepath date '_resultsPCA_AlsRroiRot.mat'],'coeff_alsRroirot', 'score_alsRroirot', 'latent_alsRroirot', 'explained_alsRroirot');

        if tablegraphmode == 1
            for c = 1:10
                name = sprintf('C%.0f',c);
                figure('Name','Score Pairwise RroiRot')
                plot(score_pairRroirot(idG1(1:end-1,1),c),1,'r+',score_pairRroirot(idG2(1:end-1,1),c),1,'bo')
                label = [name ' Principal Component'];
                xlabel(label);
                titre = ['Score on ' name ' Pairwise RroiRot'];
                title(titre)
                savefig([savepath date 'PairRroiRot ' name]);
                exportgraphics(gcf,[savepath date 'PairRroiRot ' name '.png'])
            end

            clear label titre c name

            for c = 1:6
                name = sprintf('C%.0f',c);
                figure('Name','Score ALS RroiRot')
                plot(score_alsRroirot(idG1(1:end-1,1),c),1,'r+',score_alsRroirot(idG2(1:end-1,1),c),1,'bo')
                label = [name ' Principal Component'];
                xlabel(label);
                titre = ['Score on ' name ' ALS RroiRot'];
                title(titre)
                savefig([savepath date 'ALSRroiRot ' name]);
                exportgraphics(gcf,[savepath date 'ALSRroiRot ' name '.png'])
            end

            clear label titre c name
        end
        toc

%%%%%%%%%%%%%%%%%%ANCOVA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         tic
         disp('Computing ANCOVA on Initial Rroi Data, Rotated Components')
         p = 0.07;

         %%Pas la peine de faire des ANOVA sur les données rotatées, car trop de NAN avec les deux méthodes%%%
         %%%Pairwise rotated
         x = 1;
         p_pcapairRroirot = zeros (2,10);
         res_pcapairRroirot = [];
         for c = 1:10
             [p_pcapairRroirot(:,x),res] = anovan(score_pairRroirot(:,c),{gr(1:end-1,1),ses(1:end-1,1)},'Continuous',2,'varnames',{'Group','SES'});
             close hidden;
             res_pcapairRroirot = [res_pcapairRroirot, res];
             x = x + 1;
         end

         n_sig = sum(p_pcapairRroirot(1,:) <= p);
         fprintf('%d tests are significant p<=%.2f for ANCOVA on Pairwise Rotated Initial Rroi components\n',n_sig, p)

         clear r x res n_sig

         %%% ALS rotated%%%%%%%%%%%%
         x = 1;
         p_pcaalsRroirot = zeros (2,6);
         res_pcaalsRroirot = [];
         for c = 1:6
             [p_pcaalsRroirot(:,x),res] = anovan(score_alsRroirot(:,c),{gr(1:end-1,1),ses(1:end-1,1)},'Continuous',2,'varnames',{'Group','SES'});
             close hidden;
             res_pcaalsRroirot = [res_pcaalsRroirot, res];
             x = x + 1;
         end

         n_sig = sum(p_pcaalsRroirot(1,:) <= p);
         fprintf('%d tests are significant p<=%.2f for ANCOVA on ALS Rotated Initial Rroi components\n',n_sig, p)

         clear r x res n_sig
        end
    end
    toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%PCA sur données imputées%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if imputemode == 1
         tic
         disp('Computing PCA on Imputed Rroi Data, Unrotated Components')

     %%%%Imputations%%%%%%%%

         %load ('MDItoolbox_resultsRroi.mat')
         %varargout = NumberComponents(varargin);
         %expvarRroi = NumberComponents(dataRroi);

         %%On ne voit que les 10 premières composantes, devrait être de 8/12 comme établi plus haut%%%%
         %%%% 70% = 8 composantes, 80% = 14 composantes VS 7 et 10 pour pairwise Rroi%
         
         %[impdata,impmean,impcov,recdata] = pcambda(data,Nchains,Lchains,Ncomp);
         %[impdataRroi,impmeanRroi,impcovRroi,recdataRroi] = pcambda(dataRroi,10,100,10); %matrice not positive definite%

         %[impdata,impmean,impcov,It,diff,recdata]=pcambia(data,Ncomp,maxiter,tolerance)
         [dataimpRroi,meanimpRroi,covimpRroi,ItRroi,diffRroi,recdataRroi] = pcambia(dataRroi,8,10000,10E-10);

         crdataimpRroi = dataimpRroi - repmat(nanmean(dataimpRroi),size(dataimpRroi,1),1);

     %%%%%PCA%%%%%%%
         [coeff_pairimpRroi,score_pairimpRroi,latent_pairimpRroi,tsquared_pairimpRroi,explained_pairimpRroi,mu_pairimpRroi] = pca(dataimpRroi,'Rows','pairwise');
         opt = statset('pca'); opt.MaxIter = 10000;
         [coeff_alsimpRroi,score_alsimpRroi,latent_alsimpRroi,tsquared_alsimpRroi,explained_alsimpRroi,mu_alsimpRroi] = pca(dataimpRroi,'Algorithm','als','Options',opt);

         save([savepath date '_resultsPCA_PairImpRroi.mat'],'coeff_pairimpRroi', 'score_pairimpRroi', 'latent_pairimpRroi', 'tsquared_pairimpRroi', 'explained_pairimpRroi', 'mu_pairimpRroi');
         save([savepath date '_resultsPCA_AlsImpRroi.mat'],'coeff_alsimpRroi', 'score_alsimpRroi', 'latent_alsimpRroi', 'tsquared_alsimpRroi', 'explained_alsimpRroi', 'mu_alsimpRroi');
        
         if norotationmode == 1

             %plot explained variance of each component
             figure('Name','Explained Variance Pairwise ImpRroi')
             hold on
             pareto(explained_pairimpRroi)
             yline(70)
             yline(80)
             xlabel('Principal Component')
             ylabel('Variance Explained (%)')
             title('Variance of each component Pairwise ImpRroi')
             savefig([savepath date 'ExpVarPairImpRroi']);

             figure('Name','Explained Variance ALS ImpRroi')
             hold on
             pareto(explained_alsimpRroi)
             yline(70)
             yline(80)
             xlabel('Principal Component')
             ylabel('Variance Explained (%)')
             title('Variance of each component ALS ImpRroi')
             savefig([savepath date 'ExpVarALSImpRroi']);

             %Pair: 70% variance expliquée = 4 composantes, 80% = 6 composantes
             %ALS: 70% variance expliquée = 4 composantes, 80% = 6 composantes %% À VALIDER

             % estimcrdataRroi_pair = score_pair*coeff_pair';
             % estimcrdataRroi_als = score_als*coeff_als';

             if tablegraphmode == 1
                 for c = 1:6
                     name = sprintf('C%.0f',c);
                     figure('Name','Score Pairwise ImpRroi')
                     plot(score_pairimpRroi(idG1(1:end-1,1),c),1,'r+',score_pairimpRroi(idG2(1:end-1,1),c),1,'bo')
                     label = [name ' Principal Component'];
                     xlabel(label);
                     titre = ['Score on ' name ' Pairwise ImpRroi'];
                     title(titre)
                     savefig([savepath date 'PairImpRroi ' name]);
                     exportgraphics(gcf,[savepath date 'PairImpRroi ' name '.png'])
                 end

                 clear label titre c name

                 for c = 1:6
                     name = sprintf('C%.0f',c);
                     figure('Name','Score ALS ImpRroi')
                     plot(score_alsimpRroi(idG1(1:end-1,1),c),1,'r+',score_alsimpRroi(idG2(1:end-1,1),c),1,'bo')
                     label = [name ' Principal Component'];
                     xlabel(label);
                     titre = ['Score on ' name ' ALS ImpRroi'];
                     title(titre)
                     savefig([savepath date 'ALSImpRroi ' name]);
                     exportgraphics(gcf,[savepath date 'ALSImpRroi ' name '.png'])
                 end

                 clear label titre c name

                 %sort tsquared and extract extreme values
                 %[st2_pairimpRroi,index_pairimpRroi] = sort(tsquared_pairimpRroi,'descend'); % sort in descending order
                 %extreme_pairimpRroi = index_pairimpRroi(1);
                 %labelroi{1,R}(:,extreme_pairimpRroi)

                 %[st2_alsimpRroi,index_alsimpRroi] = sort(tsquared_alsimpRroi,'descend'); % sort in descending order
                 %extreme_alsimpRroi = index_alsimpRroi(1);
                 %labelroi{1,R}(:,extreme_alsimpRroi)
             end
            toc
             
%%%%%%%%%%%%%%%%%%ANCOVA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         disp('Computing ANCOVA on Imputed Rroi Data, Unrotated Components')
         tic
         p = 0.07;

         %Pairwise Not Rotated%%%%%
         x = 1;
         p_pcapairimpRroi = zeros (2,6);
         res_pcapairimpRroi = [];
         for c = 1:6
             [p_pcapairimpRroi(:,x),res] = anovan(score_pairimpRroi(:,c),{gr(1:end-1,1),ses(1:end-1,1)},'Continuous',2,'varnames',{'Group','SES'});
             close hidden;
             res_pcapairimpRroi = [res_pcapairimpRroi, res];
             x = x + 1;
         end

         n_sig = sum(p_pcapairimpRroi(1,:) <= p);
         fprintf('%d tests are significant p<=%.2f for ANCOVA on Pairwise NonRotated Imputed Rroi components\n',n_sig, p)

         clear r x res n_sig

         %ALS Not Rotated%%%%%
         x = 1;
         p_pcaalsimpRroi = zeros (2,6);
         res_pcaalsimpRroi = [];
         for c = 1:6
             [p_pcaalsimpRroi(:,x),res] = anovan(score_alsimpRroi(:,c),{gr(1:end-1,1),ses(1:end-1,1)},'Continuous',2,'varnames',{'Group','SES'});
             close hidden;
             res_pcaalsimpRroi = [res_pcaalsimpRroi, res];
             x = x + 1;
         end

         n_sig = sum(p_pcaalsimpRroi(1,:) <= p);
         fprintf('%d tests are significant p<=%.2f for ANCOVA on ALS NonRotated Imputed Rroi components\n',n_sig, p)

         clear r x res n_sig
         end
         toc
             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%Rotation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         tic
         disp('Computing PCA on Imputed Rroi Data, Rotated Components')
         if rotationmode == 1
             [coeff_pairimpRroirot, matrix_pairimpRroirot] = rotatefactors(coeff_pairimpRroi, 'Method', 'promax','Maxit',1000);
             [coeff_alsimpRroirot, matrix_alsimpRroirot] = rotatefactors(coeff_alsimpRroi, 'Method', 'promax','Maxit',1000);

             score_pairimpRroirot = crdataimpRroi/coeff_pairimpRroirot';
             score_alsimpRroirot = crdataimpRroi/coeff_alsimpRroirot';

             latent_pairimpRroirot = sort(nanvar(score_pairimpRroirot),'descend')';
             latent_alsimpRroirot = sort(nanvar(score_alsimpRroirot),'descend')';

             explained_pairimpRroirot = (latent_pairimpRroirot*100/sum(latent_pairimpRroirot));
             explained_alsimpRroirot = (latent_alsimpRroirot*100/sum(latent_alsimpRroirot));

             save([savepath date '_resultsPCA_PairImpRroiRot.mat'],'coeff_pairimpRroirot', 'score_pairimpRroirot', 'latent_pairimpRroirot', 'explained_pairimpRroirot');
             save([savepath date '_resultsPCA_AlsImpRroiRot.mat'],'coeff_alsimpRroirot', 'score_alsimpRroirot', 'latent_alsimpRroirot', 'explained_alsimpRroirot');

             if tablegraphmode == 1

                 for c = 1:6
                     name = sprintf('C%.0f',c);
                     figure('Name','Score Pairwise ImpRroiRot')
                     plot(score_pairimpRroirot(idG1(1:end-1,1),c),1,'r+',score_pairimpRroirot(idG2(1:end-1,1),c),1,'bo')
                     label = [name ' Principal Component'];
                     xlabel(label);
                     titre = ['Score on ' name ' Pairwise ImpRroiRot'];
                     title(titre)
                     savefig([savepath date 'PairImpRroiRot ' name]);
                     exportgraphics(gcf,[savepath date 'PairImpRroiRot ' name '.png'])
                 end

                 clear label titre c name

                 for c = 1:6
                     name = sprintf('C%.0f',c);
                     figure('Name','Score ALS ImpRroiRot')
                     plot(score_alsimpRroirot(idG1(1:end-1,1),c),1,'r+',score_alsimpRroirot(idG2(1:end-1,1),c),1,'bo')
                     label = [name ' Principal Component'];
                     xlabel(label);
                     titre = ['Score on ' name ' ALS ImpRroiRot'];
                     title(titre)
                     savefig([savepath date 'ALSImpRroiRot ' name]);
                     exportgraphics(gcf,[savepath date 'ALSImpRroiRot ' name '.png'])
                 end

                 clear label titre c name
             end
         toc


%%%%%%%%%%%%%%%%%%ANCOVA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         disp('Computing ANCOVA on Imputed Rroi Data, Rotated Components')
         tic
         p = 0.07;

         %Pairwise RroiRot%%%%%
         x = 1;
         p_pcapairimpRroirot = zeros (2,6);
         res_pcapairimpRroirot = [];
         for c = 1:6
             [p_pcapairimpRroirot(:,x),res] = anovan(score_pairimpRroirot(:,c),{gr(1:end-1,1),ses(1:end-1,1)},'Continuous',2,'varnames',{'Group','SES'});
             close hidden;
             res_pcapairimpRroirot = [res_pcapairimpRroirot, res];
             x = x + 1;
         end

         n_sig = sum(p_pcapairimpRroirot(1,:) <= p);
         fprintf('%d tests are significant p<=%.2f for ANCOVA on Pairwise Rotated Imputed Rroi components\n',n_sig, p)

         clear r x res n_sig

         %ALS RroiRot%%%%%
         x = 1;
         p_pcaalsimpRroirot = zeros (2,6);
         res_pcaalsimpRroirot = [];
         for c = 1:6
             [p_pcaalsimpRroirot(:,x),res] = anovan(score_alsimpRroirot(:,c),{gr(1:end-1,1),ses(1:end-1,1)},'Continuous',2,'varnames',{'Group','SES'});
             close hidden;
             res_pcaalsimpRroirot = [res_pcaalsimpRroirot, res];
             x = x + 1;
         end

         n_sig = sum(p_pcaalsimpRroirot(1,:) <= p);
         fprintf('%d tests are significant p<=%.2f for ANCOVA on ALS Rotated Imputed Rroi components\n',n_sig, p)

         clear r x res n_sig
        end
    end
end
toc