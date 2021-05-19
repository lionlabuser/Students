%%%%%%%%%%%%%%ANCOVA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
disp('Computing ANCOVA')

datapath = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Stats\CORR0,01_0,08\Physio\';
load ([datapath 'workspace.mat'])
load ([datapath 'workspacemat.mat'])

savepathinitial='C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Stats\CORR0,01_0,08\Physio\ANCOVA\';
if ~isfolder(savepathinitial)
    mkdir(savepathinitial)
end
channelmode = 1;
importedROImode = 0;
calculatedROImode = 1;

tablegraphmode = 1;
%fileorderconnectogram  = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Mixte.txt';
fdrmode = 1;
nointeractionmode = 0;
interactionmode = 1;
p = 0.07;

%%
%%ANCOVA SANS INTERACTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nointeractionmode == 1
disp('Computing ANCOVA without interaction')
savepath = [savepathinitial 'NoInteraction\'];
if ~isfolder(savepath)
    mkdir(savepath)
end

%%CHANNELS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if channelmode ==1
   
       %%ANOVA%%%%%%%%%%%%%%%%%%%
        x = 1;
        resch = [];
        for r = 1:(width(tblch)- 5)
            [pch(:,x),res] = anovan(table2array(tblch(:,r+5)),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
            close hidden;
            resch = [resch, res];
            x = x + 1;
        end
        % delete(findall(0));
        
        labelfactors = resch(2:end-2,1);
        for f = 1:numel(labelfactors)
            n_sig = sum(pch(f,:) <= p);
            fprintf('%d %s effects are significant p<=%.2f without correction for channels\n',n_sig,labelfactors{f,1}, p)
        end
        
        clear r x res n_sig

        %%fdr%%%%%%%%%%%%
        if fdrmode == 1
           for f = 1:numel(labelfactors)
             % FDR according to Storey 2002%%%
             [fdr_ch,q_ch] = mafdr(pch(f,:)');
             n_sig = sum(q_ch <= p);
             fprintf('%d %s effects are significant p<=%.2f using FDR correction (Storey2002) for channels\n',n_sig,labelfactors{f,1}, p)

            % FDR according to Benjamini & Hochberg(1995) de Groppe%%%%
            [h_ch1, crit_p_ch, adj_ci_cvrg_ch, adj_p_ch] = fdr_bh(pch(f,:), 0.05, 'pdep', 'no' );
            n_sig = sum(adj_p_ch(1,:) <= p);
            fprintf('%d %s effects are significant p<=%.2f using FDR correction (Benjamini1995) for channels\n',n_sig,labelfactors{f,1}, p)

            % % FDR according to Benjamini & Hochberg(1995) de Gerber%%%%
            % [ind, thres] = FDR(pch(f,:));

            %Bonferroni-Holm correction by Groppe%%% Corrected_p are greater than 1...
            %doesn't make sense
            [corrected_p_ch, h_ch2] = bonf_holm(pch(f,:));
            n_sig = sum(corrected_p_ch <= p);
            fprintf('%d %s effects are significant p<=%.2f using Bonferroni-Holm correction for channels\n',n_sig,labelfactors{f,1}, p)

            %Multcmp function to apply Holm Step Down Procedure%%
            [padjBH_ch,alphaBH_ch] = multicmp (pch(f,:)','down',0.05);
            n_sig = sum(padjBH_ch <= p);
            fprintf('%d %s effects are significant p<=%.2f using Holm Step Down Procedure for channels\n',n_sig,labelfactors{f,1}, p)

            %Multcmp function to apply Hochberg's step up procedure%%
            [padjH_ch,alphaH_ch] = multicmp (pch(f,:)','up',0.05);
            n_sig = sum(padjH_ch <= p);
            fprintf('%d %s effects are significant p<=%.2f using Hochberg Step Up Procedure for channels\n',n_sig,labelfactors{f,1}, p)

            %Multicmp function to apply FDR correction according to Benjamini & Hochberg%%
            [padjFDR_ch,alphaFDR_ch] = multicmp (pch(f,:)','fdr',0.05);
            n_sig = sum(padjFDR_ch <= p);
            fprintf('%d %s effects are significant p<=%.2f using FDR correction for channels\n',n_sig,labelfactors{f,1}, p)
           end
        end
        
        %%tableaux des résultats%%%%%%%%%%%%%%%%
        if tablegraphmode == 1
            tblresch = array2table(resch);
            for r = 1:(width(tblch)- 5)
                tblresch = mergevars(tblresch,(r:r+6));
            end
            tblresch.Properties.VariableNames = labelch;
            tblresch.Properties.RowNames = {'Source','Group','SES','Error','Total'};

            tblpch = array2table(pch);
            tblpch.Properties.VariableNames = labelch;
            tblpch.Properties.RowNames = {'Group','SES'};

            sig = tblpch.Variables <= p;
            tblpsigch = tblpch(:,any(sig));
            tblsigch = tblresch(:,tblpsigch.Properties.VariableNames);

            grsig = tblpch{{'Group'},:} <= p;
            tblgrpsigch = tblpch(1,grsig);
            tblgrsigch = tblresch(:, tblgrpsigch.Properties.VariableNames);

            disp(tblgrpsigch)

            %writetable(tblsigch,[savepath date '_sigresultsch.xls']); trop gros
            writetable(tblpsigch,[savepath date '_psigresultsch.xls'],'WriteRowNames',true);
            %writetable(tblgrsigch,[savepath date '_grsigresultsch.xls']);
            save([savepath date '_resultsch.mat'],'tblresch','tblpch','tblpsigch','tblsigch','tblgrpsigch','tblgrsigch');

            clear r

            %%%% graphique%%%%%%%%%%%%%%%%%%%%
            if find(grsig)
                figure
                A = meanG1ch(:,grsig);
                B = meanG2ch(:,grsig);
                C = [A' B'];
                X = categorical(tblgrpsigch.Properties.VariableNames);
                p1 = bar(X,C);
                p1(1).FaceColor = 'r';
                p1(2).FaceColor = 'b';
                ylabel('Pearson correlation');
                savefig([savepath date '_SigCh']);
                exportgraphics(gcf,[savepath 'SigCh.png'])

            else
            end

            clear A B C X p1 sig grsig

            %%%% Connectogramme%%%%%%%%%%%%
%             %MATmeanG2G1 = MATmeanG2-MATmeanG1; %Calculer matrice G2-G1
%             MATmeanG1G2 = MATmeanG1 - MATmeanG2;

            for c=1:size(MATall,1)
                for cc =(c + 1):size(MATall,1)
                    x = [num2str(c) '-' num2str(cc)];
                    tf = strcmp(x, labelch);
                    idx = find(tf);
                    MATpch(c,cc)= pch(1,idx); %Mettre les  p sous forme de matrice
                end
            end

            MATpch(length(MATpch),length(MATpch)) = 0; 
            MATpch = MATpch + triu(MATpch,1)';
            tf = MATpch <=p & MATpch > 0;
    %         idx = find(tf);
            %MATsigG2G1ch = MATmeanG2G1;
            MATsigG1G2ch = MATmeanG1G2;
            %MATsigG2G1ch(~tf) = 0;
            MATsigG1G2ch(~tf) = 0;

            clear c cc x tf idx

            id = 1;
            MATneg = MATsigG1G2ch < 0;
            MAT = MATsigG1G2ch.*(MATneg);%DATA{id}.MAT %loader matrice
            if find(MAT)
                List = strvcat(DATA{id}.ZoneList); %liste des paires SD
                ML = DATA{id}.zone.ml; %Loader S/D/ROI/Gr
                plotLst = DATA{id}.zone.plotLst;
                label =  DATA{id}.zone.label;
                %fileorderconnectogram  = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Mixte.txt';
                plotconnectogram(fileorderconnectogram{1,1},MAT,List,label,plotLst,ML)
                savefig([savepath date '_NegConnectSigG1G2Ch']);
            else
            end

            MATpos = MATsigG1G2ch > 0;
            MAT = MATsigG1G2ch.*(MATpos);%DATA{id}.MAT %loader matrice
            if find(MAT)
                List = strvcat(DATA{id}.ZoneList); %liste des paires SD
                ML = DATA{id}.zone.ml; %Loader S/D/ROI/Gr
                plotLst = DATA{id}.zone.plotLst;
                label =  DATA{id}.zone.label;
                %fileorderconnectogram  = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Mixte.txt';
                plotconnectogram(fileorderconnectogram{1,1},MAT,List,label,plotLst,ML)
                savefig([savepath date '_PosConnectSigG1G2Ch']);
            else
            end

        end
    clear r x res sig A B C X p1
   end

      if importedROImode
   
       %%ANOVA%%%%%%%%%%%%%%%%%%%
        x = 1;
        resroi = [];
        for r = 1:(width(tblroi)- 5)
            [proi(:,x),res] = anovan(table2array(tblroi(:,r+5)),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
            close hidden;
            resroi = [resroi, res];
            x = x + 1;
        end
        % delete(findall(0));
        
        labelfactors = resroi(2:end-2,1);
        for f = 1:numel(labelfactors)
            n_sig = sum(proi(f,:) <= p);
            fprintf('%d %s effects are significant p<=%.2f without correction for ROIs\n',n_sig,labelfactors{f,1}, p)
        end
 
        clear r x res n_sig

        %%fdr%%%%%%%%%%%%
        if fdrmode == 1
            for f = 1:numel(labelfactors)
                 % FDR according to Storey 2002%%%
                 [fdr_roi,q_roi] = mafdr(proi(f,:)');
                 n_sig = sum(q_roi <= p);
                 fprintf('%d %s effects are significant p<=%.2f using FDR correction (Storey2002) for ROIs\n',n_sig,labelfactors{f,1}, p)

                % FDR according to Benjamini & Hochberg(1995) de Groppe%%%%
                [h_roi1, crit_p_roi, adj_ci_cvrg_roi, adj_p_roi] = fdr_bh(proi(f,:), 0.05, 'pdep', 'no' );
                n_sig = sum(adj_p_roi(1,:) <= p);
                fprintf('%d %s effects are significant p<=%.2f using FDR correction (Benjamini1995) for ROIs\n',n_sig,labelfactors{f,1}, p)

                % % FDR according to Benjamini & Hochberg(1995) de Gerber%%%%
                % [ind, thres] = FDR(pch(f,:));

                %Bonferroni-Holm correction by Groppe%%% Corrected_p are greater than 1...
                %doesn't make sense
                [corrected_p_roi, h_roi2] = bonf_holm(proi(f,:));
                n_sig = sum(corrected_p_roi <= p);
                fprintf('%d %s effects are significant p<=%.2f using Bonferroni-Holm correction for ROIs\n',n_sig,labelfactors{f,1}, p)

                %Multcmp function to apply Holm Step Down Procedure%%
                [padjBH_roi,alphaBH_roi] = multicmp (proi(f,:)','down',0.05);
                n_sig = sum(padjBH_roi <= p);
                fprintf('%d %s effects are significant p<=%.2f using Holm Step Down Procedure for ROIs\n',n_sig,labelfactors{f,1}, p)

                %Multcmp function to apply Hochberg's step up procedure%%
                [padjH_roi,alphaH_roi] = multicmp (proi(f,:)','up',0.05);
                n_sig = sum(padjH_roi <= p);
                fprintf('%d %s effects are significant p<=%.2f using Hochberg Step Up Procedure for ROIs\n',n_sig,labelfactors{f,1}, p)

                %Multicmp function to apply FDR correction according to Benjamini & Hochberg%%
                [padjFDR_roi,alphaFDR_roi] = multicmp (proi(f,:)','fdr',0.05);
                n_sig = sum(padjFDR_roi <= p);
                fprintf('%d %s effects are significant p<=%.2f using FDR correction for ROIs\n',n_sig,labelfactors{f,1}, p)
            end
        end
        
        %%tableaux des résultats%%%%%%%%%%%%%%%%
        if tablegraphmode == 1
            tblresroi = array2table(resroi);
            for r = 1:(width(tblroi)- 5)
                tblresroi = mergevars(tblresroi,(r:r+6));
            end
            tblresroi.Properties.VariableNames = labelroi;
            tblresroi.Properties.RowNames = {'Source','Group','SES','Error','Total'};

            tblproi = array2table(proi);
            tblproi.Properties.VariableNames = labelroi;
            tblproi.Properties.RowNames = {'Group','SES'};

            sig = tblproi.Variables <= p;
            tblpsigroi = tblproi(:,any(sig));
            tblsigroi = tblresroi(:,tblpsigroi.Properties.VariableNames);

            grsig = tblproi{{'Group'},:} <= p;
            tblgrpsigroi = tblproi(1,grsig);
            tblgrsigroi = tblresroi(:, tblgrpsigroi.Properties.VariableNames);

            disp(tblgrpsigroi)

            %writetable(tblsigroi,[savepath date '_sigresultsroi.xls']); trop gros
            writetable(tblpsigroi,[savepath date '_psigresultsroi.xls'],'WriteRowNames',true);
            %writetable(tblgrsigroi,[savepath date '_grsigresultsroi.xls']);
            save([savepath date '_resultsroi.mat'],'tblresroi','tblproi','tblpsigroi','tblsigroi','tblgrpsigroi','tblgrsigroi');

            clear r

            %%%% graphique%%%%%%%%%%%%%%%%%%%%
            if find(grsig)
                figure
                A = meanG1roi(:,grsig);
                B = meanG2roi(:,grsig);
                C = [A' B'];
                X = categorical(tblgrpsigroi.Properties.VariableNames);
                p1 = bar(X,C);
                p1(1).FaceColor = 'r';
                p1(2).FaceColor = 'b';
                ylabel('Pearson correlation');
                savefig([savepath date '_SigROI']);
                exportgraphics(gcf,[savepath 'SigROI.png'])

            else
            end

            clear A B C X p1 sig grsig
            
            %%%% Connectogramme%%%%%%%%%%%%
            %meanG2G1roi = meanG2roi-meanG1roi; %Calculer matrice G2-G1
            meanG1G2roi = meanG1roi-meanG2roi; %Calculer matrice G1-G2

            for r = 1:numel(roi(2,:))
                for rr = (r + 1):numel(roi(2,:))
                    x = ['R' num2str(r) '-' 'R' num2str(rr)];
                    tf = strcmp(x, labelroiALL(1,:));
                    idx = find(tf);
                    MATmeanG1G2roi(r,rr) = meanG1G2roi(1,idx); %Mettre les  données sous forme de matrice
                end
            end
            MATmeanG1G2roi(length(MATmeanG1G2roi),length(MATmeanG1G2roi)) = 0; %Ajouter le dernier 0 de la diagonale
            MATmeanG1G2roi = MATmeanG1G2roi + triu(MATmeanG1G2roi,1)'; %Répliquer la moitié inférieure de la matrice

            for r = 1:numel(roi(2,:))
                for rr = (r + 1):numel(roi(2,:))
                    x = ['R' num2str(r) '-' 'R' num2str(rr)];
                    tf = strcmp(x, labelroiALL(1,:));
                    idx = find(tf);
                    MATproi(r,rr)= proi(1,idx); %Mettre les  p sous forme de matrice
                end
            end

            MATproi(length(MATproi),length(MATproi)) = 0; %Ajouter le dernier 0 de la diagonale
            MATproi = MATproi + triu(MATproi,1)'; %Répliquer la moitié inférieure de la matrice
            tf = MATproi <=p & MATproi > 0; %Identifier les p significatifs
            MATsigG1G2roi = MATmeanG1G2roi; %Reprendre la matrice de connectivité du ROI
            MATsigG1G2roi(~tf) = 0; % Retirer les valeurs non sig

            clear r rr x tf idx

            MATneg = MATsigG1G2roi < 0;
            MAT = MATsigG1G2roi.*(MATneg);%DATA{id}.MAT %loader matrice
            if find(MAT)
            plotLst = roi(2,:);
            label =  roi(1,:);
            plotconnectogramroi(fileorderconnectogram,MAT,label,plotLst)
            savefig([savepath date '_NegConnectSigG1G2ROI']);
            else
            end     

            MATpos = MATsigG1G2roi > 0;
            MAT = MATsigG1G2roi.*(MATpos);%DATA{id}.MAT %loader matrice
            if find(MAT)
            plotLst = roi(2,:);
            label =  roi(1,:);
            plotconnectogramroi(fileorderconnectogram,MAT,label,plotLst)
            savefig([savepath date '_PosConnectSigG1G2ROI']);
            else
            end  
            
            clear r x res sig A B C X p1
        end
      end
      
   if calculatedROImode
        %%MROI%%%%%%%%%%%%%%%%%%%%%%%
        %%ANOVA%%%%%%%%%%%%%%%%%%%
        R = 1;
        x = 1;
        resMroi = [];
        for r = 1:(width(tbldataMroi)- 5)
            [pMroi(:,x),res] = anovan(table2array(tbldataMroi(:,r+5)),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
            close hidden;
            resMroi = [resMroi, res];
            x = x + 1;
        end
        % delete(findall(0));

        labelfactors = resMroi(2:end-2,1);
        for f = 1:numel(labelfactors)
            n_sig = sum(pMroi(f,:) <= p);
            fprintf('%d %s effects are significant p<=%.2f without correction for Mroi\n',n_sig,labelfactors{f,1}, p)
        end
 
        %%fdr%%%%%%%%%%%%
        if fdrmode == 1
           for f = 1:numel(labelfactors)
            % % FDR according to Storey 2002%%%
            [fdr_Mroi,q_Mroi] = mafdr(pMroi(f,:)');
            n_sig = sum(q_Mroi <= p);
            fprintf('%d %s effects are significant p <=%.2f using FDR correction for Mroi\n',n_sig,labelfactors{f,1}, p)

            % FDR according to Benjamini & Hochberg(1995) de Groppe%%%%
            [h_Mroi1, crit_p_Mroi, adj_ci_cvrg_Mroi, adj_p_Mroi] = fdr_bh(pMroi(f,:), 0.05, 'pdep', 'no' );
            n_sig = sum(adj_p_Mroi <= p);
            fprintf('%d %s effects are significant p<=%.2f using FDR correction for Mroi\n',n_sig,labelfactors{f,1}, p)

            % % FDR according to Benjamini & Hochberg(1995) de Gerber%%%%
            % [ind, thres] = FDR(pMroi(f,:));

            %Bonferroni-Holm correction by Groppe%%% Corrected_p are greater than 1...
            %doesn't make sense
            [corrected_p_Mroi, h_Mroi2] = bonf_holm(pMroi(f,:));
            n_sig = sum(corrected_p_Mroi <= p);
            fprintf('%d %s effects are significant p<=%.2f using Bonferroni-Holm correction for Mroi\n',n_sig,labelfactors{f,1}, p)

            %Multcmp function to apply Holm Step Down Procedure%%
            [padjBH_Mroi,alphaBH_Mroi] = multicmp (pMroi(f,:)','down',0.05);
            n_sig = sum(padjBH_Mroi <= p);
            fprintf('%d %s effects are significant p<=%.2f using Holm Step Down Procedure for Mroi\n',n_sig,labelfactors{f,1}, p)

            %Multcmp function to apply Hochberg's step up procedure%%
            [padjH_Mroi,alphaH_Mroi] = multicmp (pMroi(f,:)','up',0.05);
            n_sig = sum(padjH_Mroi <= p);
            fprintf('%d %s effects are significant p<=%.2f using Hochberg Step Up Procedure for Mroi\n',n_sig,labelfactors{f,1}, p)

            %Multicmp function to apply FDR correction according to Benjamini & Hochberg%%
            [padjFDR_Mroi,alphaFDR_Mroi] = multicmp (pMroi(f,:)','fdr',0.05);
            n_sig = sum(padjFDR_Mroi <= p);
            fprintf('%d %s effects are significant p<=%.2f using FDR correction for Mroi\n',n_sig,labelfactors{f,1}, p)
           end
        end

        %%tableaux des résultats%%%%%%%%%%%%%%%%
        if tablegraphmode == 1
            tblresMroi = array2table(resMroi);
            for r = 1:(width(tbldataMroi)- 5)
                tblresMroi = mergevars(tblresMroi,(r:r+6));
            end
            tblresMroi.Properties.VariableNames = labelroiALL{1,R};
            tblresMroi.Properties.RowNames = {'Source','Group','SES','Error','Total'};

            tblpMroi = array2table(pMroi);
            tblpMroi.Properties.VariableNames = labelroiALL{1,R};
            tblpMroi.Properties.RowNames = {'Group','SES'};

            sig = tblpMroi.Variables <= p;
            tblpsigMroi = tblpMroi(:,any(sig));
            tblsigMroi = tblresMroi(:,tblpsigMroi.Properties.VariableNames);

            grsig = tblpMroi{{'Group'},:} <= p;
            tblgrpsigMroi = tblpMroi(1,grsig);
            tblgrsigMroi = tblresMroi(:, tblgrpsigMroi.Properties.VariableNames);

            disp(tblgrpsigMroi)
            writetable(tblsigMroi,[savepath date '_sigresultsMroi.xls']);
            writetable(tblgrsigMroi,[savepath date '_grsigresultsMroi.xls']);
            save([savepath date '_resultsMroi.mat'],'tblresMroi','tblpMroi','tblpsigMroi', 'tblsigMroi','tblgrpsigMroi','tblgrsigMroi');

            %%%% graphique%%%%%%%%%%%%%%%%%%%%
            if find(grsig)
                figure
                A = meanroiG1r{1,R}(:,grsig);
                B = meanroiG1r{1,R}(:,grsig);
                C = [A' B'];
                X = categorical(tblgrpsigMroi.Properties.VariableNames);
                p1 = bar(X,C);
                p1(1).FaceColor = 'r';
                p1(2).FaceColor = 'b';
                ylabel('Pearson correlation');
                savefig([savepath date '_SigMroi']);
                exportgraphics(gcf,[savepath 'SigMroi.png'])
            else
            end

            %%%% Connectogramme%%%%%%%%%%%%
%             %meanG2G1Mroi = meanG2Mroi-meanG1Mroi; %Calculer matrice G2-G1
%             meanG1G2Mroi = meanG1Mroi-meanG2Mroi; %Calculer matrice G1-G2
% 
%             for r = 1:numel(roi{1,R}(2,:))
%                 for rr = (r + 1):numel(roi{1,R}(2,:))
%                     x = [lgndroi{2,R} num2str(r) '-' lgndroi{2,R} num2str(rr)];
%                     tf = strcmp(x, labelroiALL{1,R}(1,:));
%                     idx = find(tf);
%                     MATmeanG1G2Mroi(r,rr) = meanG1G2Mroi(1,idx); %Mettre les  p sous forme de matrice
%                 end
%             end
%             MATmeanG1G2Mroi(length(MATmeanG1G2Mroi),length(MATmeanG1G2Mroi)) = 0; %Ajouter le dernier 0 de la diagonale
%             MATmeanG1G2Mroi = MATmeanG1G2Mroi + triu(MATmeanG1G2Mroi,1)'; %Répliquer la moitié inférieure de la matrice

            for r = 1:numel(roi{1,R}(2,:))
                for rr = (r + 1):numel(roi{1,R}(2,:))
                    x = [lgndroi{2,R} num2str(r) '-' lgndroi{2,R} num2str(rr)];
                    tf = strcmp(x, labelroiALL{1,R}(1,:));
                    idx = find(tf);
                    MATpMroi(r,rr)= pMroi(1,idx); %Mettre les  p sous forme de matrice
                end
            end

            MATpMroi(length(MATpMroi),length(MATpMroi)) = 0; %Ajouter le dernier 0 de la diagonale
            MATpMroi = MATpMroi + triu(MATpMroi,1)'; %Répliquer la moitié inférieure de la matrice
            tf = MATpMroi <=p & MATpMroi > 0; %Identifier les p significatifs
            MATsigG1G2Mroi = MATroimeanG1G2{1,R}; %Reprendre la matrice de connectivité du ROI
            MATsigG1G2Mroi(~tf) = 0; % Retirer les valeurs non sig

            clear r rr x tf idx

            MATneg = MATsigG1G2Mroi < 0;
            MAT = MATsigG1G2Mroi.*(MATneg);%DATA{id}.MAT %loader matrice
            if find(MAT)
            plotLst = roi{1,R}(2,:);
            label =  roi{1,R}(1,:);
            %fileorderconnectogramMroi  = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Mixte.txt';
            plotconnectogramroi(fileorderconnectogram{1,R},MAT,label,plotLst)
            savefig([savepath date '_NegConnectSigG1G2Mroi']);
            else
            end     

            MATpos = MATsigG1G2Mroi > 0;
            MAT = MATsigG1G2Mroi.*(MATpos);%DATA{id}.MAT %loader matrice
            if find(MAT)
            plotLst = roi{1,R}(2,:);
            label =  roi{1,R}(1,:);
            %fileorderconnectogramMroi  = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Mixte.txt';
            plotconnectogramroi(fileorderconnectogram{1,R},MAT,label,plotLst)
            savefig([savepath date '_PosConnectSigG1G2Mroi']);
            else
            end    
        end

        clear r x res sig A B C X p1


        %%RROI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%ANOVA%%%%%%%%%%%%%%%%%%%
        R = 2;
        x = 1;
        resRroi = [];
        for r = 1:(width(tbldataRroi)- 5)
            [pRroi(:,x),res] = anovan(table2array(tbldataRroi(:,r+5)),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
            close hidden;
            resRroi = [resRroi, res];
            x = x + 1;
        end
        % delete(findall(0));

        labelfactors = resRroi(2:end-2,1);
        for f = 1:numel(labelfactors)
            n_sig = sum(pRroi(f,:) <= p);
            fprintf('%d %s effects are significant p<=%.2f without correction for Rrois\n',n_sig,labelfactors{f,1}, p)
        end
        
        %%fdr%%%%%%%%%%%%
        if fdrmode == 1
            for f = 1:numel(labelfactors)
            % FDR according to Storey 2002%%%
            [fdr_Rroi,q_Rroi] = mafdr(pRroi(f,:)');
            n_sig = sum(q_Rroi <= p);
            fprintf('%d %s effects are significant p<=%.2f using FDR correction for Rroi\n',n_sig,labelfactors{f,1}, p)

            % FDR according to Benjamini & Hochberg(1995) de Groppe%%%%
            [h_Rroi1, crit_p_Rroi, adj_ci_cvrg_Rroi, adj_p_Rroi] = fdr_bh(pRroi(f,:), 0.05, 'pdep', 'no' );
            n_sig = sum(adj_p_Rroi <= p);
            fprintf('%d %s effects are significant p<=%.2f using FDR correction for Rroi\n',n_sig,labelfactors{f,1}, p)
    
            % % FDR according to Benjamini & Hochberg(1995) de Gerber%%%%
            % [ind, thres] = FDR(pRroi(f,:));

            %Bonferroni-Holm correction by Groppe%%% Corrected_p are greater than 1...
            %doesn't make sense
            [corrected_p_Rroi, h_Rroi2] = bonf_holm(pRroi(f,:));
            n_sig = sum(corrected_p_Rroi <= p);
            fprintf('%d %s effects are significant p<=%.2f using Bonferroni-Holm correction for Rroi\n',n_sig,labelfactors{f,1}, p)

            %Multcmp function to apply Holm Step Down Procedure%%
            [padjBH_Rroi,alphaBH_Rroi] = multicmp (pRroi(f,:)','down',0.05);
            n_sig = sum(padjBH_Rroi <= p);
            fprintf('%d %s effects are significant p<=%.2f using Holm Step Down Procedure for Rroi\n',n_sig,labelfactors{f,1}, p)

            %Multcmp function to apply Hochberg's step up procedure%%
            [padjH_Rroi,alphaH_Rroi] = multicmp (pRroi(f,:)','up',0.05);
            n_sig = sum(padjH_Rroi <= p);
            fprintf('%d %s effects are significant p<=%.2f using Hochberg Step Up Procedure for Rroi\n',n_sig,labelfactors{f,1}, p)

            %Multicmp function to apply FDR correction according to Benjamini & Hochberg%%
            [padjFDR_Rroi,alphaFDR_Rroi] = multicmp (pRroi(f,:)','fdr',0.05);
            n_sig = sum(padjFDR_Rroi <= p);
            fprintf('%d %s effects are significant p<=%.2f using FDR correction for Rroi\n',n_sig,labelfactors{f,1}, p)
            end
        end

        %%tableaux des résultats%%%%%%%%%%%%%%%%
        if tablegraphmode == 1
            tblresRroi = array2table(resRroi);
            for r = 1:(width(tbldataRroi)- 5)
                tblresRroi = mergevars(tblresRroi,(r:r+6));
            end
            tblresRroi.Properties.VariableNames = labelroiALL{1,2};
            tblresRroi.Properties.RowNames = {'Source','Group','SES','Error','Total'};

            tblpRroi = array2table(pRroi);
            tblpRroi.Properties.VariableNames = labelroiALL{1,2};
            tblpRroi.Properties.RowNames = {'Group','SES'};

            sig = tblpRroi.Variables <= p;
            tblpsigRroi = tblpRroi(:,any(sig));
            tblsigRroi = tblresRroi(:,tblpsigRroi.Properties.VariableNames);

            grsig = tblpRroi{{'Group'},:} <= p;
            tblgrpsigRroi = tblpRroi(1,grsig);
            tblgrsigRroi = tblresRroi(:, tblgrpsigRroi.Properties.VariableNames);

            disp(tblgrpsigRroi)
            %writetable(tblsigRroi,[savepath date '_sigresultsRroi.xls']); trop gros
            writetable(tblpsigRroi,[savepath date '_psigresultsRroi.xls'],'WriteRowNames',true);
            writetable(tblgrsigRroi,[savepath date '_grsigresultsRroi.xls']);
            save([savepath date '_resultsRroi.mat'],'tblresRroi','tblpRroi','tblpsigRroi','tblsigRroi','tblgrpsigRroi','tblgrsigRroi');

            %%%% graphique%%%%%%%%%%%%%%%%%%%%
            if find(grsig)
                figure
                A = meanroiG1r{1,R}(:,grsig);
                B = meanroiG2r{1,R}(:,grsig);
                C = [A' B'];
                X = categorical(tblgrpsigRroi.Properties.VariableNames);
                p1 = bar(X,C);
                p1(1).FaceColor = 'r';
                p1(2).FaceColor = 'b';
                ylabel('Pearson correlation');
                savefig([savepath date '_SigRroi']);
                exportgraphics(gcf,[savepath 'SigRroi.png'])
            else
            end

            %%%% Connectogramme%%%%%%%%%%%%
%             %meanG2G1Rroi = meanG2Rroi-meanG1Rroi; %Calculer matrice G2-G1
%             meanG1G2Rroi = meanG1Rroi-meanG2Rroi; %Calculer matrice G1-G2
% 
%             for r = 1:numel(roi{1,R}(2,:))
%                 for rr = (r + 1):numel(roi{1,R}(2,:))
%                     x = [lgndroi{2,R} num2str(r) '-' lgndroi{2,R} num2str(rr)];
%                     tf = strcmp(x, labelroiALL{1,R}(1,:));
%                     idx = find(tf);
%                     MATmeanG1G2Rroi(r,rr) = meanG1G2Rroi(1,idx); %Mettre les  p sous forme de matrice
%                 end
%             end
%             MATmeanG1G2Rroi(length(MATmeanG1G2Rroi),length(MATmeanG1G2Rroi)) = 0; %Ajouter le dernier 0 de la diagonale
%             MATmeanG1G2Rroi = MATmeanG1G2Rroi + triu(MATmeanG1G2Rroi,1)'; %Répliquer la moitié inférieure de la matrice

            for r = 1:numel(roi{1,R}(2,:))
                for rr = (r + 1):numel(roi{1,R}(2,:))
                    x = [lgndroi{2,R} num2str(r) '-' lgndroi{2,R} num2str(rr)];
                    tf = strcmp(x, labelroiALL{1,R}(1,:));
                    idx = find(tf);
                    MATpRroi(r,rr)= pRroi(1,idx); %Mettre les  p sous forme de matrice
                end
            end

            MATpRroi(length(MATpRroi),length(MATpRroi)) = 0; %Ajouter le dernier 0 de la diagonale
            MATpRroi = MATpRroi + triu(MATpRroi,1)'; %Répliquer la moitié inférieure de la matrice
            tf = MATpRroi <=p & MATpRroi > 0; %Identifier les p significatifs
            MATsigG1G2Rroi = MATroimeanG1G2{1,R}; %Reprendre la matrice de connectivité du ROI
            MATsigG1G2Rroi(~tf) = 0; % Retirer les valeurs non sig

            clear r rr x tf idx

            MATneg = MATsigG1G2Rroi < 0;
            MAT = MATsigG1G2Rroi.*(MATneg);%DATA{id}.MAT %loader matrice
            if find(MAT)
            plotLst = roi{1,R}(2,:);
            label =  roi{1,R}(1,:);
            %fileorderconnectogramRroi  = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Region.txt';
            plotconnectogramroi(fileorderconnectogram{1,R},MAT,label,plotLst)
            savefig([savepath date '_NegConnectSigG1G2Rroi']);
            else
            end     

            MATpos = MATsigG1G2Rroi > 0;
            MAT = MATsigG1G2Rroi.*(MATpos);%DATA{id}.MAT %loader matrice
            if find(MAT)
            plotLst = roi{1,R}(2,:);
            label =  roi{1,R}(1,:);
            %fileorderconnectogramRroi  = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Region.txt';
            plotconnectogramroi(fileorderconnectogram{1,R},MAT,label,plotLst)
            savefig([savepath date '_PosConnectSigG1G2Rroi']);
            else
            end  

        end

        clear r x res sig A B C X p1


        %%FROI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%ANOVA%%%%%%%%%%%%%%%%%%%
        R = 3;
        x = 1;
        resFroi = [];
        for r = 1:(width(tbldataFroi)- 5)
            [pFroi(:,x),res] = anovan(table2array(tbldataFroi(:,r+5)),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
            close hidden;
            resFroi = [resFroi, res];
            x = x + 1;
        end
        % delete(findall(0));

        labelfactors = resFroi(2:end-2,1);
        for f = 1:numel(labelfactors)
            n_sig = sum(pFroi(f,:) <= p);
            fprintf('%d %s effects are significant p<=%.2f without correction for Froi\n',n_sig,labelfactors{f,1}, p)
        end
        
        %%fdr%%%%%%%%%%%%
        if fdrmode == 1
            for f = 1:numel(labelfactors)
            % FDR according to Storey 2002%%%
            [fdr_Froi,q_Froi] = mafdr(pFroi(f,:)');
            n_sig = sum(q_Froi <= p);
            fprintf('%d %s effects are significant p<=%.2f using FDR correction for Froi\n',n_sig,labelfactors{f,1}, p)

            % FDR according to Benjamini & Hochberg(1995) de Groppe%%%%
            [h_Froi1, crit_p_Froi, adj_ci_cvrg_Froi, adj_p_Froi] = fdr_bh(pFroi(f,:), 0.05, 'pdep', 'no' );
            n_sig = sum(adj_p_Froi <= p);
            fprintf('%d %s effects are significant p<=%.2f using FDR correction for Froi\n',n_sig,labelfactors{f,1}, p)

            % % FDR according to Benjamini & Hochberg(1995) de Gerber%%%%
            % [ind, thres] = FDR(pFroi(f,:));

            %Bonferroni-Holm correction by Groppe%%% Corrected_p are greater than 1...
            %doesn't make sense
            [corrected_p_Froi, h_Froi2] = bonf_holm(pFroi(f,:));
            n_sig = sum(corrected_p_Froi <= p);
            fprintf('%d %s effects are significant p<=%.2f using Bonferroni-Holm correction for Froi\n',n_sig,labelfactors{f,1}, p)

            %Multcmp function to apply Holm Step Down Procedure%%
            [padjBH_Froi,alphaBH_Froi] = multicmp (pFroi(f,:)','down',0.05);
            n_sig = sum(padjBH_Froi <= p);
            fprintf('%d %s effects are significant p<=%.2f using Holm Step Down Procedure for Froi\n',n_sig,labelfactors{f,1}, p)

            %Multcmp function to apply Hochberg's step up procedure%%
            [padjH_Froi,alphaH_Froi] = multicmp (pFroi(f,:)','up',0.05);
            n_sig = sum(padjH_Froi <= p);
            fprintf('%d %s effects are significant p<=%.2f using Hochberg Step Up Procedure for Froi\n',n_sig,labelfactors{f,1}, p)

            %Multicmp function to apply FDR correction according to Benjamini & Hochberg%%
            [padjFDR_Froi,alphaFDR_Froi] = multicmp (pFroi(f,:)','fdr',0.05);
            n_sig = sum(padjFDR_Froi <= p);
            fprintf('%d %s effects are significant p<=%.2f using FDR correction for Froi\n',n_sig,labelfactors{f,1}, p)
            end
        end

        %%tableaux des résultats%%%%%%%%%%%%%%%%
        if tablegraphmode == 1
            tblresFroi = array2table(resFroi);
            for r = 1:(width(tbldataFroi)- 5)
                tblresFroi = mergevars(tblresFroi,(r:r+6));
            end
            tblresFroi.Properties.VariableNames = labelroiALL{1,R};
            tblresFroi.Properties.RowNames = {'Source','Group','SES','Error','Total'};

            tblpFroi = array2table(pFroi);
            tblpFroi.Properties.VariableNames = labelroiALL{1,R};
            tblpFroi.Properties.RowNames = {'Group','SES'};

            sig = tblpFroi.Variables <= p;
            tblpsigFroi = tblpFroi(:,any(sig));
            tblsigFroi = tblresFroi(:,tblpsigFroi.Properties.VariableNames);

            grsig = tblpFroi{{'Group'},:} <= p;
            tblgrpsigFroi = tblpFroi(1,grsig);
            tblgrsigFroi = tblresFroi(:, tblgrpsigFroi.Properties.VariableNames);

            disp(tblgrpsigFroi)
            writetable(tblsigFroi,[savepath date '_sigresultsFroi.xls']);
            writetable(tblgrsigFroi,[savepath date '_grsigresultsFroi.xls']);
            save([savepath date '_resultsFroi.mat'],'tblresFroi','tblpFroi','tblpsigFroi','tblsigFroi','tblgrpsigFroi','tblgrsigFroi');

            %%%% graphique%%%%%%%%%%%%%%%%%%%%
            if find(grsig)
                figure
                A = meanroiG1r{1,R}(:,grsig);
                B = meanroiG2r{1,R}(:,grsig);
                C = [A' B'];
                X = categorical(tblgrpsigFroi.Properties.VariableNames);
                p1 = bar(X,C);
                p1(1).FaceColor = 'r';
                p1(2).FaceColor = 'b';
                ylabel('Pearson correlation');
                savefig([savepath date '_SigFroi']);
                exportgraphics(gcf,[savepath 'SigFroi.png'])
            else
            end

            %%%% Connectogramme%%%%%%%%%%%%
%             %meanG2G1Froi = meanroiG2r{1,R}-meanroiG1r{1,R}; %Calculer matrice G2-G1
%             meanG1G2Froi = meanroiG1r{1,R}-meanroiG2r{1,R}; %Calculer matrice G1-G2
% 
%             for r = 1:numel(roi{1,R}(2,:))
%                 for rr = (r + 1):numel(roi{1,R}(2,:))
%                     x = [lgndroi{2,R} num2str(r) '-' lgndroi{2,R} num2str(rr)];
%                     tf = strcmp(x, labelroiALL{1,R}(1,:));
%                     idx = find(tf);
%                     matroimeanG1G2{1,R}(r,rr) = meanG1G2Froi(1,idx); %Mettre sous forme de matrice
%                 end
%             end
%             matroimeanG1G2{1,R}(length(matroimeanG1G2{1,R}),length(matroimeanG1G2{1,R})) = 0; %Ajouter le dernier 0 de la diagonale
%             matroimeanG1G2{1,R} = matroimeanG1G2{1,R} + triu(matroimeanG1G2{1,R},1)'; %Répliquer la moitié inférieure de la matrice

            for r = 1:numel(roi{1,R}(2,:))
                for rr = (r + 1):numel(roi{1,R}(2,:))
                    x = [lgndroi{2,R} num2str(r) '-' lgndroi{2,R} num2str(rr)];
                    tf = strcmp(x, labelroiALL{1,R}(1,:));
                    idx = find(tf);
                    MATpFroi(r,rr)= pFroi(1,idx); %Mettre les  p sous forme de matrice
                end
            end

            MATpFroi(length(MATpFroi),length(MATpFroi)) = 0; %Ajouter le dernier 0 de la diagonale
            MATpFroi = MATpFroi + triu(MATpFroi,1)'; %Répliquer la moitié inférieure de la matrice
            tf = MATpFroi <=p & MATpFroi > 0; %Identifier les p significatifs
            MATsigG1G2Froi = MATroimeanG1G2{1,R}; %Reprendre la matrice de connectivité du ROI
            MATsigG1G2Froi(~tf) = 0; % Retirer les valeurs non sig

            clear r rr x tf idx

            MATneg = MATsigG1G2Froi < 0;
            MAT = MATsigG1G2Froi.*(MATneg);%DATA{id}.MAT %loader matrice
            if find(MAT)
            plotLst = roi{1,R}(2,:);
            label =  roi{1,R}(1,:);
            %fileorderconnectogramFroi  = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Fonction.txt';
            plotconnectogramroi(fileorderconnectogram{1,R},MAT,label,plotLst)
            savefig([savepath date '_NegConnectSigG1G2Froi']);
            else
            end     

            MATpos = MATsigG1G2Froi > 0;
            MAT = MATsigG1G2Froi.*(MATpos);%DATA{id}.MAT %loader matrice
            if find(MAT)
            plotLst = roi{1,R}(2,:);
            label =  roi{1,R}(1,:);
            %fileorderconnectogramFroi  = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Fonction.txt';
            plotconnectogramroi(fileorderconnectogram{1,R},MAT,label,plotLst)
            savefig([savepath date '_PosConnectSigG1G2Froi']);
            else
            end  

        end

        clear r x res sig A B C X p1


        %%AROI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%ANOVA%%%%%%%%%%%%%%%%%%%
        R = 4;
        x = 1;
        resAroi = [];
        for r = 1:(width(tbldataAroi)- 5)
            [pAroi(:,x),res] = anovan(table2array(tbldataAroi(:,r+5)),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
            close hidden;
            resAroi = [resAroi, res];
            x = x + 1;
        end
        % delete(findall(0));

        labelfactors = resAroi(2:end-2,1);
        for f = 1:numel(labelfactors)
            n_sig = sum(pAroi(f,:) <= p);
            fprintf('%d %s effects are significant p<=%.2f without correction for Aroi\n',n_sig,labelfactors{f,1}, p)
        end

        %%%%fdr%%%%%%%%%%%%%
        if fdrmode == 1
            for f = 1:numel(labelfactors)
            % FDR according to Storey 2002%%%
            [fdr_Aroi,q_Aroi] = mafdr(pAroi(f,:)');
            n_sig = sum(q_Aroi <= p);
            fprintf('%d %s effects are significant p<=%.2f using FDR correction for Aroi\n',n_sig,labelfactors{f,1}, p)

            % FDR according to Benjamini & Hochberg(1995) de Groppe%%%%
            [h_Aroi1, crit_p_Aroi, adj_ci_cvrg_Aroi, adj_p_Aroi] = fdr_bh(pAroi(f,:)', 0.05, 'pdep', 'no' );
            n_sig = sum(adj_p_Aroi <= p);
            fprintf('%d %s effects are significant p<=%.2f using FDR correction for Aroi\n',n_sig,labelfactors{f,1}, p)

            % % FDR according to Benjamini & Hochberg(1995) de Gerber%%%%
            %[ind, thres] = FDR(pAroi(f,:));

            %Bonferroni-Holm correction by Groppe%%% Corrected_p are greater than 1... doesn't make sense
            [corrected_p_Aroi, h_Aroi2] = bonf_holm(pAroi(f,:)', 0.05);
            n_sig = sum(corrected_p_Aroi <= p);
            fprintf('%d %s effects are significant p<=%.2f using Bonferroni-Holm correction for Aroi\n',n_sig,labelfactors{f,1}, p)

            %Multicmp function to apply Holm Step Down Procedure%%
            [padjBH_Aroi,alphaBH_Aroi] = multicmp (pAroi(f,:)','down',0.05);
            n_sig = sum(padjBH_Aroi <= p);
            fprintf('%d %s effects are significant p<=%.2f using Holm Step Down Procedure for Aroi\n',n_sig,labelfactors{f,1}, p)

            %Multicmp function to apply Hochberg's step up procedure%%
            [padjH_Aroi,alphaH_Aroi] = multicmp (pAroi(f,:)','up',0.05);
            n_sig = sum(padjH_Aroi <= p);
            fprintf('%d %s effects are significant p<=%.2f using Hochberg Step Up Procedure for Aroi\n',n_sig,labelfactors{f,1}, p)

            %Multicmp function to apply FDR correction according to Benjamini & Hochberg%%
            [padjFDR_Aroi,alphaFDR_Aroi] = multicmp (pAroi(f,:)','fdr',0.05);
            n_sig = sum(padjFDR_Aroi <= p);
            fprintf('%d %s effects are significant p<=%.2f using FDR correction for Aroi\n',n_sig,labelfactors{f,1}, p)
            end
        end

        %%tableaux des résultats%%%%%%%%%%%%%%%%
        if tablegraphmode == 1
            tblresAroi = array2table(resAroi);
            for r = 1:(width(tbldataAroi)- 5)
                tblresAroi = mergevars(tblresAroi,(r:r+6));
            end
            tblresAroi.Properties.VariableNames = labelroiALL{1,R};
            tblresAroi.Properties.RowNames = {'Source','Group','SES','Error','Total'};

            tblpAroi = array2table(pAroi);
            tblpAroi.Properties.VariableNames = labelroiALL{1,R};
            tblpAroi.Properties.RowNames = {'Group','SES'};

            sig = tblpAroi.Variables <= p;
            tblpsigAroi = tblpAroi(:,any(sig));
            tblsigAroi = tblresAroi(:,tblpsigAroi.Properties.VariableNames);

            grsig = tblpAroi{{'Group'},:} <= p;
            tblgrpsigAroi = tblpAroi(1,grsig);
            tblgrsigAroi = tblresAroi(:, tblgrpsigAroi.Properties.VariableNames);

            disp(tblgrpsigAroi)
            writetable(tblsigAroi,[savepath date '_sigresultsAroi.xls']);
            writetable(tblgrsigAroi,[savepath date '_grsigresultsAroi.xls']);
            save([savepath date '_resultsAroi.mat'],'tblresAroi','tblpAroi','tblpsigAroi','tblsigAroi','tblgrpsigAroi','tblgrsigAroi');

            %%%% graphique%%%%%%%%%%%%%%%%%%%%
            if find(grsig)
                figure
                A = meanroiG1r{1,R}(:,grsig);
                B = meanroiG2r{1,R}(:,grsig);
                C = [A' B'];
                X = categorical(tblgrpsigAroi.Properties.VariableNames);
                p1 = bar(X,C);
                p1(1).FaceColor = 'r';
                p1(2).FaceColor = 'b';
                ylabel('Pearson correlation');
                savefig([savepath date '_SigAroi']);
                exportgraphics(gcf,[savepath 'SigAroi.png'])
            else
            end

                %%%% Connectogramme%%%%%%%%%%%%
%             %meanG2G1Aroi = meanroiG2r{1,R}-meanroiG1r{1,R}; %Calculer matrice G2-G1
%             meanG1G2Aroi = meanroiG1r{1,R}-meanroiG2r{1,R}; %Calculer matrice G2-G1
% 
%             for r = 1:numel(roi{1,R}(2,:))
%                 for rr = (r + 1):numel(roi{1,R}(2,:))
%                     x = [lgndroi{2,R} num2str(r) '-' lgndroi{2,R} num2str(rr)];
%                     tf = strcmp(x, labelroiALL{1,R}(1,:));
%                     idx = find(tf);
%                     MATroimeanG1G2{1,R}(r,rr) = meanG1G2Aroi(1,idx); %Mettre les données sous forme de matrice
%                 end
%             end
%             MATroimeanG1G2{1,R}(length(MATroimeanG1G2{1,R}),length(MATroimeanG1G2{1,R})) = 0; %Ajouter le dernier 0 de la diagonale
%             MATroimeanG1G2{1,R} = MATroimeanG1G2{1,R} + triu(MATroimeanG1G2{1,R},1)'; %Répliquer la moitié inférieure de la matrice

            for r = 1:numel(roi{1,R}(2,:))
                for rr = (r + 1):numel(roi{1,R}(2,:))
                    x = [lgndroi{2,R} num2str(r) '-' lgndroi{2,R} num2str(rr)];
                    tf = strcmp(x, labelroiALL{1,R}(1,:));
                    idx = find(tf);
                    MATpAroi(r,rr)= pAroi(1,idx); %Mettre les  p sous forme de matrice
                end
            end

            MATpAroi(length(MATpAroi),length(MATpAroi)) = 0; %Ajouter le dernier 0 de la diagonale
            MATpAroi = MATpAroi + triu(MATpAroi,1)'; %Répliquer la moitié inférieure de la matrice
            tf = MATpAroi <=p & MATpAroi > 0; %Identifier les p significatifs
            MATsigG1G2Aroi = MATroimeanG1G2{1,R}; %Reprendre la matrice de connectivité du ROI
            MATsigG1G2Aroi(~tf) = 0; % Retirer les valeurs non sig

            clear r rr x tf idx

            MATneg = MATsigG1G2Aroi < 0;
            MAT = MATsigG1G2Aroi.*(MATneg);%DATA{id}.MAT %loader matrice
            if find(MAT)
            plotLst = roi{1,R}(2,:);
            label =  roi{1,R}(1,:);
            %fileorderconnectogramAroi  = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Aire.txt';
            plotconnectogramroi(fileorderconnectogram{1,R},MAT,label,plotLst)
            savefig([savepath date '_NegConnectSigG1G2Aroi']);
            else
            end     

            MATpos = MATsigG1G2Aroi > 0;
            MAT = MATsigG1G2Aroi.*(MATpos);%DATA{id}.MAT %loader matrice
            if find(MAT)
            plotLst = roi{1,R}(2,:);
            label =  roi{1,R}(1,:);
            %fileorderconnectogramAroi  = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Aire.txt';
            plotconnectogramroi(fileorderconnectogram{1,R},MAT,label,plotLst)
            savefig([savepath date '_PosConnectSigG1G2Aroi']);
            else
            end  

        end

        clear r x res sig A B C X p1

        if tablegraphmode == 1
            tblgrpsigALL = [tblgrpsigMroi, tblgrpsigRroi, tblgrpsigFroi, tblgrpsigAroi];
            writetable(tblgrpsigALL,[savepath date '_grpsigresultsALL.xls'],'WriteRowNames',true);
        end
   end
X = ['Results saved in ', savepath];
disp(X)
clear X
end

clearvars -except datapath savepathinitial channelmode importedROImode calculatedROImode tablegraphmode fileorderconnectogram fdrmode nointeractionmode interactionmode p
%load ([datapath 'workspace.mat'])
%load ([datapath 'workspacemat.mat'])

%%
%ANCOVA AVEC INTERACTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if interactionmode == 1
disp('Computing ANCOVA with interaction')
savepath = [savepathinitial 'Interaction\'];
if ~isfolder(savepath)
    mkdir(savepath)
end

%%CHANNELS%%%%%%%%%%%%%%%%%%%%%
    if channelmode ==1
        x = 1;
        resch = [];
        for r = 1:(width(tblch)- 5)
            [pch(:,x),res] = anovan(table2array(tblch(:,r+5)),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
            close hidden;
            resch = [resch, res];
            x = x + 1;
        end
        % delete(findall(0));

        labelfactors = resch(2:end-2,1);
        for f = 1:numel(labelfactors)
            n_sig = sum(pch(f,:) <= p);
            fprintf('%d %s effects are significant p<=%.2f without correction for channels\n',n_sig,labelfactors{f,1}, p)
        end

        clear r x res n_sig
       
       %%fdr%%%%%%%%%%%%
        if fdrmode == 1
            for f = 1:numel(labelfactors)
                % % FDR according to Storey 2002%%%
                [fdr_ch,q_ch] = mafdr(pch(f,:)');
                n_sig = sum(q_ch <= p);
                fprintf('%d %s effects are significant p<=%.2f using FDR correction for channels\n',n_sig,labelfactors{f,1}, p)

                % FDR according to Benjamini & Hochberg(1995) de Groppe%%%%
                [h_ch1, crit_p_ch, adj_ci_cvrg_ch, adj_p_ch] = fdr_bh(pch(f,:), 0.05, 'pdep', 'no' );
                n_sig = sum(adj_p_ch(1,:) <= p);
                fprintf('%d %s effects are significant p<=%.2f using FDR correction for channels\n',n_sig,labelfactors{f,1}, p)

                % % FDR according to Benjamini & Hochberg(1995) de Gerber%%%%
                % [ind, thres] = FDR(pch(f,:));

                %Bonferroni-Holm correction by Groppe%%% Corrected_p are greater than 1...
                %doesn't make sense
                [corrected_p_ch, h_ch2] = bonf_holm(pch(f,:));
                n_sig = sum(corrected_p_ch <= p);
                fprintf('%d %s effects are significant p<=%.2f using Bonferroni-Holm correction for channels\n',n_sig,labelfactors{f,1}, p)

                %Multcmp function to apply Holm Step Down Procedure%%
                [padjBH_ch,alphaBH_ch] = multicmp (pch(f,:)','down',0.05);
                n_sig = sum(padjBH_ch <= p);
                fprintf('%d %s effects are significant p<=%.2f using Holm Step Down Procedure for channels\n',n_sig,labelfactors{f,1}, p)

                %Multcmp function to apply Hochberg's step up procedure%%
                [padjH_ch,alphaH_ch] = multicmp (pch(f,:)','up',0.05);
                n_sig = sum(padjH_ch <= p);
                fprintf('%d %s effects are significant p<=%.2f using Hochberg Step Up Procedure for channels\n',n_sig,labelfactors{f,1}, p)

                %Multicmp function to apply FDR correction according to Benjamini & Hochberg%%
                [padjFDR_ch,alphaFDR_ch] = multicmp (pch(f,:)','fdr',0.05);
                n_sig = sum(padjFDR_ch <= p);
                fprintf('%d %s effects are significant p<=%.2f using FDR correction for channels\n',n_sig,labelfactors{f,1}, p)
            end 
        end
        
        %%tableaux des résultats%%%%%%%%%%%%%%%%
        if tablegraphmode == 1
            tblresch = array2table(resch);
            for r = 1:(width(tblch)- 5)
                tblresch = mergevars(tblresch,(r:r+6));
            end
            tblresch.Properties.VariableNames = labelch;
            tblresch.Properties.RowNames = {'Source','Group','SES','Group*SES','Error','Total'};

            tblpch = array2table(pch);
            tblpch.Properties.VariableNames = labelch;
            tblpch.Properties.RowNames = {'Group','SES','Group*SES'};

            sig = tblpch.Variables <= p;
            tblpsigch = tblpch(:,any(sig));
            tblsigch = tblresch(:,tblpsigch.Properties.VariableNames);

            grsig = tblpch{{'Group'},:} <= p;
            tblgrpsigch = tblpch(1,grsig);
            tblgrsigch = tblresch(:, tblgrpsigch.Properties.VariableNames);

            grsessig = tblpch{{'Group*SES'},:} <= p;
            tblgrsespsigch = tblpch(3,grsessig);
            tblgrsessigch = tblresch(:, tblgrsespsigch.Properties.VariableNames);

            disp(tblgrpsigch)
            disp(tblgrsespsigch)

            %writetable(tblsigch,[savepath date '_sigresultsch.xls']); trop gros
            writetable(tblpsigch,[savepath date '_psigresultsch.xls'],'WriteRowNames',true);
            %writetable(tblgrsigch,[savepath date '_grsigresultsch.xls']);
            save([savepath date '_resultsch.mat'],'tblresch','tblpch','tblpsigch','tblsigch','tblgrpsigch','tblgrsigch');

            clear r

            %%%% graphique%%%%%%%%%%%%%%%%%%%%
            if find(grsig)
                figure
                A = meanG1ch(:,grsig);
                B = meanG2ch(:,grsig);
                C = [A' B'];
                X = categorical(tblgrpsigch.Properties.VariableNames);
                p1 = bar(X,C);
                p1(1).FaceColor = 'r';
                p1(2).FaceColor = 'b';
                ylabel('Pearson correlation');
                savefig([savepath date '_SigCh']);
                exportgraphics(gcf,[savepath 'SigCh.png'])

            % figure
            % p1 = bar(meanG1ch(:,grsig),'r');
            % hold on
            % p2 = bar(meanG2ch(:,grsig),'b');
            % xlabel(tblgrpsigch.Properties.VariableNames);
            % ylabel('Pearson correlation');
            % hold off
            else
            end

            clear A B C X p1 sig grsig

            %%%% Connectogramme%%%%%%%%%%%%
%             %MATmeanG2G1 = MATmeanG2-MATmeanG1; %Calculer matrice G2-G1
%             MATmeanG1G2 = MATmeanG1 - MATmeanG2;

            for c=1:46
                for cc =(c + 1):46
                    x = [num2str(c) '-' num2str(cc)];
                    tf = strcmp(x, labelch);
                    idx = find(tf);
                    MATpch(c,cc)= pch(1,idx); %Mettre les  p sous forme de matrice
                end
            end

            MATpch(length(MATpch),length(MATpch)) = 0; 
            MATpch = MATpch + triu(MATpch,1)';
            tf = MATpch <=p & MATpch > 0;
    %         idx = find(tf);
            %MATsigG2G1ch = MATmeanG2G1;
            MATsigG1G2ch = MATmeanG1G2;
            %MATsigG2G1ch(~tf) = 0;
            MATsigG1G2ch(~tf) = 0;

            clear c cc x tf idx

            id = 1;
    %         MAT = MATsigch; %DATA{id}.MAT %loader matrice
    %         List = strvcat(DATA{id}.ZoneList); %liste des paires SD
    %         ML = DATA{id}.zone.ml; %Loader S/D/ROI/Gr
    %         plotLst = DATA{id}.zone.plotLst;
    %         label =  DATA{id}.zone.label;
    %         fileorderconnectogram  = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Mixte.txt';
    %         plotconnectogram(fileorderconnectogram,MAT,List,label,plotLst,ML)
    %         savefig([savepath date '_ConnectSigCh']);

            MATneg = MATsigG1G2ch < 0;
            MAT = MATsigG1G2ch.*(MATneg);%DATA{id}.MAT %loader matrice
            if find(MAT)
                List = strvcat(DATA{id}.ZoneList); %liste des paires SD
                ML = DATA{id}.zone.ml; %Loader S/D/ROI/Gr
                plotLst = DATA{id}.zone.plotLst;
                label =  DATA{id}.zone.label;
                plotconnectogram(fileorderconnectogram{1,1},MAT,List,label,plotLst,ML)
                savefig([savepath date '_NegConnectSigG1G2Ch']);
            else
            end

            MATpos = MATsigG1G2ch > 0;
            MAT = MATsigG1G2ch.*(MATpos);%DATA{id}.MAT %loader matrice
            if find(MAT)
                List = strvcat(DATA{id}.ZoneList); %liste des paires SD
                ML = DATA{id}.zone.ml; %Loader S/D/ROI/Gr
                plotLst = DATA{id}.zone.plotLst;
                label =  DATA{id}.zone.label;
                plotconnectogram(fileorderconnectogram{1,1},MAT,List,label,plotLst,ML)
                savefig([savepath date '_PosConnectSigG1G2Ch']);
            else
            end
        end
        clear r x res sig A B C X p1
    end

         if importedROImode
           %%ANOVA%%%%%%%%%%%%%%%%%%%
            x = 1;
            resroi = [];
            for r = 1:(width(tblroi)- 5)
                [proi(:,x),res] = anovan(table2array(tblroi(:,r+5)),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
                close hidden;
                resroi = [resroi, res];
                x = x + 1;
            end
            % delete(findall(0));

            labelfactors = resroi(2:end-2,1);
            for f = 1:numel(labelfactors)
                n_sig = sum(proi(f,:) <= p);
                fprintf('%d %s effects are significant p<=%.2f without correction for ROIs\n',n_sig,labelfactors{f,1}, p)
            end

            clear r x res n_sig

            %%fdr%%%%%%%%%%%%
            if fdrmode == 1;
                for f = 1:numel(labelfactors)
                     % FDR according to Storey 2002%%%
                     [fdr_roi,q_roi] = mafdr(proi(f,:)');
                     n_sig = sum(q_roi <= p);
                     fprintf('%d %s effects are significant p<=%.2f using FDR correction (Storey2002) for ROIs\n',n_sig,labelfactors{f,1}, p)

                    % FDR according to Benjamini & Hochberg(1995) de Groppe%%%%
                    [h_roi1, crit_p_roi, adj_ci_cvrg_roi, adj_p_roi] = fdr_bh(proi(f,:), 0.05, 'pdep', 'no' );
                    n_sig = sum(adj_p_roi(1,:) <= p);
                    fprintf('%d %s effects are significant p<=%.2f using FDR correction (Benjamini1995) for ROIs\n',n_sig,labelfactors{f,1}, p)

                    % % FDR according to Benjamini & Hochberg(1995) de Gerber%%%%
                    % [ind, thres] = FDR(pch(f,:));

                    %Bonferroni-Holm correction by Groppe%%% Corrected_p are greater than 1...
                    %doesn't make sense
                    [corrected_p_roi, h_roi2] = bonf_holm(proi(f,:));
                    n_sig = sum(corrected_p_roi <= p);
                    fprintf('%d %s effects are significant p<=%.2f using Bonferroni-Holm correction for ROIs\n',n_sig,labelfactors{f,1}, p)

                    %Multcmp function to apply Holm Step Down Procedure%%
                    [padjBH_roi,alphaBH_roi] = multicmp (proi(f,:)','down',0.05);
                    n_sig = sum(padjBH_roi <= p);
                    fprintf('%d %s effects are significant p<=%.2f using Holm Step Down Procedure for ROIs\n',n_sig,labelfactors{f,1}, p)

                    %Multcmp function to apply Hochberg's step up procedure%%
                    [padjH_roi,alphaH_roi] = multicmp (proi(f,:)','up',0.05);
                    n_sig = sum(padjH_roi <= p);
                    fprintf('%d %s effects are significant p<=%.2f using Hochberg Step Up Procedure for ROIs\n',n_sig,labelfactors{f,1}, p)

                    %Multicmp function to apply FDR correction according to Benjamini & Hochberg%%
                    [padjFDR_roi,alphaFDR_roi] = multicmp (proi(f,:)','fdr',0.05);
                    n_sig = sum(padjFDR_roi <= p);
                    fprintf('%d %s effects are significant p<=%.2f using FDR correction for ROIs\n',n_sig,labelfactors{f,1}, p)
                end
            end
            
            %%tableaux des résultats%%%%%%%%%%%%%%%%
            if tablegraphmode == 1
                tblresroi = array2table(resroi);
                for r = 1:(width(tblroi)- 5)
                    tblresroi = mergevars(tblresroi,(r:r+6));
                end
                tblresroi.Properties.VariableNames = labelroiALL;
                tblresroi.Properties.RowNames = {'Source','Group','SES','Group*SES','Error','Total'};

                tblproi = array2table(proi);
                tblproi.Properties.VariableNames = labelroiALL;
                tblproi.Properties.RowNames = {'Group','SES','Group*SES'};

                sig = tblproi.Variables <= p;
                tblpsigroi = tblproi(:,any(sig));
                tblsigroi = tblresroi(:,tblpsigroi.Properties.VariableNames);

                grsig = tblproi{{'Group'},:} <= p;
                tblgrpsigroi = tblproi(1,grsig);
                tblgrsigroi = tblresroi(:, tblgrpsigroi.Properties.VariableNames);

                disp(tblgrpsigroi)

                grsessig = tblproi{{'Group*SES'},:} <= p;
                tblgrsespsigroi = tblproi(3,grsessig);
                tblgrsessigroi = tblresroi(:, tblgrsespsigroi.Properties.VariableNames);

                disp(tblgrsespsigroi)
                
                %writetable(tblsigroi,[savepath date '_sigresultsroi.xls']); trop gros
                writetable(tblpsigroi,[savepath date '_psigresultsroi.xls'],'WriteRowNames',true);
                %writetable(tblgrsigroi,[savepath date '_grsigresultsroi.xls']);
                save([savepath date '_resultsroi.mat'],'tblresroi','tblproi','tblpsigroi','tblsigroi','tblgrpsigroi','tblgrsigroi');

                clear r

                %%%% graphique%%%%%%%%%%%%%%%%%%%%
                if find(grsig)
                    figure
                    A = meanG1roi(:,grsig);
                    B = meanG2roi(:,grsig);
                    C = [A' B'];
                    X = categorical(tblgrpsigroi.Properties.VariableNames);
                    p1 = bar(X,C);
                    p1(1).FaceColor = 'r';
                    p1(2).FaceColor = 'b';
                    ylabel('Pearson correlation');
                    savefig([savepath date '_SigROI']);
                    exportgraphics(gcf,[savepath 'SigROI.png'])

                else
                end

                clear A B C X p1 sig grsig

           %%%% Connectogramme%%%%%%%%%%%%
            %meanG2G1roi = meanG2roi-meanG1roi; %Calculer matrice G2-G1
            meanG1G2roi = meanG1roi-meanG2roi; %Calculer matrice G1-G2

            for r = 1:numel(roi(2,:))
                for rr = (r + 1):numel(roi(2,:))
                    x = ['R' num2str(r) '-' 'R' num2str(rr)];
                    tf = strcmp(x, labelroiALL(1,:));
                    idx = find(tf);
                    MATmeanG1G2roi(r,rr) = meanG1G2roi(1,idx); %Mettre les  données sous forme de matrice
                end
            end
            MATmeanG1G2roi(length(MATmeanG1G2roi),length(MATmeanG1G2roi)) = 0; %Ajouter le dernier 0 de la diagonale
            MATmeanG1G2roi = MATmeanG1G2roi + triu(MATmeanG1G2roi,1)'; %Répliquer la moitié inférieure de la matrice

            for r = 1:numel(roi(2,:))
                for rr = (r + 1):numel(roi(2,:))
                    x = ['R' num2str(r) '-' 'R' num2str(rr)];
                    tf = strcmp(x, labelroiALL(1,:));
                    idx = find(tf);
                    MATproi(r,rr)= proi(1,idx); %Mettre les  p sous forme de matrice
                end
            end

            MATproi(length(MATproi),length(MATproi)) = 0; %Ajouter le dernier 0 de la diagonale
            MATproi = MATproi + triu(MATproi,1)'; %Répliquer la moitié inférieure de la matrice
            tf = MATproi <=p & MATproi > 0; %Identifier les p significatifs
            MATsigG1G2roi = MATmeanG1G2roi; %Reprendre la matrice de connectivité du ROI
            MATsigG1G2roi(~tf) = 0; % Retirer les valeurs non sig

            clear r rr x tf idx

            MATneg = MATsigG1G2roi < 0;
            MAT = MATsigG1G2roi.*(MATneg);%DATA{id}.MAT %loader matrice
            if find(MAT)
            plotLst = roi(2,:);
            label =  roi(1,:);
            plotconnectogramroi(fileorderconnectogram,MAT,label,plotLst)
            savefig([savepath date '_NegConnectSigG1G2ROI']);
            else
            end     

            MATpos = MATsigG1G2roi > 0;
            MAT = MATsigG1G2roi.*(MATpos);%DATA{id}.MAT %loader matrice
            if find(MAT)
            plotLst = roi(2,:);
            label =  roi(1,:);
            plotconnectogramroi(fileorderconnectogram,MAT,label,plotLst)
            savefig([savepath date '_PosConnectSigG1G2ROI']);
            else
            end  

            end
        clear r x res sig A B C X p1
         end

    if calculatedROImode   
    
        %%MROI%%%%%%%%%%%%%%%%%%%%%%%
        R = 1;
        x = 1;
        resMroi = [];
        for r = 1:(width(tbldataMroi)- 5)
            [pMroi(:,x),res] = anovan(table2array(tbldataMroi(:,r+5)),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
            close hidden;
            resMroi = [resMroi, res];
            x = x + 1;
        end
        % delete(findall(0));

        labelfactors = resMroi(2:end-2,1);
        for f = 1:numel(labelfactors)
            n_sig = sum(pMroi(f,:) <= p);
            fprintf('%d %s effects are significant p<=%.2f without correction for Mroi\n',n_sig,labelfactors{f,1}, p)
        end
        
        %%fdr%%%%%%%%%%%%
        if fdrmode == 1
            for f = 1:numel(labelfactors)
            % FDR according to Storey 2002%%%
            [fdr_Mroi,q_Mroi] = mafdr(pMroi(f,:)');
            n_sig = sum(q_Mroi <= p);
            fprintf('%d %s effects are significant p <=%.2f using FDR correction for Mroi\n',n_sig,labelfactors{f,1}, p)

            % FDR according to Benjamini & Hochberg(1995) de Groppe%%%%
            [h_Mroi1, crit_p_Mroi, adj_ci_cvrg_Mroi, adj_p_Mroi] = fdr_bh(pMroi(f,:), 0.05, 'pdep', 'no' );
            n_sig = sum(adj_p_Mroi <= p);
            fprintf('%d %s effects are significant p<=%.2f using FDR correction for Mroi\n',n_sig,labelfactors{f,1}, p)

            % % FDR according to Benjamini & Hochberg(1995) de Gerber%%%%
            % [ind, thres] = FDR(pMroi(f,:));

            %Bonferroni-Holm correction by Groppe%%% Corrected_p are greater than 1...
            %doesn't make sense
            [corrected_p_Mroi, h_Mroi2] = bonf_holm(pMroi(f,:));
            n_sig = sum(corrected_p_Mroi <= p);
            fprintf('%d %s effects are significant p<=%.2f using Bonferroni-Holm correction for Mroi\n',n_sig,labelfactors{f,1}, p)

            %Multcmp function to apply Holm Step Down Procedure%%
            [padjBH_Mroi,alphaBH_Mroi] = multicmp (pMroi(f,:)','down',0.05);
            n_sig = sum(padjBH_Mroi <= p);
            fprintf('%d %s effects are significant p<=%.2f using Holm Step Down Procedure for Mroi\n',n_sig,labelfactors{f,1}, p)

            %Multcmp function to apply Hochberg's step up procedure%%
            [padjH_Mroi,alphaH_Mroi] = multicmp (pMroi(f,:)','up',0.05);
            n_sig = sum(padjH_Mroi <= p);
            fprintf('%d %s effects are significant p<=%.2f using Hochberg Step Up Procedure for Mroi\n',n_sig,labelfactors{f,1}, p)

            %Multicmp function to apply FDR correction according to Benjamini & Hochberg%%
            [padjFDR_Mroi,alphaFDR_Mroi] = multicmp (pMroi(f,:)','fdr',0.05);
            n_sig = sum(padjFDR_Mroi <= p);
            fprintf('%d %s effects are significant p<=%.2f using FDR correction for Mroi\n',n_sig,labelfactors{f,1}, p)
            end
        end

        %%tableaux des résultats%%%%%%%%%%%%%%%%
        if tablegraphmode == 1
            tblresMroi = array2table(resMroi);
            for r = 1:(width(tbldataMroi)- 5)
                tblresMroi = mergevars(tblresMroi,(r:r+6));
            end
            tblresMroi.Properties.VariableNames = labelroiALL{1,R};
            tblresMroi.Properties.RowNames = {'Source','Group','SES','Group*SES','Error','Total'};

            tblpMroi = array2table(pMroi);
            tblpMroi.Properties.VariableNames = labelroiALL{1,R};
            tblpMroi.Properties.RowNames = {'Group','SES','Group*SES'};

            sig = tblpMroi.Variables <= p;
            tblpsigMroi = tblpMroi(:,any(sig));
            tblsigMroi = tblresMroi(:,tblpsigMroi.Properties.VariableNames);

            grsig = tblpMroi{{'Group'},:} <= p;
            tblgrpsigMroi = tblpMroi(1,grsig);
            tblgrsigMroi = tblresMroi(:, tblgrpsigMroi.Properties.VariableNames);

            disp(tblgrpsigMroi)
            
            grsessig = tblpMroi{{'Group*SES'},:} <= p;
            tblgrsespsigMroi = tblpMroi(3,grsessig);
            tblgrsessigMroi = tblresMroi(:, tblgrsespsigMroi.Properties.VariableNames);

            disp(tblgrsespsigMroi)
            
            writetable(tblsigMroi,[savepath date '_sigresultsMroi.xls']);
            writetable(tblgrsigMroi,[savepath date '_grsigresultsMroi.xls']);
            save([savepath date '_resultsMroi.mat'],'tblresMroi','tblpMroi','tblpsigMroi', 'tblsigMroi','tblgrpsigMroi','tblgrsigMroi');

            %%%% graphique%%%%%%%%%%%%%%%%%%%%
            if find(grsig)
                figure
                A = meanroiG1r{1,R}(:,grsig);
                B = meanroiG2r{1,R}(:,grsig);
                C = [A' B'];
                X = categorical(tblgrpsigMroi.Properties.VariableNames);
                p1 = bar(X,C);
                p1(1).FaceColor = 'r';
                p1(2).FaceColor = 'b';
                ylabel('Pearson correlation');
                savefig([savepath date '_SigMroi']);
                exportgraphics(gcf,[savepath 'SigMroi.png'])
            else
            end
            % figure
            % p1 = bar(meanroiG1r{1,R}(:,grsig),'r');
            % hold on
            % p2 = bar(meanroiG2r{1,R}(:,grsig),'b');
            % xlabel(tblgrpsigMroi.Properties.VariableNames);
            % ylabel('Pearson correlation');
            % hold off

            %%%% Connectogramme%%%%%%%%%%%%
%             %meanG2G1Mroi = meanroiG2r{1,R}-meanroiG1r{1,R}; %Calculer matrice G2-G1
%             meanG1G2Mroi = meanroiG1r{1,R}-meanroiG2r{1,R}; %Calculer matrice G1-G2
% 
%             for r = 1:numel(roi{1,R}(2,:))
%                 for rr = (r + 1):numel(roi{1,R}(2,:))
%                     x = [lgndroi{2,R} num2str(r) '-' lgndroi{2,R} num2str(rr)];
%                     tf = strcmp(x, labelroiALL{1,R}(1,:));
%                     idx = find(tf);
%                     MATroimeanG1G2{1,R}(r,rr) = meanG1G2Mroi(1,idx); %Mettre les données sous forme de matrice
%                 end
%             end
%             MATroimeanG1G2{1,R}(length(MATroimeanG1G2{1,R}),length(MATroimeanG1G2{1,R})) = 0; %Ajouter le dernier 0 de la diagonale
%             MATroimeanG1G2{1,R} = MATroimeanG1G2{1,R} + triu(MATroimeanG1G2{1,R},1)'; %Répliquer la moitié inférieure de la matrice

            for r = 1:numel(roi{1,R}(2,:))
                for rr = (r + 1):numel(roi{1,R}(2,:))
                    x = [lgndroi{2,R} num2str(r) '-' lgndroi{2,R} num2str(rr)];
                    tf = strcmp(x, labelroiALL{1,R}(1,:));
                    idx = find(tf);
                    MATpMroi(r,rr)= pMroi(1,idx); %Mettre les  p sous forme de matrice
                end
            end

            MATpMroi(length(MATpMroi),length(MATpMroi)) = 0; %Ajouter le dernier 0 de la diagonale
            MATpMroi = MATpMroi + triu(MATpMroi,1)'; %Répliquer la moitié inférieure de la matrice
            tf = MATpMroi <=p & MATpMroi > 0; %Identifier les p significatifs
            MATsigG1G2Mroi = MATroimeanG1G2{1,R}; %Reprendre la matrice de connectivité du ROI
            MATsigG1G2Mroi(~tf) = 0; % Retirer les valeurs non sig

            clear r rr x tf idx

            MATneg = MATsigG1G2Mroi < 0;
            MAT = MATsigG1G2Mroi.*(MATneg);%DATA{id}.MAT %loader matrice
            if find(MAT)
            plotLst = roi{1,R}(2,:);
            label =  roi{1,R}(1,:);
            %fileorderconnectogramMroi  = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Mixte.txt';
            plotconnectogramroi(fileorderconnectogram{1,R},MAT,label,plotLst)
            savefig([savepath date '_NegConnectSigG1G2Mroi']);
            else
            end     

            MATpos = MATsigG1G2Mroi > 0;
            MAT = MATsigG1G2Mroi.*(MATpos);%DATA{id}.MAT %loader matrice
            if find(MAT)
            plotLst = roi{1,R}(2,:);
            label =  roi{1,R}(1,:);
            %fileorderconnectogramMroi  = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Mixte.txt';
            plotconnectogramroi(fileorderconnectogram{1,R},MAT,label,plotLst)
            savefig([savepath date '_PosConnectSigG1G2Mroi']);
            else
            end    
        end

        clear r x res sig A B C X p1


        %%RROI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        R = 2;
        x = 1;
        resRroi = [];
        for r = 1:(width(tbldataRroi)- 5)
            [pRroi(:,x),res] = anovan(table2array(tbldataRroi(:,r+5)),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
            close hidden;
            resRroi = [resRroi, res];
            x = x + 1;
        end
        % delete(findall(0));

         labelfactors = resRroi(2:end-2,1);
        for f = 1:numel(labelfactors)
            n_sig = sum(pRroi(f,:) <= p);
            fprintf('%d %s effects are significant p<=%.2f without correction for Rroi\n',n_sig,labelfactors{f,1}, p)
        end
        
        %%fdr%%%%%%%%%%%%
        if fdrmode == 1
            for f = 1:numel(labelfactors)
                % % FDR according to Storey 2002%%%
                [fdr_Rroi,q_Rroi] = mafdr(pRroi(f,:)');
                n_sig = sum(q_Rroi <= p);
                fprintf('%d %s effects are significant p<=%.2f using FDR correction for Rroi\n',n_sig,labelfactors{f,1}, p)

                % FDR according to Benjamini & Hochberg(1995) de Groppe%%%%
                [h_Rroi1, crit_p_Rroi, adj_ci_cvrg_Rroi, adj_p_Rroi] = fdr_bh(pRroi(f,:), 0.05, 'pdep', 'no' );
                n_sig = sum(adj_p_Rroi <= p);
                fprintf('%d %s effects are significant p<=%.2f using FDR correction for Rroi\n',n_sig,labelfactors{f,1}, p)

                % % FDR according to Benjamini & Hochberg(1995) de Gerber%%%%
                % [ind, thres] = FDR(pRroi(f,:));

                %Bonferroni-Holm correction by Groppe%%% Corrected_p are greater than 1...
                %doesn't make sense
                [corrected_p_Rroi, h_Rroi2] = bonf_holm(pRroi(f,:));
                n_sig = sum(corrected_p_Rroi <= p);
                fprintf('%d %s effects are significant p<=%.2f using Bonferroni-Holm correction for Rroi\n',n_sig,labelfactors{f,1}, p)

                %Multcmp function to apply Holm Step Down Procedure%%
                [padjBH_Rroi,alphaBH_Rroi] = multicmp (pRroi(f,:)','down',0.05);
                n_sig = sum(padjBH_Rroi <= p);
                fprintf('%d %s effects are significant p<=%.2f using Holm Step Down Procedure for Rroi\n',n_sig,labelfactors{f,1}, p)

                %Multcmp function to apply Hochberg's step up procedure%%
                [padjH_Rroi,alphaH_Rroi] = multicmp (pRroi(f,:)','up',0.05);
                n_sig = sum(padjH_Rroi <= p);
                fprintf('%d %s effects are significant p<=%.2f using Hochberg Step Up Procedure for Rroi\n',n_sig,labelfactors{f,1}, p)

                %Multicmp function to apply FDR correction according to Benjamini & Hochberg%%
                [padjFDR_Rroi,alphaFDR_Rroi] = multicmp (pRroi(f,:)','fdr',0.05);
                n_sig = sum(padjFDR_Rroi <= p);
                fprintf('%d %s effects are significant p<=%.2f using FDR correction for Rroi\n',n_sig,labelfactors{f,1}, p)
            end
        end
        
        %%tableaux des résultats%%%%%%%%%%%%%%%%
        if tablegraphmode == 1
            tblresRroi = array2table(resRroi);
            for r = 1:(width(tbldataRroi)- 5)
                tblresRroi = mergevars(tblresRroi,(r:r+6));
            end
            tblresRroi.Properties.VariableNames = labelroiALL{1,2};
            tblresRroi.Properties.RowNames = {'Source','Group','SES','Group*SES','Error','Total'};

            tblpRroi = array2table(pRroi);
            tblpRroi.Properties.VariableNames = labelroiALL{1,2};
            tblpRroi.Properties.RowNames = {'Group','SES','Group*SES'};

            sig = tblpRroi.Variables <= p;
            tblpsigRroi = tblpRroi(:,any(sig));
            tblsigRroi = tblresRroi(:,tblpsigRroi.Properties.VariableNames);

            grsig = tblpRroi{{'Group'},:} <= p;
            tblgrpsigRroi = tblpRroi(1,grsig);
            tblgrsigRroi = tblresRroi(:, tblgrpsigRroi.Properties.VariableNames);

            grsessig = tblpRroi{{'Group*SES'},:} <= p;
            tblgrsespsigRroi = tblpRroi(3,grsessig);
            tblgrsessigRroi = tblresRroi(:, tblgrsespsigRroi.Properties.VariableNames);

            disp(tblgrpsigRroi)
            disp(tblgrsespsigRroi)

            %writetable(tblsigRroi,[savepath date '_sigresultsRroi.xls']); trop gros
            writetable(tblpsigRroi,[savepath date '_psigresultsRroi.xls'],'WriteRowNames',true);
            writetable(tblgrsigRroi,[savepath date '_grsigresultsRroi.xls']);
            save([savepath date '_resultsRroi.mat'],'tblresRroi','tblpRroi','tblpsigRroi','tblsigRroi','tblgrpsigRroi','tblgrsigRroi');

            %%%% graphique%%%%%%%%%%%%%%%%%%%%
            if find(grsig)
                figure
                A = meanroiG1r{1,R}(:,grsig);
                B = meanroiG2r{1,R}(:,grsig);
                C = [A' B'];
                X = categorical(tblgrpsigRroi.Properties.VariableNames);
                p1 = bar(X,C);
                p1(1).FaceColor = 'r';
                p1(2).FaceColor = 'b';
                ylabel('Pearson correlation');
                savefig([savepath date '_SigRroi']);
                exportgraphics(gcf,[savepath 'SigRroi.png'])
            else
            end
            % figure
            % p1 = bar(meanroiG1r{1,R}(:,grsig),'r');
            % hold on
            % p2 = bar(meanroiG2r{1,R}(:,grsig),'b');
            % xlabel(tblgrpsigRroi.Properties.VariableNames);
            % ylabel('Pearson correlation');
            % hold off

            %%%% Connectogramme%%%%%%%%%%%%
%             %meanG2G1Rroi = meanroiG2r{1,R}-meanroiG1r{1,R}; %Calculer matrice G2-G1
%             meanG1G2Rroi = meanroiG1r{1,R}-meanroiG2r{1,R}; %Calculer matrice G1-G2
% 
%             for r = 1:numel(roi{1,R}(2,:))
%                 for rr = (r + 1):numel(roi{1,R}(2,:))
%                     x = [lgndroi{2,R} num2str(r) '-' lgndroi{2,R} num2str(rr)];
%                     tf = strcmp(x, labelroiALL{1,R}(1,:));
%                     idx = find(tf);
%                     MATroimeanG1G2{1,R}(r,rr) = meanG1G2Rroi(1,idx); %Mettre les  données sous forme de matrice
%                 end
%             end
%             MATroimeanG1G2{1,R}(length(MATroimeanG1G2{1,R}),length(MATroimeanG1G2{1,R})) = 0; %Ajouter le dernier 0 de la diagonale
%             MATroimeanG1G2{1,R} = MATroimeanG1G2{1,R} + triu(MATroimeanG1G2{1,R},1)'; %Répliquer la moitié inférieure de la matrice

            for r = 1:numel(roi{1,R}(2,:))
                for rr = (r + 1):numel(roi{1,R}(2,:))
                    x = [lgndroi{2,R} num2str(r) '-' lgndroi{2,R} num2str(rr)];
                    tf = strcmp(x, labelroiALL{1,R}(1,:));
                    idx = find(tf);
                    MATpRroi(r,rr)= pRroi(1,idx); %Mettre les  p sous forme de matrice
                end
            end

            MATpRroi(length(MATpRroi),length(MATpRroi)) = 0; %Ajouter le dernier 0 de la diagonale
            MATpRroi = MATpRroi + triu(MATpRroi,1)'; %Répliquer la moitié inférieure de la matrice
            tf = MATpRroi <=p & MATpRroi > 0; %Identifier les p significatifs
            MATsigG1G2Rroi = MATroimeanG1G2{1,R}; %Reprendre la matrice de connectivité du ROI
            MATsigG1G2Rroi(~tf) = 0; % Retirer les valeurs non sig

            clear r rr x tf idx

            MATneg = MATsigG1G2Rroi < 0;
            MAT = MATsigG1G2Rroi.*(MATneg);%DATA{id}.MAT %loader matrice
            if find(MAT)
            plotLst = roi{1,R}(2,:);
            label =  roi{1,R}(1,:);
            %fileorderconnectogramRroi  = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Region.txt';
            plotconnectogramroi(fileorderconnectogram{1,R},MAT,label,plotLst)
            savefig([savepath date '_NegConnectSigG1G2Rroi']);
            else
            end     

            MATpos = MATsigG1G2Rroi > 0;
            MAT = MATsigG1G2Rroi.*(MATpos);%DATA{id}.MAT %loader matrice
            if find(MAT)
            plotLst = roi{1,R}(2,:);
            label =  roi{1,R}(1,:);
            %fileorderconnectogramRroi  = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Region.txt';
            plotconnectogramroi(fileorderconnectogram{1,R},MAT,label,plotLst)
            savefig([savepath date '_PosConnectSigG1G2Rroi']);
            else
            end  
        end

        clear r x res sig A B C X p1


        %%FROI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        R = 3;
        x = 1;
        resFroi = [];
        for r = 1:(width(tbldataFroi)- 5)
            [pFroi(:,x),res] = anovan(table2array(tbldataFroi(:,r+5)),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
            close hidden;
            resFroi = [resFroi, res];
            x = x + 1;
        end
        % delete(findall(0));

        labelfactors = resFroi(2:end-2,1);
        for f = 1:numel(labelfactors)
            n_sig = sum(pFroi(f,:) <= p);
            fprintf('%d %s effects are significant p<=%.2f without correction for Froi\n',n_sig,labelfactors{f,1}, p)
        end

        %%fdr%%%%%%%%%%%%
        if fdrmode == 1
            for f = 1:numel(labelfactors)
            % FDR according to Storey 2002%%%
            [fdr_Froi,q_Froi] = mafdr(pFroi(f,:)');
            n_sig = sum(q_Froi <= p);
            fprintf('%d %s effects are significant p<=%.2f using FDR correction for Froi\n',n_sig,labelfactors{f,1}, p)

            % FDR according to Benjamini & Hochberg(1995) de Groppe%%%%
            [h_Froi1, crit_p_Froi, adj_ci_cvrg_Froi, adj_p_Froi] = fdr_bh(pFroi(f,:), 0.05, 'pdep', 'no' );
            n_sig = sum(adj_p_Froi <= p);
            fprintf('%d %s effects are significant p<=%.2f using FDR correction for Froi\n',n_sig,labelfactors{f,1}, p)

            % % FDR according to Benjamini & Hochberg(1995) de Gerber%%%%
            % [ind, thres] = FDR(pFroi(f,:));

            %Bonferroni-Holm correction by Groppe%%% Corrected_p are greater than 1...
            %doesn't make sense
            [corrected_p_Froi, h_Froi2] = bonf_holm(pFroi(f,:));
            n_sig = sum(corrected_p_Froi <= p);
            fprintf('%d %s effects are significant p<=%.2f using Bonferroni-Holm correction for Froi\n',n_sig,labelfactors{f,1}, p)

            %Multcmp function to apply Holm Step Down Procedure%%
            [padjBH_Froi,alphaBH_Froi] = multicmp (pFroi(f,:)','down',0.05);
            n_sig = sum(padjBH_Froi <= p);
            fprintf('%d %s effects are significant p<=%.2f using Holm Step Down Procedure for Froi\n',n_sig,labelfactors{f,1}, p)

            %Multcmp function to apply Hochberg's step up procedure%%
            [padjH_Froi,alphaH_Froi] = multicmp (pFroi(f,:)','up',0.05);
            n_sig = sum(padjH_Froi <= p);
            fprintf('%d %s effects are significant p<=%.2f using Hochberg Step Up Procedure for Froi\n',n_sig,labelfactors{f,1}, p)

            %Multicmp function to apply FDR correction according to Benjamini & Hochberg%%
            [padjFDR_Froi,alphaFDR_Froi] = multicmp (pFroi(f,:)','fdr',0.05);
            n_sig = sum(padjFDR_Froi <= p);
            fprintf('%d %s effects are significant p<=%.2f using FDR correction for Froi\n',n_sig,labelfactors{f,1}, p)
            end
        end

        %%tableaux des résultats%%%%%%%%%%%%%%%%
        if tablegraphmode == 1
            tblresFroi = array2table(resFroi);
            for r = 1:(width(tbldataFroi)- 5)
                tblresFroi = mergevars(tblresFroi,(r:r+6));
            end
            tblresFroi.Properties.VariableNames = labelroiALL{1,3};
            tblresFroi.Properties.RowNames = {'Source','Group','SES','Group*SES','Error','Total'};

            tblpFroi = array2table(pFroi);
            tblpFroi.Properties.VariableNames = labelroiALL{1,3};
            tblpFroi.Properties.RowNames = {'Group','SES','Group*SES'};

            sig = tblpFroi.Variables <= p;
            tblpsigFroi = tblpFroi(:,any(sig));
            tblsigFroi = tblresFroi(:,tblpsigFroi.Properties.VariableNames);

            grsig = tblpFroi{{'Group'},:} <= p;
            tblgrpsigFroi = tblpFroi(1,grsig);
            tblgrsigFroi = tblresFroi(:, tblgrpsigFroi.Properties.VariableNames);

            disp(tblgrpsigFroi)
            
            grsessig = tblpFroi{{'Group*SES'},:} <= p;
            tblgrsespsigFroi = tblpFroi(3,grsessig);
            tblgrsessigFroi = tblresFroi(:, tblgrsespsigFroi.Properties.VariableNames);

            disp(tblgrsespsigFroi)
            
            writetable(tblsigFroi,[savepath date '_sigresultsFroi.xls']);
            writetable(tblgrsigFroi,[savepath date '_grsigresultsFroi.xls']);
            save([savepath date '_resultsFroi.mat'],'tblresFroi','tblpFroi','tblpsigFroi','tblsigFroi','tblgrpsigFroi','tblgrsigFroi');

            %%%% graphique%%%%%%%%%%%%%%%%%%%%
            if find(grsig)
                figure
                A = meanroiG1r{1,R}(:,grsig);
                B = meanroiG2r{1,R}(:,grsig);
                C = [A' B'];
                X = categorical(tblgrpsigFroi.Properties.VariableNames);
                p1 = bar(X,C);
                p1(1).FaceColor = 'r';
                p1(2).FaceColor = 'b';
                ylabel('Pearson correlation');
                savefig([savepath date '_SigFroi']);
                exportgraphics(gcf,[savepath 'SigFroi.png'])
            else
            end
            % figure
            % p1 = bar(meanroiG1r{1,R}(:,grsig),'r');
            % hold on
            % p2 = bar(meanroiG2r{1,R}(:,grsig),'b');
            % xlabel(tblgrpsigFroi.Properties.VariableNames);
            % ylabel('Pearson correlation');
            % hold off
            % savefig([savepath date 'SigFroi']);

            %%%% Connectogramme%%%%%%%%%%%%
%             %meanG2G1Froi = meanroiG2r{1,R}-meanroiG1r{1,R}; %Calculer matrice G2-G1
%             meanG1G2Froi = meanroiG1r{1,R}-meanroiG2r{1,R}; %Calculer matrice G1-G2
% 
%             for r = 1:numel(roi{1,R}(2,:))
%                 for rr = (r + 1):numel(roi{1,R}(2,:))
%                     x = [lgndroi{2,R} num2str(r) '-' lgndroi{2,R} num2str(rr)];
%                     tf = strcmp(x, labelroiALL{1,R}(1,:));
%                     idx = find(tf);
%                     MATroimeanG1G2{1,R}(r,rr) = meanG1G2Froi(1,idx); %Mettre les  données sous forme de matrice
%                 end
%             end
%             MATroimeanG1G2{1,R}(length(MATroimeanG1G2{1,R}),length(MATroimeanG1G2{1,R})) = 0; %Ajouter le dernier 0 de la diagonale
%             MATroimeanG1G2{1,R} = MATroimeanG1G2{1,R} + triu(MATroimeanG1G2{1,R},1)'; %Répliquer la moitié inférieure de la matrice

            for r = 1:numel(roi{1,R}(2,:))
                for rr = (r + 1):numel(roi{1,R}(2,:))
                    x = [lgndroi{2,R} num2str(r) '-' lgndroi{2,R} num2str(rr)];
                    tf = strcmp(x, labelroiALL{1,R}(1,:));
                    idx = find(tf);
                    MATpFroi(r,rr)= pFroi(1,idx); %Mettre les  p sous forme de matrice
                end
            end

            MATpFroi(length(MATpFroi),length(MATpFroi)) = 0; %Ajouter le dernier 0 de la diagonale
            MATpFroi = MATpFroi + triu(MATpFroi,1)'; %Répliquer la moitié inférieure de la matrice
            tf = MATpFroi <=p & MATpFroi > 0; %Identifier les p significatifs
            MATsigG1G2Froi = MATroimeanG1G2{1,R}; %Reprendre la matrice de connectivité du ROI
            MATsigG1G2Froi(~tf) = 0; % Retirer les valeurs non sig

            clear r rr x tf idx

            MATneg = MATsigG1G2Froi < 0;
            MAT = MATsigG1G2Froi.*(MATneg);%DATA{id}.MAT %loader matrice
            if find(MAT)
            plotLst = roi{1,R}(2,:);
            label =  roi{1,R}(1,:);
            %fileorderconnectogramFroi  = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Fonction.txt';
            plotconnectogramroi(fileorderconnectogram{1,R},MAT,label,plotLst)
            savefig([savepath date '_NegConnectSigG1G2Froi']);
            else
            end     

            MATpos = MATsigG1G2Froi > 0;
            MAT = MATsigG1G2Froi.*(MATpos);%DATA{id}.MAT %loader matrice
            if find(MAT)
            plotLst = roi{1,R}(2,:);
            label =  roi{1,R}(1,:);
            %fileorderconnectogramFroi  = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Fonction.txt';
            plotconnectogramroi(fileorderconnectogram{1,R},MAT,label,plotLst)
            savefig([savepath date '_PosConnectSigG1G2Froi']);
            else
            end  

        end

        clear r x res sig A B C X p1


        %%AROI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        R = 4;
        x = 1;
        resAroi = [];
        for r = 1:(width(tbldataAroi)- 5)
            [pAroi(:,x),res] = anovan(table2array(tbldataAroi(:,r+5)),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
            close hidden;
            resAroi = [resAroi, res];
            x = x + 1;
        end
        % delete(findall(0));

        labelfactors = resAroi(2:end-2,1);
        for f = 1:numel(labelfactors)
            n_sig = sum(pAroi(f,:) <= p);
            fprintf('%d %s effects are significant p<=%.2f without correction for Aroi\n',n_sig,labelfactors{f,1}, p)
        end

        %%%%fdr%%%%%%%%%%%%%
        if fdrmode == 1
            for f = 1:numel(labelfactors)
            % FDR according to Storey 2002%%%
            [fdr_Aroi,q_Aroi] = mafdr(pAroi(f,:)');
            n_sig = sum(q_Aroi <= p);
            fprintf('%d %s effects are significant p<=%.2f using FDR correction for Aroi\n',n_sig,labelfactors{f,1}, p)

            % FDR according to Benjamini & Hochberg(1995) de Groppe%%%%
            [h_Aroi1, crit_p_Aroi, adj_ci_cvrg_Aroi, adj_p_Aroi] = fdr_bh(pAroi(f,:)', 0.05, 'pdep', 'no' );
            n_sig = sum(adj_p_Aroi <= p);
            fprintf('%d %s effects are significant p<=%.2f using FDR correction for Aroi\n',n_sig,labelfactors{f,1}, p)

            % % FDR according to Benjamini & Hochberg(1995) de Gerber%%%%
            %[ind, thres] = FDR(pAroi(f,:));

            %Bonferroni-Holm correction by Groppe%%% Corrected_p are greater than 1... doesn't make sense
            [corrected_p_Aroi, h_Aroi2] = bonf_holm(pAroi(f,:)', 0.05);
            n_sig = sum(corrected_p_Aroi <= p);
            fprintf('%d %s effects are significant p<=%.2f using Bonferroni-Holm correction for Aroi\n',n_sig,labelfactors{f,1}, p)

            %Multicmp function to apply Holm Step Down Procedure%%
            [padjBH_Aroi,alphaBH_Aroi] = multicmp (pAroi(f,:)','down',0.05);
            n_sig = sum(padjBH_Aroi <= p);
            fprintf('%d %s effects are significant p<=%.2f using Holm Step Down Procedure for Aroi\n',n_sig,labelfactors{f,1}, p)

            %Multicmp function to apply Hochberg's step up procedure%%
            [padjH_Aroi,alphaH_Aroi] = multicmp (pAroi(f,:)','up',0.05);
            n_sig = sum(padjH_Aroi <= p);
            fprintf('%d %s effects are significant p<=%.2f using Hochberg Step Up Procedure for Aroi\n',n_sig,labelfactors{f,1}, p)

            %Multicmp function to apply FDR correction according to Benjamini & Hochberg%%
            [padjFDR_Aroi,alphaFDR_Aroi] = multicmp (pAroi(f,:)','fdr',0.05);
            n_sig = sum(padjFDR_Aroi <= p);
            fprintf('%d %s effects are significant p<=%.2f using FDR correction for Aroi\n',n_sig,labelfactors{f,1}, p)
            end
        end

        %%tableaux des résultats%%%%%%%%%%%%%%%%
        if tablegraphmode == 1
            tblresAroi = array2table(resAroi);
            for r = 1:(width(tbldataAroi)- 5)
                tblresAroi = mergevars(tblresAroi,(r:r+6));
            end
            tblresAroi.Properties.VariableNames = labelroiALL{1,4};
            tblresAroi.Properties.RowNames = {'Source','Group','SES','Group*SES','Error','Total'};

            tblpAroi = array2table(pAroi);
            tblpAroi.Properties.VariableNames = labelroiALL{1,4};
            tblpAroi.Properties.RowNames = {'Group','SES','Group*SES'};

            sig = tblpAroi.Variables <= p;
            tblpsigAroi = tblpAroi(:,any(sig));
            tblsigAroi = tblresAroi(:,tblpsigAroi.Properties.VariableNames);

            grsig = tblpAroi{{'Group'},:} <= p;
            tblgrpsigAroi = tblpAroi(1,grsig);
            tblgrsigAroi = tblresAroi(:, tblgrpsigAroi.Properties.VariableNames);

            disp(tblgrpsigAroi)
            
            grsessig = tblpAroi{{'Group*SES'},:} <= p;
            tblgrsespsigAroi = tblpAroi(3,grsessig);
            tblgrsessigAroi = tblresAroi(:, tblgrsespsigAroi.Properties.VariableNames);

            disp(tblgrsespsigAroi)
            
            writetable(tblsigAroi,[savepath date '_sigresultsAroi.xls']);
            writetable(tblgrsigAroi,[savepath date '_grsigresultsAroi.xls']);
            save([savepath date '_resultsAroi.mat'],'tblresAroi','tblpAroi','tblpsigAroi','tblsigAroi','tblgrpsigAroi','tblgrsigAroi');

            %%%% graphique%%%%%%%%%%%%%%%%%%%%
            if find(grsig)
                figure
                A = meanroiG1r{1,R}(:,grsig);
                B = meanroiG2r{1,R}(:,grsig);
                C = [A' B'];
                X = categorical(tblgrpsigAroi.Properties.VariableNames);
                p1 = bar(X,C);
                p1(1).FaceColor = 'r';
                p1(2).FaceColor = 'b';
                ylabel('Pearson correlation');
                savefig([savepath date '_SigAroi']);
                exportgraphics(gcf,[savepath 'SigAroi.png'])
            else
            end
            % p1 = bar(meanroiG1r{1,R}(:,grsig),'r');
            % hold on
            % p2 = bar(meanroiG2r{1,R}(:,grsig),'b');
            % xlabel(tblgrpsigAroi.Properties.VariableNames);
            % ylabel('Pearson correlation');
            % hold off
            % savefig([savepath date 'SigAroi']);

                %%%% Connectogramme%%%%%%%%%%%%
%             %meanG2G1Aroi = meanroiG2r{1,R}-meanroiG1r{1,R}; %Calculer matrice G2-G1
%             meanG1G2Aroi = meanroiG1r{1,R}-meanroiG2r{1,R}; %Calculer matrice G2-G1
% 
%             for r = 1:numel(roi{1,R}(2,:))
%                 for rr = (r + 1):numel(roi{1,R}(2,:))
%                     x = [lgndroi{2,R} num2str(r) '-' lgndroi{2,R} num2str(rr)];
%                     tf = strcmp(x, labelroiALL{1,R}(1,:));
%                     idx = find(tf);
%                     MATroimeanG1G2{1,R}(r,rr) = meanG1G2Aroi(1,idx); %Mettre les  données sous forme de matrice
%                 end
%             end
%             MATroimeanG1G2{1,R}(length(MATroimeanG1G2{1,R}),length(MATroimeanG1G2{1,R})) = 0; %Ajouter le dernier 0 de la diagonale
%             MATroimeanG1G2{1,R} = MATroimeanG1G2{1,R} + triu(MATroimeanG1G2{1,R},1)'; %Répliquer la moitié inférieure de la matrice

            for r = 1:numel(roi{1,R}(2,:))
                for rr = (r + 1):numel(roi{1,R}(2,:))
                    x = [lgndroi{2,R} num2str(r) '-' lgndroi{2,R} num2str(rr)];
                    tf = strcmp(x, labelroiALL{1,R}(1,:));
                    idx = find(tf);
                    MATpAroi(r,rr)= pAroi(1,idx); %Mettre les  p sous forme de matrice
                end
            end

            MATpAroi(length(MATpAroi),length(MATpAroi)) = 0; %Ajouter le dernier 0 de la diagonale
            MATpAroi = MATpAroi + triu(MATpAroi,1)'; %Répliquer la moitié inférieure de la matrice
            tf = MATpAroi <=p & MATpAroi > 0; %Identifier les p significatifs
            MATsigG1G2Aroi = MATroimeanG1G2{1,R}; %Reprendre la matrice de connectivité du ROI
            MATsigG1G2Aroi(~tf) = 0; % Retirer les valeurs non sig

            clear r rr x tf idx

            MATneg = MATsigG1G2Aroi < 0;
            MAT = MATsigG1G2Aroi.*(MATneg);%DATA{id}.MAT %loader matrice
            if find(MAT)
            plotLst = roi{1,R}(2,:);
            label =  roi{1,R}(1,:);
            %fileorderconnectogramAroi  = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Aire.txt';
            plotconnectogramroi(fileorderconnectogram{1,R},MAT,label,plotLst)
            savefig([savepath date '_NegConnectSigG1G2Aroi']);
            else
            end     

            MATpos = MATsigG1G2Aroi > 0;
            MAT = MATsigG1G2Aroi.*(MATpos);%DATA{id}.MAT %loader matrice
            if find(MAT)
            plotLst = roi{1,R}(2,:);
            label =  roi{1,R}(1,:);
            %fileorderconnectogramAroi  = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Aire.txt';
            plotconnectogramroi(fileorderconnectogram{1,R},MAT,label,plotLst)
            savefig([savepath date '_PosConnectSigG1G2Aroi']);
            else
            end  

        end

        clear r x res sig A B C X p1

        if tablegraphmode == 1
            tblgrpsigALL = [tblgrpsigMroi, tblgrpsigRroi, tblgrpsigFroi, tblgrpsigAroi];
            writetable(tblgrpsigALL,[savepath date '_grpsigresultsALL.xls'],'WriteRowNames',true);
        end
    end
    
 X = ['Results saved in ', savepath];
disp(X)
clear X   
end



toc