%%%%%%%%%%%%%%ANOVA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
disp('Computing T-test')

data = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Stats\CORR0,01_0,08\';
load ([data 'workspace.mat'])

savepath='C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Stats\CORR0,01_0,08\Ttest\';
channelmode = 1;
tablegraphmode = 1;
p = 0.05;

%%CHANNELS%%%%%%%%%%%%%%%%%%%%%
if channelmode ==1
    x = 1;
    for r = 1:(width(tblch)- 5)
        [hch(:,x), pch(:,x)] = ttest2(table2array(tblch(idG1,r+5)),table2array(tblch(idG2,r+5)));
        close hidden;
        x = x + 1;
    end
    % delete(findall(0));

    n_sig = sum(pch(1,:) <= p);
    fprintf('%d tests are significant p<=%.2f without correction for channels\n',n_sig, p)

    clear r x res n_sig
    
    %%fdr%%%%%%%%%%%%

    % % FDR according to Storey 2002%%%
    % [fdr,q] = mafdr(pch(1,:)');
    % n_sig = sum(q <= 0.07);
    % fprintf('%d tests are significant p<=0.07 using FDR correction for channels\n',n_sig)

    % FDR according to Benjamini & Hochberg(1995) de Groppe%%%%
    [h_ch1, crit_p_ch, adj_ci_cvrg_ch, adj_p_ch] = fdr_bh(pch(1,:), 0.05, 'pdep', 'no' );
    n_sig = sum(adj_p_ch(1,:) <= 0.07);
    fprintf('%d tests are significant p<=%.2f using FDR correction for channels\n',n_sig, p)

    % % FDR according to Benjamini & Hochberg(1995) de Gerber%%%%
    % [ind, thres] = FDR(pch(1,:));

    %Bonferroni-Holm correction by Groppe%%% Corrected_p are greater than 1...
    %doesn't make sense
    [corrected_p_ch, h_ch2] = bonf_holm(pch(1,:));
    n_sig = sum(corrected_p_ch <= 0.07);
    fprintf('%d tests are significant p<=%.2f using Bonferroni-Holm correction for channels\n',n_sig, p)

    %Multcmp function to apply Bonferroni-Holm correction%%
    [padjBH_ch,alphaBH_ch] = multicmp (pch(1,:)','down',0.05);
    n_sig = sum(padjBH_ch <= 0.07);
    fprintf('%d tests are significant p<=%.2f using Bonferroni-Holm correction for channels\n',n_sig, p)
    
    %Multcmp function to apply Hochberg's step up procedure%%
    [padjH_ch,alphaH_ch] = multicmp (pch(1,:)','up',0.05);
    n_sig = sum(padjH_ch <= 0.07);
    fprintf('%d tests are significant p<=%.2f using Hochberg Step Up Procedure for channels\n',n_sig, p)
    
    %%tableaux des résultats%%%%%%%%%%%%%%%%
    if tablegraphmode == 1
        tblpch = array2table(pch);
        tblpch.Properties.VariableNames = labelch;
        tblpch.Properties.RowNames = {'Group'};

        sig = tblpch.Variables <= p;
        tblpsigch = tblpch(:,sig);

        disp(tblpsigch)

        %writetable(tblsigch,[savepath date '_sigresultsch.xls']); trop gros
        writetable(tblpsigch,[savepath date '_psig0,05resultsch.xls'],'WriteRowNames',true);
        save([savepath date '_results0,05ch.mat'],'tblpch','tblpsigch');

        %%%% graphique%%%%%%%%%%%%%%%%%%%%
        figure
        A = meanG1ch(:,sig);
        B = meanG2ch(:,sig);
        C = [A' B'];
        X = categorical(tblpsigch.Properties.VariableNames);
        p1 = bar(X,C);
        p1(1).FaceColor = 'r';
        p1(2).FaceColor = 'b';
        ylabel('Pearson correlation');
        savefig([savepath date '_Sig0,05Ch']);
        
        clear A B C X p1 sig grsig
        
        %%%% Connectogramme%%%%%%%%%%%%
        MATmeanG2G1 = MATmeanG2-MATmeanG1; %Calculer matrice G2-G1
        %MATmeanG1G2 = MATmeanG1 - MATmeanG2;
        
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
        MATsigG2G1ch = MATmeanG2G1;
        %MATsigG1G2ch = MATmeanG1G2;
        MATsigG2G1ch(~tf) = 0;
        %MATsigG1G2ch(~tf) = 0;
        
        clear c cc x tf idx
        
        id = 1;
        
        MATneg = MATsigG2G1ch < 0;
        MAT = MATsigG2G1ch.*(MATneg);%DATA{id}.MAT %loader matrice
        if find(MAT)
            List = strvcat(DATA{id}.ZoneList); %liste des paires SD
            ML = DATA{id}.zone.ml; %Loader S/D/ROI/Gr
            plotLst = DATA{id}.zone.plotLst;
            label =  DATA{id}.zone.label;
            fileorderconnectogram  = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Mixte.txt';
            plotconnectogram(fileorderconnectogram,MAT,List,label,plotLst,ML)
            savefig([savepath date '_NegConnectSig0,05G2G1Ch']);
        else
        end
        
        MATpos = MATsigG2G1ch > 0;
        MAT = MATsigG2G1ch.*(MATpos);%DATA{id}.MAT %loader matrice
        if find(MAT)
            List = strvcat(DATA{id}.ZoneList); %liste des paires SD
            ML = DATA{id}.zone.ml; %Loader S/D/ROI/Gr
            plotLst = DATA{id}.zone.plotLst;
            label =  DATA{id}.zone.label;
            fileorderconnectogram  = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Mixte.txt';
            plotconnectogram(fileorderconnectogram,MAT,List,label,plotLst,ML)
            savefig([savepath date '_PosConnectSig0,05G2G1Ch']);
        else
        end
    end
end

clear r x res sig A B C X p1


%%MROI%%%%%%%%%%%%%%%%%%%%%%%
R = 1;
x = 1;
for r = 1:(width(tblmeanMroi)- 5)
    [hMroi(:,x),pMroi(:,x)] = ttest2(table2array(tblmeanMroi(idG1,r+5)),table2array(tblmeanMroi(idG2,r+5)));
    close hidden;
    x = x + 1;
end
% delete(findall(0));

n_sig = sum(pMroi(1,:) <= 0.07);
fprintf('%d tests are significant p<=0.07 without correction for Mroi\n',n_sig)

%%fdr%%%%%%%%%%%%

% % FDR according to Storey 2002%%%
% [fdr,q] = mafdr(pMroi(1,:)');
% n_sig = sum(q <= 0.07);
% fprintf('%d tests are significant p <=0.07 using FDR correction for Mroi\n',n_sig)

% FDR according to Benjamini & Hochberg(1995) de Groppe%%%%
[h_Mroi1, crit_p_Mroi, adj_ci_cvrg_Mroi, adj_p_Mroi] = fdr_bh(pMroi(1,:), 0.05, 'pdep', 'no' );
n_sig = sum(adj_p_Mroi <= 0.07);
fprintf('%d tests are significant p<=0.07 using FDR correction for Mroi\n',n_sig)

% % FDR according to Benjamini & Hochberg(1995) de Gerber%%%%
% [ind, thres] = FDR(pMroi(1,:));

%Bonferroni-Holm correction by Groppe%%% Corrected_p are greater than 1...
%doesn't make sense
[corrected_p_Mroi, h_Mroi2] = bonf_holm(pMroi(1,:));
n_sig = sum(corrected_p_Mroi <= 0.07);
fprintf('%d tests are significant p<=0.07 using Bonferroni-Holm correction for Mroi\n',n_sig)

%Multcmp function to apply Bonferroni-Holm correction%%
[padjBH_Mroi,alphaBH_Mroi] = multicmp (pMroi(1,:)','down',0.05);
n_sig = sum(padjBH_Mroi <= 0.07);
fprintf('%d tests are significant p<=0.07 using Bonferroni-Holm correction for Mroi\n',n_sig)
    
%Multcmp function to apply Hochberg's step up procedure%%
[padjH_Mroi,alphaH_Mroi] = multicmp (pMroi(1,:)','up',0.05);
n_sig = sum(padjH_Mroi <= 0.07);
fprintf('%d tests are significant p<=0.07 using Hochberg Step Up Procedure for Mroi\n',n_sig)

%%tableaux des résultats%%%%%%%%%%%%%%%%
if tablegraphmode == 1
    tblpMroi = array2table(pMroi);
    tblpMroi.Properties.VariableNames = labelroi{1,1};
    tblpMroi.Properties.RowNames = {'Group'};

    sig = tblpMroi.Variables <= 0.07;
    tblpsigMroi = tblpMroi(:,sig);

    disp(tblpsigMroi)
    writetable(tblpsigMroi,[savepath date '_psigresultsMroi.xls']);
    save([savepath date '_resultsMroi.mat'],'tblpMroi','tblpsigMroi');

    %%%% graphique%%%%%%%%%%%%%%%%%%%%
    figure
    A = meanG1Mroi(:,sig);
    B = meanG2Mroi(:,sig);
    C = [A' B'];
    X = categorical(tblpsigMroi.Properties.VariableNames);
    p1 = bar(X,C);
    p1(1).FaceColor = 'r';
    p1(2).FaceColor = 'b';
    ylabel('Pearson correlation');
    savefig([savepath date '_SigMroi']);

        clear A B C X p1 sig grsig
        
    %%%% Connectogramme%%%%%%%%%%%%
    meanG2G1Mroi = meanG2Mroi-meanG1Mroi; %Calculer matrice G2-G1
    %meanG1G2Mroi = meanG1Mroi-meanG2Mroi; %Calculer matrice G1-G2
    
    for r = 1:numel(roi{1,R}(2,:))
        for rr = (r + 1):numel(roi{1,R}(2,:))
            x = [lgndroi{2,R} num2str(r) '-' lgndroi{2,R} num2str(rr)];
            tf = strcmp(x, labelroi{1,R}(1,:));
            idx = find(tf);
            MATmeanG2G1Mroi(r,rr) = meanG2G1Mroi(1,idx); %Mettre les  p sous forme de matrice
        end
    end
    MATmeanG2G1Mroi(length(MATmeanG2G1Mroi),length(MATmeanG2G1Mroi)) = 0; %Ajouter le dernier 0 de la diagonale
    MATmeanG2G1Mroi = MATmeanG2G1Mroi + triu(MATmeanG2G1Mroi,1)'; %Répliquer la moitié inférieure de la matrice
    
    for r = 1:numel(roi{1,R}(2,:))
        for rr = (r + 1):numel(roi{1,R}(2,:))
            x = [lgndroi{2,R} num2str(r) '-' lgndroi{2,R} num2str(rr)];
            tf = strcmp(x, labelroi{1,R}(1,:));
            idx = find(tf);
            MATpMroi(r,rr)= pMroi(1,idx); %Mettre les  p sous forme de matrice
        end
    end
    
    MATpMroi(length(MATpMroi),length(MATpMroi)) = 0; %Ajouter le dernier 0 de la diagonale
    MATpMroi = MATpMroi + triu(MATpMroi,1)'; %Répliquer la moitié inférieure de la matrice
    tf = MATpMroi <=p & MATpMroi > 0; %Identifier les p significatifs
    MATsigG2G1Mroi = MATmeanG2G1Mroi; %Reprendre la matrice de connectivité du ROI
    MATsigG2G1Mroi(~tf) = 0; % Retirer les valeurs non sig
    
    clear r rr x tf idx
    
    MATneg = MATsigG2G1Mroi < 0;
    MAT = MATsigG2G1Mroi.*(MATneg);%DATA{id}.MAT %loader matrice
    if find(MAT)
    plotLst = roi{1,R}(2,:);
    label =  roi{1,R}(1,:);
    fileorderconnectogram  = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Mixte.txt';
    plotconnectogramroi(fileorderconnectogram,MAT,label,plotLst)
    savefig([savepath date '_NegConnectSigG2G1Mroi']);
    else
    end     
    
    MATpos = MATsigG2G1Mroi > 0;
    MAT = MATsigG2G1Mroi.*(MATpos);%DATA{id}.MAT %loader matrice
    if find(MAT)
    plotLst = roi{1,R}(2,:);
    label =  roi{1,R}(1,:);
    fileorderconnectogram  = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Mixte.txt';
    plotconnectogramroi(fileorderconnectogram,MAT,label,plotLst)
    savefig([savepath date '_PosConnectSigG2G1Mroi']);
    else
    end    
end

clear r x res sig A B C X p1


%%RROI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = 2;
x = 1;
for r = 1:(width(tblmeanRroi)- 5)
    [hRroi(:,x),pRroi(:,x)] = ttest2(table2array(tblmeanRroi(idG1,r+5)),table2array(tblmeanRroi(idG2,r+5)));
    close hidden;
    x = x + 1;
end
% delete(findall(0));

n_sig = sum(pRroi(1,:) <= 0.07);
fprintf('%d tests are significant p<=0.07 without correction for Rroi\n',n_sig)

%%fdr%%%%%%%%%%%%

% % FDR according to Storey 2002%%%
% [fdr,q] = mafdr(pRroi(1,:)');
% n_sig = sum(q <= 0.07);
% fprintf('%d tests are significant p<=0.07 using FDR correction for Rroi\n',n_sig)

% FDR according to Benjamini & Hochberg(1995) de Groppe%%%%
[h_Rroi1, crit_p_Rroi, adj_ci_cvrg_Rroi, adj_p_Rroi] = fdr_bh(pRroi(1,:), 0.05, 'pdep', 'no' );
n_sig = sum(adj_p_Rroi <= 0.07);
fprintf('%d tests are significant p<=0.07 using FDR correction for Rroi\n',n_sig)

% % FDR according to Benjamini & Hochberg(1995) de Gerber%%%%
% [ind, thres] = FDR(pRroi(1,:));

%Bonferroni-Holm correction by Groppe%%% Corrected_p are greater than 1...
%doesn't make sense
[corrected_p_Rroi, h_Rroi2] = bonf_holm(pRroi(1,:));
n_sig = sum(corrected_p_Rroi <= 0.07);
fprintf('%d tests are significant p<=0.07 using Bonferroni-Holm correction for Rroi\n',n_sig)

%Multcmp function to apply Bonferroni-Holm correction%%
[padjBH_Rroi,alphaBH_Rroi] = multicmp (pRroi(1,:)','down',0.05);
n_sig = sum(padjBH_Rroi <= 0.07);
fprintf('%d tests are significant p<=0.07 using Bonferroni-Holm correction for Rroi\n',n_sig)
    
%Multcmp function to apply Hochberg's step up procedure%%
[padjH_Rroi,alphaH_Rroi] = multicmp (pRroi(1,:)','up',0.05);
n_sig = sum(padjH_Rroi <= 0.07);
fprintf('%d tests are significant p<=0.07 using Hochberg Step Up Procedure for Rroi\n',n_sig)

%%tableaux des résultats%%%%%%%%%%%%%%%%
if tablegraphmode == 1
    tblpRroi = array2table(pRroi);
    tblpRroi.Properties.VariableNames = labelroi{1,2};
    tblpRroi.Properties.RowNames = {'Group'};

    sig = tblpRroi.Variables <= 0.07;
    tblpsigRroi = tblpRroi(:,sig);

    disp(tblpsigRroi)
    %writetable(tblsigRroi,[savepath date '_sigresultsRroi.xls']); trop gros
    writetable(tblpsigRroi,[savepath date '_psigresultsRroi.xls'],'WriteRowNames',true);
    save([savepath date '_resultsRroi.mat'],'tblpRroi','tblpsigRroi');

    %%%% graphique%%%%%%%%%%%%%%%%%%%%
    figure
    A = meanG1Rroi(:,sig);
    B = meanG2Rroi(:,sig);
    C = [A' B'];
    X = categorical(tblpsigRroi.Properties.VariableNames);
    p1 = bar(X,C);
    p1(1).FaceColor = 'r';
    p1(2).FaceColor = 'b';
    ylabel('Pearson correlation');
    savefig([savepath date '_SigRroi']);

     %%%% Connectogramme%%%%%%%%%%%%
    meanG2G1Rroi = meanG2Rroi-meanG1Rroi; %Calculer matrice G2-G1
    %meanG1G2Rroi = meanG1Rroi-meanG2Rroi; %Calculer matrice G1-G2
    
    for r = 1:numel(roi{1,R}(2,:))
        for rr = (r + 1):numel(roi{1,R}(2,:))
            x = [lgndroi{2,R} num2str(r) '-' lgndroi{2,R} num2str(rr)];
            tf = strcmp(x, labelroi{1,R}(1,:));
            idx = find(tf);
            MATmeanG2G1Rroi(r,rr) = meanG2G1Rroi(1,idx); %Mettre les  p sous forme de matrice
        end
    end
    MATmeanG2G1Rroi(length(MATmeanG2G1Rroi),length(MATmeanG2G1Rroi)) = 0; %Ajouter le dernier 0 de la diagonale
    MATmeanG2G1Rroi = MATmeanG2G1Rroi + triu(MATmeanG2G1Rroi,1)'; %Répliquer la moitié inférieure de la matrice
    
    for r = 1:numel(roi{1,R}(2,:))
        for rr = (r + 1):numel(roi{1,R}(2,:))
            x = [lgndroi{2,R} num2str(r) '-' lgndroi{2,R} num2str(rr)];
            tf = strcmp(x, labelroi{1,R}(1,:));
            idx = find(tf);
            MATpRroi(r,rr)= pRroi(1,idx); %Mettre les  p sous forme de matrice
        end
    end
    
    MATpRroi(length(MATpRroi),length(MATpRroi)) = 0; %Ajouter le dernier 0 de la diagonale
    MATpRroi = MATpRroi + triu(MATpRroi,1)'; %Répliquer la moitié inférieure de la matrice
    tf = MATpRroi <=p & MATpRroi > 0; %Identifier les p significatifs
    MATsigG2G1Rroi = MATmeanG2G1Rroi; %Reprendre la matrice de connectivité du ROI
    MATsigG2G1Rroi(~tf) = 0; % Retirer les valeurs non sig
    
    clear r rr x tf idx
    
    MATneg = MATsigG2G1Rroi < 0;
    MAT = MATsigG2G1Rroi.*(MATneg);%DATA{id}.MAT %loader matrice
    if find(MAT)
    plotLst = roi{1,R}(2,:);
    label =  roi{1,R}(1,:);
    fileorderconnectogram  = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Region.txt';
    plotconnectogramroi(fileorderconnectogram,MAT,label,plotLst)
    savefig([savepath date '_NegConnectSigG2G1Rroi']);
    else
    end     
    
    MATpos = MATsigG2G1Rroi > 0;
    MAT = MATsigG2G1Rroi.*(MATpos);%DATA{id}.MAT %loader matrice
    if find(MAT)
    plotLst = roi{1,R}(2,:);
    label =  roi{1,R}(1,:);
    fileorderconnectogram  = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Region.txt';
    plotconnectogramroi(fileorderconnectogram,MAT,label,plotLst)
    savefig([savepath date '_PosConnectSigG2G1Rroi']);
    else
    end  
end

clear r x res sig A B C X p1


%%FROI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = 3;
x = 1;
for r = 1:(width(tblmeanFroi)- 5)
    [hFroi(:,x), pFroi(:,x)] = ttest2(table2array(tblmeanFroi(idG1,r+5)),table2array(tblmeanFroi(idG2,r+5)));
    close hidden;
    x = x + 1;
end
% delete(findall(0));

n_sig = sum(pFroi(1,:) <= 0.07);
fprintf('%d tests are significant p<=0.07 without correction for Froi\n',n_sig)

%%fdr%%%%%%%%%%%%

% FDR according to Storey 2002%%%
% [fdr,q] = mafdr(pFroi(1,:)');
% n_sig = sum(q <= 0.07);
% fprintf('%d tests are significant p<=0.07 using FDR correction for Froi\n',n_sig)

% FDR according to Benjamini & Hochberg(1995) de Groppe%%%%
[h_Froi1, crit_p_Froi, adj_ci_cvrg_Froi, adj_p_Froi] = fdr_bh(pFroi(1,:), 0.05, 'pdep', 'no' );
n_sig = sum(adj_p_Froi <= 0.07);
fprintf('%d tests are significant p<=0.07 using FDR correction for Froi\n',n_sig)

% % FDR according to Benjamini & Hochberg(1995) de Gerber%%%%
% [ind, thres] = FDR(pFroi(1,:));

%Bonferroni-Holm correction by Groppe%%% Corrected_p are greater than 1...
%doesn't make sense
[corrected_p_Froi, h_Froi2] = bonf_holm(pFroi(1,:));
n_sig = sum(corrected_p_Froi <= 0.07);
fprintf('%d tests are significant p<=0.07 using Bonferroni-Holm correction for Froi\n',n_sig)

%Multcmp function to apply Bonferroni-Holm correction%%
[padjBH_Froi,alphaBH_Froi] = multicmp (pFroi(1,:)','down',0.05);
n_sig = sum(padjBH_Froi <= 0.07);
fprintf('%d tests are significant p<=0.07 using Bonferroni-Holm correction for Froi\n',n_sig)
    
%Multcmp function to apply Hochberg's step up procedure%%
[padjH_Froi,alphaH_Froi] = multicmp (pFroi(1,:)','up',0.05);
n_sig = sum(padjH_Froi <= 0.07);
fprintf('%d tests are significant p<=0.07 using Hochberg Step Up Procedure for Froi\n',n_sig)

%%tableaux des résultats%%%%%%%%%%%%%%%%
if tablegraphmode == 1
    tblpFroi = array2table(pFroi);
    tblpFroi.Properties.VariableNames = labelroi{1,3};
    tblpFroi.Properties.RowNames = {'Group'};

    sig = tblpFroi.Variables <= 0.07;
    tblpsigFroi = tblpFroi(:,sig);
    
    disp(tblpsigFroi)
    writetable(tblpsigFroi,[savepath date '_psigresultsFroi.xls']);
    save([savepath date '_resultsFroi.mat'],'tblpFroi','tblpsigFroi');

    %%%% graphique%%%%%%%%%%%%%%%%%%%%
    figure
    A = meanG1Froi(:,sig);
    B = meanG2Froi(:,sig);
    C = [A' B'];
    X = categorical(tblpsigFroi.Properties.VariableNames);
    p1 = bar(X,C);
    p1(1).FaceColor = 'r';
    p1(2).FaceColor = 'b';
    ylabel('Pearson correlation');
    savefig([savepath date '_SigFroi']);

   %%%% Connectogramme%%%%%%%%%%%%
    meanG2G1Froi = meanG2Froi-meanG1Froi; %Calculer matrice G2-G1
    %meanG1G2Froi = meanG1Froi-meanG2Froi; %Calculer matrice G1-G2
    
    for r = 1:numel(roi{1,R}(2,:))
        for rr = (r + 1):numel(roi{1,R}(2,:))
            x = [lgndroi{2,R} num2str(r) '-' lgndroi{2,R} num2str(rr)];
            tf = strcmp(x, labelroi{1,R}(1,:));
            idx = find(tf);
            MATmeanG2G1Froi(r,rr) = meanG2G1Froi(1,idx); %Mettre les  p sous forme de matrice
        end
    end
    MATmeanG2G1Froi(length(MATmeanG2G1Froi),length(MATmeanG2G1Froi)) = 0; %Ajouter le dernier 0 de la diagonale
    MATmeanG2G1Froi = MATmeanG2G1Froi + triu(MATmeanG2G1Froi,1)'; %Répliquer la moitié inférieure de la matrice
    
    for r = 1:numel(roi{1,R}(2,:))
        for rr = (r + 1):numel(roi{1,R}(2,:))
            x = [lgndroi{2,R} num2str(r) '-' lgndroi{2,R} num2str(rr)];
            tf = strcmp(x, labelroi{1,R}(1,:));
            idx = find(tf);
            MATpFroi(r,rr)= pFroi(1,idx); %Mettre les  p sous forme de matrice
        end
    end
    
    MATpFroi(length(MATpFroi),length(MATpFroi)) = 0; %Ajouter le dernier 0 de la diagonale
    MATpFroi = MATpFroi + triu(MATpFroi,1)'; %Répliquer la moitié inférieure de la matrice
    tf = MATpFroi <=p & MATpFroi > 0; %Identifier les p significatifs
    MATsigG2G1Froi = MATmeanG2G1Froi; %Reprendre la matrice de connectivité du ROI
    MATsigG2G1Froi(~tf) = 0; % Retirer les valeurs non sig
    
    clear r rr x tf idx
    
    MATneg = MATsigG2G1Froi < 0;
    MAT = MATsigG2G1Froi.*(MATneg);%DATA{id}.MAT %loader matrice
    if find(MAT)
    plotLst = roi{1,R}(2,:);
    label =  roi{1,R}(1,:);
    fileorderconnectogram  = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Fonction.txt';
    plotconnectogramroi(fileorderconnectogram,MAT,label,plotLst)
    savefig([savepath date '_NegConnectSigG2G1Froi']);
    else
    end     
    
    MATpos = MATsigG2G1Froi > 0;
    MAT = MATsigG2G1Froi.*(MATpos);%DATA{id}.MAT %loader matrice
    if find(MAT)
    plotLst = roi{1,R}(2,:);
    label =  roi{1,R}(1,:);
    fileorderconnectogram  = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Fonction.txt';
    plotconnectogramroi(fileorderconnectogram,MAT,label,plotLst)
    savefig([savepath date '_PosConnectSigG2G1Froi']);
    else
    end  
end

clear r x res sig A B C X p1


%%AROI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = 4;
x = 1;
for r = 1:(width(tblmeanAroi)- 5)
    [hAroi(:,x),pAroi(:,x)] = ttest2(table2array(tblmeanAroi(idG1,r+5)),table2array(tblmeanAroi(idG2,r+5)));
    close hidden;
    x = x + 1;
end
% delete(findall(0));

n_sig = sum(pAroi(1,:) <= 0.07);
fprintf('%d tests are significant p<=0.07 without correction for Aroi\n',n_sig)

%%%%fdr%%%%%%%%%%%%%

% % FDR according to Storey 2002%%%
[fdr,q] = mafdr(pAroi(1,:)');
n_sig = sum(q <= 0.07);
fprintf('%d tests are significant p<=0.07 using FDR correction for Aroi\n',n_sig)

% FDR according to Benjamini & Hochberg(1995) de Groppe%%%%
[h_Aroi1, crit_p_Aroi, adj_ci_cvrg_Aroi, adj_p_Aroi] = fdr_bh(pAroi(1,:)', 0.05, 'pdep', 'no' );
n_sig = sum(adj_p_Aroi <= 0.07);
fprintf('%d tests are significant p<=0.07 using FDR correction for Aroi\n',n_sig)

% % FDR according to Benjamini & Hochberg(1995) de Gerber%%%%
[ind, thres] = FDR(pAroi(1,:));

%Bonferroni-Holm correction by Groppe%%% Corrected_p are greater than 1... doesn't make sense
%[corrected_p_Aroi, h_Aroi2] = bonf_holm(pAroi(1,:)', 0.05);
%n_sig = sum(corrected_p_Aroi <= 0.07);
%fprintf('%d tests are significant p<=0.07 using Bonferroni-Holm correction for Aroi\n',n_sig)

%Multicmp function to apply Bonferroni-Holm correction%%
[padjBH_Aroi,alphaBH_Aroi] = multicmp (pAroi(1,:)','down',0.05);
n_sig = sum(padjBH_Aroi <= 0.07);
fprintf('%d tests are significant p<=0.07 using Bonferroni-Holm correction for Aroi\n',n_sig)
    
%Multicmp function to apply Hochberg's step up procedure%%
[padjH_Aroi,alphaH_Aroi] = multicmp (pAroi(1,:)','up',0.05);
n_sig = sum(padjH_Aroi <= 0.07);
fprintf('%d tests are significant p<=0.07 using Hochberg Step Up Procedure for Aroi\n',n_sig)

%Multicmp function to apply FDR correction according to Benjamini & Hochberg%%
[padjFDR_Aroi,alphaFDR_Aroi] = multicmp (pAroi(1,:)','fdr',0.05);
n_sig = sum(padjFDR_Aroi <= 0.07);
fprintf('%d tests are significant p<=0.07 using FDR correction for Aroi\n',n_sig)

%%tableaux des résultats%%%%%%%%%%%%%%%%
if tablegraphmode == 1
    tblpAroi = array2table(pAroi);
    tblpAroi.Properties.VariableNames = labelroi{1,4};
    tblpAroi.Properties.RowNames = {'Group'};

    sig = tblpAroi.Variables <= 0.07;
    tblpsigAroi = tblpAroi(:,sig);

    disp(tblpsigAroi)
    writetable(tblpsigAroi,[savepath date '_psigresultsAroi.xls']);
    save([savepath date '_resultsAroi.mat'],'tblpAroi','tblpsigAroi');

    %%%% graphique%%%%%%%%%%%%%%%%%%%%
    figure
    A = meanG1Aroi(:,sig);
    B = meanG2Aroi(:,sig);
    C = [A' B'];
    X = categorical(tblpsigAroi.Properties.VariableNames);
    p1 = bar(X,C);
    p1(1).FaceColor = 'r';
    p1(2).FaceColor = 'b';
    ylabel('Pearson correlation');
    savefig([savepath date '_SigAroi']);

        %%%% Connectogramme%%%%%%%%%%%%
    meanG2G1Aroi = meanG2Aroi-meanG1Aroi; %Calculer matrice G2-G1
    %meanG1G2Aroi = meanG1Aroi-meanG2Aroi; %Calculer matrice G2-G1
            
    for r = 1:numel(roi{1,R}(2,:))
        for rr = (r + 1):numel(roi{1,R}(2,:))
            x = [lgndroi{2,R} num2str(r) '-' lgndroi{2,R} num2str(rr)];
            tf = strcmp(x, labelroi{1,R}(1,:));
            idx = find(tf);
            MATmeanG2G1Aroi(r,rr) = meanG2G1Aroi(1,idx); %Mettre les  p sous forme de matrice
        end
    end
    MATmeanG2G1Aroi(length(MATmeanG2G1Aroi),length(MATmeanG2G1Aroi)) = 0; %Ajouter le dernier 0 de la diagonale
    MATmeanG2G1Aroi = MATmeanG2G1Aroi + triu(MATmeanG2G1Aroi,1)'; %Répliquer la moitié inférieure de la matrice
    
    for r = 1:numel(roi{1,R}(2,:))
        for rr = (r + 1):numel(roi{1,R}(2,:))
            x = [lgndroi{2,R} num2str(r) '-' lgndroi{2,R} num2str(rr)];
            tf = strcmp(x, labelroi{1,R}(1,:));
            idx = find(tf);
            MATpAroi(r,rr)= pAroi(1,idx); %Mettre les  p sous forme de matrice
        end
    end
    
    MATpAroi(length(MATpAroi),length(MATpAroi)) = 0; %Ajouter le dernier 0 de la diagonale
    MATpAroi = MATpAroi + triu(MATpAroi,1)'; %Répliquer la moitié inférieure de la matrice
    tf = MATpAroi <=p & MATpAroi > 0; %Identifier les p significatifs
    MATsigG2G1Aroi = MATmeanG2G1Aroi; %Reprendre la matrice de connectivité du ROI
    MATsigG2G1Aroi(~tf) = 0; % Retirer les valeurs non sig
    
    clear r rr x tf idx
    
    MATneg = MATsigG2G1Aroi < 0;
    MAT = MATsigG2G1Aroi.*(MATneg);%DATA{id}.MAT %loader matrice
    if find(MAT)
    plotLst = roi{1,R}(2,:);
    label =  roi{1,R}(1,:);
    fileorderconnectogram  = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Aire.txt';
    plotconnectogramroi(fileorderconnectogram,MAT,label,plotLst)
    savefig([savepath date '_NegConnectSigG2G1Aroi']);
    else
    end     
    
    MATpos = MATsigG2G1Aroi > 0;
    MAT = MATsigG2G1Aroi.*(MATpos);%DATA{id}.MAT %loader matrice
    if find(MAT)
    plotLst = roi{1,R}(2,:);
    label =  roi{1,R}(1,:);
    fileorderconnectogram  = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Aire.txt';
    plotconnectogramroi(fileorderconnectogram,MAT,label,plotLst)
    savefig([savepath date '_PosConnectSigG2G1Aroi']);
    else
    end  
end

clear r x res sig A B C X p1

if tablegraphmode == 1
    tblpsigALL = [tblpsigMroi, tblpsigRroi, tblpsigFroi, tblpsigAroi];
    writetable(tblpsigALL,[savepath date '_psigresultsALL.xls'],'WriteRowNames',true);
end
    
X = ['Results saved in ', savepath];
disp(X)
clear X

toc