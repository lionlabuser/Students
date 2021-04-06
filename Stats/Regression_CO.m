%%%%%%%%%%%%%%REGRESSION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
disp('Computing Regression')

data = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Stats\CORR0,01_0,08\';
load ([data 'workspace.mat'])

savepath='C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Stats\CORR0,01_0,08\Regression\';
channelmode = 1;
tablegraphmode = 1;

%%CHANNELS%%%%%%%%%%%%%%%%%%%%%
if channelmode ==1
%     bch = [];
%     rch = [];
%     statsch = [];
%     x2 = [gr ses];
%     x1 = ones(size(x2,1),1);
%     X = [x1 x2];
%     for a = 1:(width(tblch)- 5)
%         [b,~,r,~,stats] = regress(table2array(tblch(:,a+5)),X);
%         bch = [bch, b];
%         rch = [rch, r];
%         statsch = [statsch, stats];
%     end
    
%     pch = statsch(1,3:4:end);
%     n_sig = sum(pch(1,:) <= 0.07);
%     fprintf('%d tests are significant p<=0.07 without correction for channels\n',n_sig)

    x = 1;
    tblresch = [];
    for a = 1:(width(tblch)- 5)
        mdl = fitlm(tblch(:,[2 5 a+5]));
        pgrch(x) = mdl.Coefficients{2,4};
        psesch(x) = mdl.Coefficients{3,4};
        tblresch = [tblresch table(mdl.Coefficients,'VariableNames',tblch.Properties.VariableNames(a+5))];
        x = x + 1;
    end
    
    n_sig = sum(pgrch(1,:) <= 0.07);
    fprintf('%d tests are significant p<=0.07 without correction for channels\n',n_sig)

    %%fdr%%%%%%%%%%%%

    % % FDR according to Storey 2002%%%
    % [fdr,q] = mafdr(pch(1,:)');
    % n_sig = sum(q <= 0.07);
    % fprintf('%d tests are significant p<=0.07 using FDR correction for channels\n',n_sig)

    % FDR according to Benjamini & Hochberg(1995) de Groppe%%%%
    [h_ch1, crit_p_ch, adj_ci_cvrg_ch, adj_p_ch] = fdr_bh(pch(1,:), 0.05, 'pdep', 'no' );
    n_sig = sum(adj_p_ch(1,:) <= 0.07);
    fprintf('%d tests are significant p<=0.07 using FDR correction for channels\n',n_sig)

    % % FDR according to Benjamini & Hochberg(1995) de Gerber%%%%
    % [ind, thres] = FDR(pch(1,:));

    %Bonferroni-Holm correction by Groppe%%% Corrected_p are greater than 1...
    %doesn't make sense
    [corrected_p_ch, h_ch2] = bonf_holm(pch(1,:));
    n_sig = sum(corrected_p_ch <= 0.07);
    fprintf('%d tests are significant p<=0.07 using Bonferroni-Holm correction for channels\n',n_sig)

    %Multcmp function to apply Bonferroni-Holm correction%%
    [padjBH_ch,alphaBH_ch] = multicmp (pch(1,:)','down',0.05);
    n_sig = sum(padjBH_ch <= 0.07);
    fprintf('%d tests are significant p<=0.07 using Bonferroni-Holm correction for channels\n',n_sig)
    
    %Multcmp function to apply Hochberg's step up procedure%%
    [padjH_ch,alphaH_ch] = multicmp (pch(1,:)','up',0.05);
    n_sig = sum(padjH_ch <= 0.07);
    fprintf('%d tests are significant p<=0.07 using Hochberg Step Up Procedure for channels\n',n_sig)
    
    
    %%tableaux des résultats%%%%%%%%%%%%%%%%
    if tablegraphmode == 1

        tblpgrch = array2table(pgrch);
        tblpgrch.Properties.VariableNames = labelch;
        tblpgrch.Properties.RowNames = {'Group'};
        
        grsig = tblpgrch{1,:} <= 0.07;
        tblpgrsigch = tblpgrch(1,grsig);
        
        tblpsesch = array2table(psesch);
        tblpsesch.Properties.VariableNames = labelch;
        tblpsesch.Properties.RowNames = {'SES'};
        
        sessig = tblpsesch{1,:} <= 0.07;
        tblpsessigch = tblpsesch(1,sessig);

        disp(tblpgrsigch)

        pch = [pgrch; psesch];
        tblpch = array2table(pch);
        tblpch.Properties.VariableNames = labelch;
        tblpch.Properties.RowNames = {'Group', 'SES'};
        
        writetable(tblresch,[savepath date '_resultsch.xls']); %trop gros
        writetable(tblpch,[savepath date '_presultsch.xls'],'WriteRowNames',true); %trop gros
        writetable(tblpgrsigch,[savepath date '_grsigresultsch.xls']);
        save([savepath date '_resultsch.mat'],'tblresch','tblpch','tblpgrsigch','tblpsessigch');

        %%%% graphique%%%%%%%%%%%%%%%%%%%%
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
    end
end

clear r x res sig A B C X p1
    
X = ['Results saved in ', savepath];
disp(X)
clear X

toc