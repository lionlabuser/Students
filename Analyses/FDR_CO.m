%%%%fdr%%%%%%%%%%%%%
p = 0,07;

% % FDR according to Storey 2002%%%
[fdr,q] = mafdr(pvalues(1,:)');
n_sig = sum(q <= p);
fprintf('%d tests are significant p<=%.2f using FDR correction (Storey, 2002)\n',n_sig, p)

% FDR according to Benjamini & Hochberg(1995) de Groppe%%%%
[h1, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvalues(1,:)', 0.05, 'pdep', 'no' );
n_sig = sum(adj_p <= p);
fprintf('%d tests are significant p<=%.2f using FDR correction (Benjamini,1995)\n',n_sig, p)

% % FDR according to Benjamini & Hochberg(1995) de Gerber%%%%
%[ind, thres] = FDR(pvalues(1,:));

%Bonferroni-Holm correction by Groppe%%% Corrected_p are greater than 1... doesn't make sense
[corrected_p, h] = bonf_holm(pvalues(1,:)', 0.05);
n_sig = sum(corrected_p <= p);
fprintf('%d tests are significant p<=%.2f using Bonferroni-Holm correction\n',n_sig, p)

%Multicmp function to apply Holm Step Down Procedure%%
[padjBH,alphaBH] = multicmp (pvalues(1,:)','down',0.05);
n_sig = sum(padjBH <= p);
fprintf('%d tests are significant p<=%.2f using Holm Step Down Procedure\n',n_sig, p)
    
%Multicmp function to apply Hochberg's step up procedure%%
[padjH,alphaH] = multicmp (pvalues(1,:)','up',0.05);
n_sig = sum(padjH <= p);
fprintf('%d tests are significant p<=%.2f using Hochberg Step Up Procedure\n',n_sig, p)

%Multicmp function to apply FDR correction according to Benjamini & Hochberg%%
[padjFDR,alphaFDR] = multicmp (pvalues(1,:)','fdr',0.05);
n_sig = sum(padjFDR <= p);
fprintf('%d tests are significant p<=%.2f using FDR correction (Benjamini, 1995)\n',n_sig, p)