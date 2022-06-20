function FDR_CO(pval, labelfactors, labels, p, dat, savepath)

sz = size(pval);
    for f = 1:numel(labelfactors)
        % FDR according to Storey 2002%%% %BAD ONE
%         [fdr,q] = mafdr(pval(f,:)');
%         n_sig = sum(q <= p);
%         fprintf('%d %s effects are significant p<=%.2f using FDR correction (Storey2002) for %s\n',n_sig,labelfactors{f,1}, p, dat)

        % FDR according to Benjamini & Hochberg(1995) de Groppe%%%%
        [h1, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pval(f,:), 0.05, 'pdep', 'no' );
        n_sig = sum(adj_p(1,:) <= p);
        fprintf('%d %s effects are significant p<=%.2f using FDR correction (Benjamini1995) for %s\n',n_sig,labelfactors{f,1}, p, dat)

        if numel(sz) == 2
            if n_sig >= 1
                isig = find(adj_p(1,:) <= p);
                for e = 1:numel(isig)
                    fprintf('\t Significant effect: p=%.3f for %s\n', adj_p(1,isig(e)),['c' num2str(labels{1,isig(e)})])
                end
            end
        elseif numel(sz) == 3
            labels = repmat(labels,1,sz(2));
            if n_sig >= 1
                isig = find(adj_p(1,:) <= p);
                for e = 1:numel(isig)
                    fprintf('\t Significant effect: p=%.3f for %s\n', adj_p(1,isig(e)),['c' num2str(labels{1,isig(e)})])
                end
            end
        end

        % % FDR according to Benjamini & Hochberg(1995) de Gerber%%%%
        % [ind, thres] = FDR(pch(f,:));

%         %Bonferroni-Holm correction by Groppe%%% Corrected_p are greater than 1...
%         %doesn't make sense
%         [corrected_p, h2] = bonf_holm(pval(f,:));
%         n_sig = sum(corrected_p <= p);
%         fprintf('%d %s effects are significant p<=%.2f using Bonferroni-Holm correction for %s\n',n_sig,labelfactors{f,1}, p, dat)
% 
%         %Multcmp function to apply Holm Step Down Procedure%%
%         [padjBH,alphaBH] = multicmp (pval(f,:)','down',0.05);
%         n_sig = sum(padjBH <= p);
%         fprintf('%d %s effects are significant p<=%.2f using Holm Step Down Procedure for %s\n',n_sig,labelfactors{f,1}, p, dat)
% 
%         %Multcmp function to apply Hochberg's step up procedure%%
%         [padjH,alphaH] = multicmp (pval(f,:)','up',0.05);
%         n_sig = sum(padjH <= p);
%         fprintf('%d %s effects are significant p<=%.2f using Hochberg Step Up Procedure for %s\n',n_sig,labelfactors{f,1}, p, dat)
% 
%         %Multicmp function to apply FDR correction according to Benjamini & Hochberg%%
%         [padjFDR,alphaFDR] = multicmp (pval(f,:)','fdr',0.05);
%         n_sig = sum(padjFDR <= p);
%         fprintf('%d %s effects are significant p<=%.2f using FDR correction for %s\n',n_sig,labelfactors{f,1}, p, dat)
    end
end