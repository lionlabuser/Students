%%%%%%%%%%%%%%%%%%%%%%%%%PERMUTATIONS/UNPAIRED TTEST%%%%%%%%%%%%%%%%%%
clear;clc;
tic

%%Permutation of Pedro's team from StatMatrices of LIONIRS toolbox%%%
%BE SURE TO HAVE THE LASTEST VERSION OF THE TOOLBOX 
%PERMUTATOOLBOX

%% change values here according to your project (changer apres le "=")

datapath = 'C:\Users\laura\OneDrive - Universite de Montreal\Downloads\'; %s'assurer de garder le dernier \
savepath='C:\Users\laura\OneDrive - Universite de Montreal\Downloads\'; %s'assurer de garder le dernier \
workspacename='workspace.mat'; %with MATall, idG1 idG2 stored

FDRmode = 1; %If you want diverse FDR correction on your results

%tablegraphmode = 1; %valider avec Kass options
%otherpermmode = 0; %valider avec Kass options

p = 0.05;
pmarginal= 0.1;

perm_parameters.ESTAD = 4;              %  ESTAD - Statistic to use
%         1 - t-student with "original" values
%         2 - sum of difference with "original" values
%         4 - t-student with absolute values
%         5 - sum of difference with absolute values
perm_parameters.INDEP = 1;              %  INDEP- Indice if groups are independent
%         0 - DEPENDENT
%         1 - INDEPENDENT
perm_parameters.NPERM = 1000;


%% permutation analysis
disp('Computing Permutations, Pedro et al. ...')
load ([datapath workspacename]);
MATallr = permute(MATall, [3 2 1]); %reorder so that participants are the first dimensions
[n1,n2,n3]=size(MATallr);
%reorganize matrice into 2D
x = 1;
for c = 1:n2
    for cc = (c + 1):n3
        labelch{x} = [num2str(c) '-' num2str(cc)];
        MATall2D(:,x) = MATallr(:,c,cc);
        x = x + 1;
    end
end

%create a matrice for each group
g1MAT = MATall2D(idG1',:);
g2MAT = MATall2D(idG2',:);

%Permutations%%%
[perm_pGlobal,perm_pChan,perm_Ptime,perm_pALL,perm_tval] = ...
    TestPermut2Grupos(perm_parameters.ESTAD,perm_parameters.INDEP,g1MAT,g2MAT,perm_parameters.NPERM);
perm_sig=perm_pALL;
perm_sig(perm_sig>p)=0;
perm_sigmarg=perm_pALL;
perm_sigmarg(perm_sigmarg>pmarginal)=0; perm_sigmarg(perm_sigmarg<p)=0;
toc

perm_Nsig = sum(perm_pALL <= p);
perm_Nsigmarg = sum(perm_pALL <= pmarginal) - perm_Nsig;
fprintf('\n************************\nSIGNIFICANT CHANNEL PAIRS p<=%.2f:\n',p)
fprintf('%d tests out of %d are significant p<=%.2f without correction\n',perm_Nsig,numel(perm_pALL), p)


y=1;
for x=find(perm_sig>0)
    fprintf('Channel pair %s, p-value: %.3f\n',labelch{x},perm_sig(x))
    perm_summary(y).ChanPair=labelch{x};
    perm_summary(y).pvalue=perm_sig(x);
    y=y+1;
end

fprintf('\n************************\nMARGINALLY SIGNIFICANT CHANNEL PAIRS p<=%.2f:\n',pmarginal)
fprintf('%d tests out of %d are marginally significant p<=%.2f without correction\n',perm_Nsigmarg,numel(perm_pALL), pmarginal)
for x=find(perm_sigmarg>0)
    fprintf('Channel pair %s, p-value: %.3f\n',labelch{x},perm_sigmarg(x))
    perm_summary(y).ChanPair=labelch{x};
    perm_summary(y).pvalue=perm_sigmarg(x);
    y=y+1;
end

x=1;
for c = 1:n2
    for cc = (c + 1):n3
        %STORE in matrice format                  % replicate data 2nd half
        FUniv(c,cc)=perm_pALL(x);                 FUniv(cc,c)=FUniv(c,cc);
        Toij(c,cc)=perm_tval(x);                  Toij(cc,c)=Toij(c,cc);
        ncb1(c,cc)=sum(~isnan(g1MAT(:,x)));       ncb1(cc,c)=ncb1(c,cc);
        npb1(c,cc)=sum(~isnan(g2MAT(:,x)));       npb1(cc,c)=npb1(c,cc);
        x = x + 1;
    end
end
helplabels={'FUniv' 'p value';'Toij' 't-student value' ; 'ncb1' 'number of valid data for G1'; 'npb1' 'number of valid data for G2'};

save([savepath date '_resultsPerm.mat'],'perm*','MATallr','g*MAT','datapath', 'workspacename',...
    'ncb1','npb1','FUniv','Toij','helplabels');
fprintf('Permutation save in:\n%s\n\n______________\n*Script completed*\n',[savepath date '_resultsPerm.mat'])




%% FDR correction mode
if FDRmode == 1
    pch=perm_pALL;
    perm_FDR(1).method='Uncorrected permutation';
    perm_FDR(1).n_sig=perm_Nsig;
    perm_FDR(1).p_adj=pch;
    
    % % FDR according to Storey 2002%%%
    try
        [~,q] = mafdr(pch');
        n_sig = sum(q <= p);
        fprintf('%d tests are significant p<=%.2f using FDR correction (Storey, 2002)\n',n_sig, p)
        perm_FDR(2).method='FDR according to Storey 2002';
        perm_FDR(2).n_sig=n_sig;
        perm_FDR(2).p_adj=q';
    catch
        fprintf('Need to download *mafdr.m* function - FDR according to Storey 2002\n')
    end
    
    % FDR according to Benjamini & Hochberg(1995) de Groppe%%%%
    try
        [~, ~, ~, q] = fdr_bh(pch', p, 'pdep', 'no' );
        n_sig = sum(q <= p);
        fprintf('%d tests are significant p<=%.2f using FDR correction (Benjamini,1995)\n',n_sig, p)
        perm_FDR(3).method='FDR according to Benjamini & Hochberg(1995) by Groppe';
        perm_FDR(3).n_sig=n_sig;
        perm_FDR(3).p_adj=q';
    catch
        fprintf('Need to download *fdr_bh.m* function - FDR according to Benjamini & Hochberg(1995) by Groppe\n')
    end
    
    %Bonferroni-Holm correction by Groppe%%% Corrected_p are greater than 1... doesn't make sense
    try
        [corrected_p, h] = bonf_holm(pch', p);
        n_sig = sum(corrected_p <= p);
        fprintf('%d tests are significant p<=%.2f using Bonferroni-Holm correction\n',n_sig, p)
        perm_FDR(4).method='BBonferroni-Holm correction by Groppe';
        perm_FDR(4).n_sig=n_sig;
        perm_FDR(4).p_adj=q';
    catch
        fprintf('Need to download *bonf_holm.m* function - Bonferroni-Holm correction by Groppe\n')
    end
    
    %Multicmp function to apply Holm Step Down Procedure%%
    try  [q,~] = multicmp (pch','down',p);
        n_sig = sum(q <= p);
        fprintf('%d tests are significant p<=%.2f using Holm Step Down Procedure\n',n_sig, p)
        perm_FDR(5).method='Multicmp function to apply Holm Step Down Procedure';
        perm_FDR(5).n_sig=n_sig;
        perm_FDR(5).p_adj=q';
        
        %Multicmp function to apply Hochberg's step up procedure%%
        
        [q,~] = multicmp (pch','up',p);
        n_sig = sum(q <= p);
        fprintf('%d tests are significant p<=%.2f using Hochberg Step Up Procedure\n',n_sig, p)
        perm_FDR(6).method='Multicmp function to apply Hochberg step up procedure';
        perm_FDR(6).n_sig=n_sig;
        perm_FDR(6).p_adj=q';
        
        
        %Multicmp function to apply FDR correction according to Benjamini & Hochberg%%
        [q,~] = multicmp (pch','fdr',p);
        n_sig = sum(q <= p);
        fprintf('%d tests are significant p<=%.2f using FDR correction (Benjamini, 1995)\n',n_sig, p)
        perm_FDR(7).method='Multicmp function to apply FDR correction according to Benjamini & Hochberg';
        perm_FDR(7).n_sig=n_sig;
        perm_FDR(7).p_adj=q';
    catch
        fprintf('Need to download *multicmp.m*\n')
    end
    
    save([savepath date '_resultsPermFDRcorr.mat'],'perm_FDR');
    fprintf('FDR correction save in:\n%s\n\n______________\n*Script completed*\n',[savepath date '_resultsPermFDRcorr.mat'])

end


%% Table graph mode
%
% if tablegraphmode == 1
%     pch=perm_pALL;
%     tblpch = array2table(pch);
%     tblpch.Properties.VariableNames = labelch;
%     tblpch.Properties.RowNames = {'Group'};
%
%     sig = tblpch.Variables <= p;
%     tblpsigch = tblpch(:,sig);
%
%     disp(tblpsigch)
%
%     %writetable(tblsigch,[savepath date '_sigresultsch.xls']); trop gros
%     writetable(tblpsigch,[savepath date '_psig0,05resultsch.xls'],'WriteRowNames',true);
%     save([savepath date '_resultsch0,05.mat'],'tblpch','tblpsigch');
%
%     %%%% graphique%%%%%%%%%%%%%%%%%%%%
%     figure
%     A = meanG1ch(:,sig);
%     B = meanG2ch(:,sig);
%     C = [A' B'];
%     X = categorical(tblpsigch.Properties.VariableNames);
%     p1 = bar(X,C);
%     p1(1).FaceColor = 'r';
%     p1(2).FaceColor = 'b';
%     ylabel('Pearson correlation');
%     savefig([savepath date '_Sig0,05Ch']);
%
%     clear A B C X p1 sig grsig
%
%     %%%% Connectogramme%%%%%%%%%%%%
%     MATmeanG2G1 = MATmeanG2-MATmeanG1; %Calculer matrice G2-G1
%     %MATmeanG1G2 = MATmeanG1 - MATmeanG2;
%
%     for c=1:46
%         for cc =(c + 1):46
%             x = [num2str(c) '-' num2str(cc)];
%             tf = strcmp(x, labelch);
%             idx = find(tf);
%             MATpch(c,cc)= pch(1,idx); %Mettre les  p sous forme de matrice
%         end
%     end
%
%     MATpch(length(MATpch),length(MATpch)) = 0;
%     MATpch = MATpch + triu(MATpch,1)';
%     tf = MATpch <=p & MATpch > 0;
%     %         idx = find(tf);
%     MATsigG2G1ch = MATmeanG2G1;
%     %MATsigG1G2ch = MATmeanG1G2;
%     MATsigG2G1ch(~tf) = 0;
%     %MATsigG1G2ch(~tf) = 0;
%
%     clear c cc x tf idx
%
%     id = 1;
%     MATneg = MATsigG2G1ch < 0;
%     MAT = MATsigG2G1ch.*(MATneg);%DATA{id}.MAT %loader matrice
%     if find(MAT)
%         List = strvcat(DATA{id}.ZoneList); %liste des paires SD
%         ML = DATA{id}.zone.ml; %Loader S/D/ROI/Gr
%         plotLst = DATA{id}.zone.plotLst;
%         label =  DATA{id}.zone.label;
%         fileorderconnectogram  = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Mixte.txt';
%         plotconnectogram(fileorderconnectogram,MAT,List,label,plotLst,ML)
%         savefig([savepath date '_NegConnectSig0,05G2G1Ch']);
%     else
%     end
%
%     MATpos = MATsigG2G1ch > 0;
%     MAT = MATsigG2G1ch.*(MATpos);%DATA{id}.MAT %loader matrice
%     if find(MAT)
%         List = strvcat(DATA{id}.ZoneList); %liste des paires SD
%         ML = DATA{id}.zone.ml; %Loader S/D/ROI/Gr
%         plotLst = DATA{id}.zone.plotLst;
%         label =  DATA{id}.zone.label;
%         fileorderconnectogram  = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Mixte.txt';
%         plotconnectogram(fileorderconnectogram,MAT,List,label,plotLst,ML)
%         savefig([savepath date '_PosConnectSig0,05G2G1Ch']);
%     else
%     end
% end
%
% clear r x res sig A B C X p1

%
% % %%Unpaired ttest%%%%%%%%%%%
% % mean_G1 = squeeze(nanmean(cb1,1));
% % mean_G2 = squeeze(nanmean(pb1,1));
% % var_G1 = squeeze(nanvar(cb1,0,1));
% % var_G2 = squeeze(nanvar(pb1,0,1));
% % n_G1 = squeeze(sum(~isnan(cb1),1));
% % n_G2 = squeeze(sum(~isnan(pb1),1));
% % df = n_G1+ n_G2-2;
% % varc = (1./n_G1+ 1./n_G2 ).*((n_G1-1).*var_G1 + (n_G2-1).* var_G2 )./df;
% % Toij= (mean_G1-mean_G2)*1./sqrt(varc);
% % ncb1 = n_G1;
% % npb1 = n_G2;
% %
% % for i=1:size(Toij,1)
% %     for j=1:size(Toij,2)
% %         try
% %             % Compute the correct p-value for the test
% %             if 1 % two-tailed test pval
% %                 FUniv(i,j) = 2 * tcdf(-abs(- Toij(i,j)), df(i,j));
% %             end
% %
% %         catch
% %             FUniv(i,j) = nan;
% %         end
% %     end
% % end
% %
% % [FDR,Q] = mafdr(FUniv(:));
% %  Q = reshape(Q,size(FUniv));
%
% % MeanG1 = squeeze(nanmean(cb1));
% % MeanG2 = squeeze(nanmean(pb1));
% % %
% % dir1 = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Stats\CORR0,01_0,08\Permutations FDR';
% % % if ~isdir(dir1)
% % %     mkdir(dir1)
% % % end
% %
% % %WRITE IN A NEW FILE
% % infonew = [{'Dir'},{'File'},{'Zone'},{'GR'}];
% % file = [name,'_',labelnode,num2str(NPERM),'permutation tstat','.mat'];
% % matcorr = Toij;
% % meancorr = Toij;
% % save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr');
% % new = [{dir1},{file}, {ZONEid},{1} ];
% % infonew = [infonew;new];
% %
% % file = [name,'_',labelnode,num2str(NPERM),'permutation 1-pval','.mat'];
% % matcorr = 1-FUniv;
% % meancorr = 1-FUniv;
% % save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr');
% % new = [{dir1},{file}, {ZONEid},{1} ];
% % infonew = [infonew;new];
% %
% %
% % file = [name,'_',labelnode,num2str(NPERM),'permutation G1-G2','.mat'];
% % matcorr = real(squeeze((nanmean(cb1,1)- nanmean(pb1,1))));
% % meancorr = real(squeeze((nanmean(cb1,1)- nanmean(pb1,1))));
% % save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr');
% % new = [{dir1},{file}, {ZONEid},{1} ];
% % infonew = [infonew;new];
% %
% %
% % file = [name,'_',labelnode,num2str(NPERM),'permutation G2-G1','.mat'];
% % matcorr = real(squeeze(nanmean(pb1,1))-squeeze(nanmean(cb1,1)));
% % meancorr = real(squeeze(nanmean(pb1,1))-squeeze(nanmean(cb1,1)));
% % save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr');
% % new = [{dir1},{file}, {ZONEid},{1} ];
% % infonew = [infonew;new];
% %
% % file = [name,'_',labelnode,num2str(NPERM),'permutation G1-G2 p05','.mat'];
% % matcorr = real(squeeze(nanmean(cb1,1))-squeeze(nanmean(pb1,1))).*double(FUniv<0.05);
% % meancorr = real(squeeze(nanmean(cb1,1))-squeeze(nanmean(pb1,1))).*double(FUniv<0.05);
% % save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr');
% % new = [{dir1},{file}, {ZONEid},{1} ];
% % infonew = [infonew;new];
% %
% % file = [name,'_',labelnode,num2str(NPERM),'permutation G2-G1 p05','.mat'];
% % matcorr = real(squeeze(nanmean(pb1,1))-squeeze(nanmean(cb1,1))).*double(FUniv<0.05);
% % meancorr = real(squeeze(nanmean(pb1,1))-squeeze(nanmean(cb1,1))).*double(FUniv<0.05);
% % save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr');
% % new = [{dir1},{file}, {ZONEid},{1} ];
% % infonew = [infonew;new];
% %
% % file = [name,'_',labelnode,num2str(NPERM),'permutation mean G1','.mat'];
% % matcorr = real(squeeze(nanmean(cb1,1)));
% % meancorr = real(squeeze(nanmean(cb1,1)));
% % save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr');
% % new = [{dir1},{file}, {ZONEid},{1} ];
% % infonew = [infonew;new];
% %
% % file = [name,'_',labelnode,num2str(NPERM),'permutation mean G2','.mat'];
% % matcorr = real(squeeze(nanmean(pb1,1)));
% % meancorr = real(squeeze(nanmean(pb1,1)));
% % save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr');
% % new = [{dir1},{file}, {ZONEid},{1} ];
% % infonew = [infonew;new];
% %
% %
% % file = [name,'_',labelnode,num2str(NPERM),'permutation N G1','.mat'];
% % matcorr = ncb1;
% % meancorr = ncb1;
% % save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr');
% % new = [{dir1},{file}, {ZONEid},{1} ];
% % infonew = [infonew;new];
% %
% %
% % file = [name,'_',labelnode,num2str(NPERM),'permutation N G2','.mat'];
% % matcorr = npb1;
% % meancorr = npb1;
% % save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr');
% % new = [{dir1},{file}, {ZONEid},{1} ];
% % infonew = [infonew;new];
% %
% % copyfile(fullfile(info{isubject,1}, ZONEid),  fullfile(dir1,  ZONEid))
% toc
%
% %%%%%%%%%%%%%Autres permutations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if otherpermmode == 1
%     %Permutations according to Blair & Karnisky, 1993; Westfall & % Young,1993
%     %de Groppe%
%     disp('Computing Permutations, Groppe')
%     [pch, t_origch, crit_tch, est_alphach, seed_statech] = mult_comp_perm_t2(table2array(tblch(idG1,6:end)),table2array(tblch(idG2,6:end)));
%     [pMroi, t_origMroi, crit_tMroi, est_alphaMroi, seed_stateMroi] = mult_comp_perm_t2(table2array(tblmeanMroi(idG1,6:end)),table2array(tblmeanMroi(idG2,6:end)));
%     [pAroi, t_origAroi, crit_tAroi, est_alphaAroi, seed_stateAroi] = mult_comp_perm_t2(table2array(tblmeanAroi(idG1,6:end)),table2array(tblmeanAroi(idG2,6:end)));
%
%     %%Resampling statistical toolkit, Delorme, 2010%%%
%     disp('Computing Permutations, Delorme')
%     datapath = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Stats\CORR0,01_0,08\';
%     load ([datapath 'workspace.mat'])
%     load ([datapath 'workspacemat.mat'])
%
%     datapath = {MATall(:,:,idG1), MATall(:,:,idG2)};
%
%     [stats, df, pvals, surrog] = statcond(datapath, 'paired','off', 'mode', 'perm');
%
%     %%% Ne semble pas fonctionner, donne 0,05 partout%%%%%%%%%%
%
%     disp('Computing Permutations, Krol')
%     datapath = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Stats\CORR0,01_0,08\';
%     load ([datapath 'workspace.mat'])
%     load ([datapath 'workspacemat.mat'])
%
%     sample1 = chALL(idG1,:);
%     sample2 = chALL(idG2,:);
%     [p, observeddifference, effectsize] = permutationTest(sample1, sample2, 500, 'plotresult',1,'showprogress',1);
%
%     %%Ne fonctionne pas pour des matrices, les données doivent être en vecteur%%%%%%
% end