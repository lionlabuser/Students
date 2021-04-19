%%%%%%%%%%%%%%%%%%%%%%%%%PERMUTATIONS/UNPAIRED TTEST%%%%%%%%%%%%%%%%%%
tic

%%Permutation of Pedro's team from StatMatrices of LIONIRS toolbox%%%
disp('Computing Permutations, Pedro et al.')
datapath = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Stats\CORR0,01_0,08\';
load ([datapath 'workspace.mat'])
load ([datapath 'workspacemat.mat'])

savepath='C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Stats\CORR0,01_0,08\Permutations FDR\';
FDRmode = 1;
tablegraphmode = 1;
otherpermmode = 0;
p = 0.05;

ESTAD = 4;
INDEP = 1;
NPERM = 500;

MATall = permute(MATall, [3 2 1]);
cb1 = MATall(idG1',:,:);
pb1 = MATall(idG2',:,:);

%initialised
Toij = zeros(size(cb1,2));
ncb1 = zeros(size(cb1,2));
npb1 = zeros(size(pb1,2));

%Permutations%%%
for i=1:size(pb1,2) % row
    for j=1:size(cb1,3) %col
        tmpcb1 = cb1(:,i,j);
        tmppb1 = pb1(:,i,j);
        idnan = find(isnan(tmpcb1));
        if ~isempty(idnan)
            tmpcb1(idnan)=[];
        end
        idnan = find(isnan(tmppb1));
        if ~isempty(idnan)
            tmppb1(idnan)=[];
        end
        if Toij(i,j)==0
            [FSupSup,FSupDeriv,FSupTime,pij,tij] = TestPermut2Grupos(ESTAD,INDEP,tmpcb1,tmppb1,NPERM);
            %apply symetric
            ncb1(i,j) = numel(tmpcb1);
            npb1(i,j) = numel(tmppb1);
            FUniv(i,j) = pij;
            Toij(i,j) = tij;
            ncb1(j,i) = numel(tmpcb1);
            npb1(j,i) = numel(tmppb1);
            FUniv(j,i) = pij;
            Toij(j,i) = tij;
            
        end
    end
end

save([savepath date '_resultsperm.mat'],'ncb1','npb1','FUniv','Toij');

%%tableaux des résultats%%%%%%%%%%%%%%%%
if tablegraphmode == 1
        
    x = 1;
    for c = 1:46
        for cc = (c + 1):46
            labelch{x} = [num2str(c) '-' num2str(cc)];
            pch(x) = FUniv(c,cc);
            x = x + 1; 
        end
    end
    
        n_sig = sum(pch <= p);
        fprintf('%d tests out of %d are significant p<=%.2f without correction\n',n_sig,numel(pch), p)
        
        tblpch = array2table(pch);
        tblpch.Properties.VariableNames = labelch;
        tblpch.Properties.RowNames = {'Group'};

        sig = tblpch.Variables <= p;
        tblpsigch = tblpch(:,sig);

        disp(tblpsigch)

        %writetable(tblsigch,[savepath date '_sigresultsch.xls']); trop gros
        writetable(tblpsigch,[savepath date '_psig0,05resultsch.xls'],'WriteRowNames',true);
        save([savepath date '_resultsch0,05.mat'],'tblpch','tblpsigch');

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

clear r x res sig A B C X p1

if FDRmode == 1
   
    %%FDR according to Storey 2002%%%
    %[fdr,q] = mafdr(FUniv(:));
    %q = reshape(q,size(FUniv));
    %n_sig = sum(sum(q <= 0.05));
    %fprintf('%d tests are significant p<=0.05 using FDR correction\n',n_sig)

    % % FDR according to Storey 2002%%%
    [fdr,q] = mafdr(pch(1,:)');
    n_sig = sum(q <= p);
    fprintf('%d tests are significant p<=%.2f using FDR correction (Storey, 2002)\n',n_sig, p)

    % FDR according to Benjamini & Hochberg(1995) de Groppe%%%%
    [h1, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pch(1,:)', 0.05, 'pdep', 'no' );
    n_sig = sum(adj_p <= p);
    fprintf('%d tests are significant p<=%.2f using FDR correction (Benjamini,1995)\n',n_sig, p)

    % % FDR according to Benjamini & Hochberg(1995) de Gerber%%%%
    %[ind, thres] = FDR(pvalues(1,:));

    %Bonferroni-Holm correction by Groppe%%% Corrected_p are greater than 1... doesn't make sense
    [corrected_p, h] = bonf_holm(pch(1,:)', 0.05);
    n_sig = sum(corrected_p <= p);
    fprintf('%d tests are significant p<=%.2f using Bonferroni-Holm correction\n',n_sig, p)

    %Multicmp function to apply Holm Step Down Procedure%%
    [padjBH,alphaBH] = multicmp (pch(1,:)','down',0.05);
    n_sig = sum(padjBH <= p);
    fprintf('%d tests are significant p<=%.2f using Holm Step Down Procedure\n',n_sig, p)

    %Multicmp function to apply Hochberg's step up procedure%%
    [padjH,alphaH] = multicmp (pch(1,:)','up',0.05);
    n_sig = sum(padjH <= p);
    fprintf('%d tests are significant p<=%.2f using Hochberg Step Up Procedure\n',n_sig, p)

    %Multicmp function to apply FDR correction according to Benjamini & Hochberg%%
    [padjFDR,alphaFDR] = multicmp (pch(1,:)','fdr',0.05);
    n_sig = sum(padjFDR <= p);
    fprintf('%d tests are significant p<=%.2f using FDR correction (Benjamini, 1995)\n',n_sig, p)
end

% %%Unpaired ttest%%%%%%%%%%%
% mean_G1 = squeeze(nanmean(cb1,1));
% mean_G2 = squeeze(nanmean(pb1,1));
% var_G1 = squeeze(nanvar(cb1,0,1));
% var_G2 = squeeze(nanvar(pb1,0,1));
% n_G1 = squeeze(sum(~isnan(cb1),1));
% n_G2 = squeeze(sum(~isnan(pb1),1));
% df = n_G1+ n_G2-2;
% varc = (1./n_G1+ 1./n_G2 ).*((n_G1-1).*var_G1 + (n_G2-1).* var_G2 )./df;
% Toij= (mean_G1-mean_G2)*1./sqrt(varc);
% ncb1 = n_G1;
% npb1 = n_G2;
% 
% for i=1:size(Toij,1)
%     for j=1:size(Toij,2)
%         try
%             % Compute the correct p-value for the test
%             if 1 % two-tailed test pval
%                 FUniv(i,j) = 2 * tcdf(-abs(- Toij(i,j)), df(i,j));
%             end
%             
%         catch
%             FUniv(i,j) = nan;
%         end
%     end
% end
%
% [FDR,Q] = mafdr(FUniv(:));
%  Q = reshape(Q,size(FUniv));

% MeanG1 = squeeze(nanmean(cb1));
% MeanG2 = squeeze(nanmean(pb1));
% % 
% dir1 = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Stats\CORR0,01_0,08\Permutations FDR';
% % if ~isdir(dir1)
% %     mkdir(dir1)
% % end
% 
% %WRITE IN A NEW FILE
% infonew = [{'Dir'},{'File'},{'Zone'},{'GR'}];
% file = [name,'_',labelnode,num2str(NPERM),'permutation tstat','.mat'];
% matcorr = Toij;
% meancorr = Toij;
% save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr');
% new = [{dir1},{file}, {ZONEid},{1} ];
% infonew = [infonew;new];
% 
% file = [name,'_',labelnode,num2str(NPERM),'permutation 1-pval','.mat'];
% matcorr = 1-FUniv;
% meancorr = 1-FUniv;
% save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr');
% new = [{dir1},{file}, {ZONEid},{1} ];
% infonew = [infonew;new];
% 
% 
% file = [name,'_',labelnode,num2str(NPERM),'permutation G1-G2','.mat'];
% matcorr = real(squeeze((nanmean(cb1,1)- nanmean(pb1,1))));
% meancorr = real(squeeze((nanmean(cb1,1)- nanmean(pb1,1))));
% save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr');
% new = [{dir1},{file}, {ZONEid},{1} ];
% infonew = [infonew;new];
% 
% 
% file = [name,'_',labelnode,num2str(NPERM),'permutation G2-G1','.mat'];
% matcorr = real(squeeze(nanmean(pb1,1))-squeeze(nanmean(cb1,1)));
% meancorr = real(squeeze(nanmean(pb1,1))-squeeze(nanmean(cb1,1)));
% save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr');
% new = [{dir1},{file}, {ZONEid},{1} ];
% infonew = [infonew;new];
% 
% file = [name,'_',labelnode,num2str(NPERM),'permutation G1-G2 p05','.mat'];
% matcorr = real(squeeze(nanmean(cb1,1))-squeeze(nanmean(pb1,1))).*double(FUniv<0.05);
% meancorr = real(squeeze(nanmean(cb1,1))-squeeze(nanmean(pb1,1))).*double(FUniv<0.05);
% save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr');
% new = [{dir1},{file}, {ZONEid},{1} ];
% infonew = [infonew;new];
% 
% file = [name,'_',labelnode,num2str(NPERM),'permutation G2-G1 p05','.mat'];
% matcorr = real(squeeze(nanmean(pb1,1))-squeeze(nanmean(cb1,1))).*double(FUniv<0.05);
% meancorr = real(squeeze(nanmean(pb1,1))-squeeze(nanmean(cb1,1))).*double(FUniv<0.05);
% save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr');
% new = [{dir1},{file}, {ZONEid},{1} ];
% infonew = [infonew;new];
% 
% file = [name,'_',labelnode,num2str(NPERM),'permutation mean G1','.mat'];
% matcorr = real(squeeze(nanmean(cb1,1)));
% meancorr = real(squeeze(nanmean(cb1,1)));
% save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr');
% new = [{dir1},{file}, {ZONEid},{1} ];
% infonew = [infonew;new];
% 
% file = [name,'_',labelnode,num2str(NPERM),'permutation mean G2','.mat'];
% matcorr = real(squeeze(nanmean(pb1,1)));
% meancorr = real(squeeze(nanmean(pb1,1)));
% save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr');
% new = [{dir1},{file}, {ZONEid},{1} ];
% infonew = [infonew;new];
% 
% 
% file = [name,'_',labelnode,num2str(NPERM),'permutation N G1','.mat'];
% matcorr = ncb1;
% meancorr = ncb1;
% save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr');
% new = [{dir1},{file}, {ZONEid},{1} ];
% infonew = [infonew;new];
% 
% 
% file = [name,'_',labelnode,num2str(NPERM),'permutation N G2','.mat'];
% matcorr = npb1;
% meancorr = npb1;
% save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr');
% new = [{dir1},{file}, {ZONEid},{1} ];
% infonew = [infonew;new];
% 
% copyfile(fullfile(info{isubject,1}, ZONEid),  fullfile(dir1,  ZONEid))
toc

%%%%%%%%%%%%%Autres permutations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if otherpermmode == 1
    %Permutations according to Blair & Karnisky, 1993; Westfall & % Young,1993
    %de Groppe%
    disp('Computing Permutations, Groppe')
    [pch, t_origch, crit_tch, est_alphach, seed_statech] = mult_comp_perm_t2(table2array(tblch(idG1,6:end)),table2array(tblch(idG2,6:end)));
    [pMroi, t_origMroi, crit_tMroi, est_alphaMroi, seed_stateMroi] = mult_comp_perm_t2(table2array(tblmeanMroi(idG1,6:end)),table2array(tblmeanMroi(idG2,6:end)));
    [pAroi, t_origAroi, crit_tAroi, est_alphaAroi, seed_stateAroi] = mult_comp_perm_t2(table2array(tblmeanAroi(idG1,6:end)),table2array(tblmeanAroi(idG2,6:end)));

    %%Resampling statistical toolkit, Delorme, 2010%%%
    disp('Computing Permutations, Delorme')
    datapath = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Stats\CORR0,01_0,08\';
    load ([datapath 'workspace.mat'])
    load ([datapath 'workspacemat.mat'])

    datapath = {MATall(:,:,idG1), MATall(:,:,idG2)};

    [stats, df, pvals, surrog] = statcond(datapath, 'paired','off', 'mode', 'perm');

    %%% Ne semble pas fonctionner, donne 0,05 partout%%%%%%%%%%

    disp('Computing Permutations, Krol')
    datapath = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Stats\CORR0,01_0,08\';
    load ([datapath 'workspace.mat'])
    load ([datapath 'workspacemat.mat'])

    sample1 = chALL(idG1,:);
    sample2 = chALL(idG2,:);
    [p, observeddifference, effectsize] = permutationTest(sample1, sample2, 500, 'plotresult',1,'showprogress',1);

    %%Ne fonctionne pas pour des matrices, les données doivent être en vecteur%%%%%%
end