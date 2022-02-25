%%%%%%%%%%%%%%ANCOVA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
disp('Computing Linear Regression')

datapath = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Stats\CORR0,01_0,08\PhysioSatRespEKG_corr\CorrPair\';
load ([datapath 'workspace.mat'])
%load ([datapath 'workspacemat.mat'])

savepath='C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Stats\CORR0,01_0,08\PhysioSatRespEKG_corr\CorrPair\Regression\';
if ~isfolder(savepath)
    mkdir(savepath)
end
channelmode = 1;
importedROImode = 0;
calculatedROImode = 1;

graphmode = 1;
%fileorderconnectogram  = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Mixte.txt';
fdrmode = 1;
%nointeractionmode = 1;
%interactionmode = 1;
p = 0.07;

%%CHANNELS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if channelmode ==1
    disp('Computing Linear Regression on channels')
    dat = lgndch{1,:};
    %%REGRESSION%%%%%%%%%%%%%%%%%%%
    x = 1;
    coefch = [];
    resch = [];
    pvalch = [];
    for c = 1:(width(tblch)- 5)
        mdl = fitlm([gr ses],table2array(tblch(:,5+c)),'interactions','CategoricalVars',1,'VarNames',{'Group','SES','Connectivity'});
        tblcoef = mdl.Coefficients;
        %col{1,1} = []; col = [col; tbl.Properties.RowNames]; coef = [col [tbl.Properties.VariableNames; table2cell(tbl)]];
        coefch{1,x} = tblcoef; %[coefch, coef];
        tblres = anova(mdl,'components');
        %row{1,1} = []; row = [row tblres.Properties.VariableNames]; res = [row; tblres.Properties.RowNames table2cell(tblres)]; 
        resch{1,x} = tblres; %[resch, tblres];
        pvalch(:,x) = tblres.pValue(1:end-1);
        x = x + 1;
        clear res coef tblcoef tbres mdl
    end
    
    labelfactors = resch{1,1}.Properties.RowNames(1:end-1,1);
    for f = 1:numel(labelfactors)
        n_sig = sum(pvalch(f,:) <= p);
        fprintf('%d %s effects are significant p<=%.2f without correction for channels\n',n_sig,labelfactors{f,1}, p)
    end
    
    clear c f x n_sig
    
    %%fdr%%%%%%%%%%%%
    if fdrmode == 1
        FDR_CO(pvalch, labelfactors, p, dat, savepath)
    end
    
    if any(contains(labelfactors,':')) %si terme interaction
        idx = find(contains(labelfactors,':')); %trouver le num du facteur
        idsiginter = find(pvalch(idx,:)<= p); %trouver les interactions sig
        coefchsiginterG1 = [];
        coefchsiginterG2 = [];
        pvalchsiginterG1 = [];
        pvalchsiginterG2 = [];
        x = 1;
        for c = 1:numel(idsiginter) %pour chaque inter sig
            id = idsiginter(c); %nom de la paire de canaux
            mdl = fitlm(ses(idG1),table2array(tblch(idG1,5+id)),'VarNames',{'SES','Connectivity'});
            tblcoefG1 = mdl.Coefficients;
            coefchsiginterG1{1,x} = tblcoefG1;
            pvalchsiginterG1(:,x) = tblcoefG1.pValue(2:end);
            clear mdl tblcoefG1
            mdl = fitlm(ses(idG2),table2array(tblch(idG2,5+id)),'VarNames',{'SES','Connectivity'});
            tblcoefG2 = mdl.Coefficients;
            coefchsiginterG2{1,x} = tblcoefG2;
            pvalchsiginterG2(:,x) = tblcoefG2.pValue(2:end);
            x = x + 1;
            clear mdl tblcoefG2
        end
        clear c x n_sig
    end
    
    %%tableaux des résultats%%%%%%%%%%%%%%%% 
    [~,~,~,~,~,~] = Tablegraph_CO(resch, pvalch, labelfactors, p, dat, lgndch, labelch, MATmeanG1G2, meanG1ch, meanG2ch, DATA, fileorderconnectogram, savepath, graphmode, coefch, coefchsiginterG1, coefchsiginterG2, pvalchsiginterG1, pvalchsiginterG2);  
end

if calculatedROImode ==1
    disp('Computing Linear Regression on ROIs')
    
     %%MROI%%%%%%%%%%%%%%%%%%%%%%%
    R = 1;
    dat = lgndroi{1,R};
    x = 1;
    
    coefMroi = [];
    resMroi = [];
    pvalMroi = [];
    for r = 1:(width(tbldataMroi)- 5)
        mdl = fitlm([gr ses],table2array(tblch(:,5+r)),'interactions','CategoricalVars',1,'VarNames',{'Group','SES','Connectivity'});
        tblcoef = mdl.Coefficients;
        coefMroi{1,x} = tblcoef; %[coefMroi, coef];
        tblres = anova(mdl,'components');
        resMroi{1,x} = tblres; %[resMroi, res];
        pvalMroi(:,x) = tblres.pValue(1:end-1);
        x = x + 1;
        clear res coef tblcoef tbres mdl
    end
    
    labelfactors = resMroi{1,1}.Properties.RowNames(1:end-1,1);
    for f = 1:numel(labelfactors)
        n_sig = sum(pvalMroi(f,:) <= p);
        fprintf('%d %s effects are significant p<=%.2f without correction for Mroi\n',n_sig,labelfactors{f,1}, p)
    end
    
    clear r x f n_sig
    
    %%fdr%%%%%%%%%%%%
    if fdrmode == 1
        FDR_CO(pvalMroi, labelfactors, p, dat, savepath)
    end
    
    if any(contains(labelfactors,':')) %si terme interaction
        idx = find(contains(labelfactors,':')); %trouver le num du facteur
        idsiginter = find(pvalMroi(idx,:)<= p); %trouver les interactions sig
        coefMroisiginterG1 = [];
        coefMroisiginterG2 = [];
        pvalMroisiginterG1 = [];
        pvalMroisiginterG2 = [];
        x = 1;
        for r = 1:numel(idsiginter) %pour chaque inter sig
            id = idsiginter(r); %nom de la paire de canaux
            mdl = fitlm(ses(idG1),table2array(tbldataMroi(idG1,5+id)),'VarNames',{'SES','Connectivity'});
            tblcoefG1 = mdl.Coefficients;
            coefMroisiginterG1{1,x} = tblcoefG1;
            pvalMroisiginterG1(:,x) = tblcoefG1.pValue(2:end);
            clear mdl tblcoefG1
            mdl = fitlm(ses(idG2),table2array(tbldataMroi(idG2,5+id)),'VarNames',{'SES','Connectivity'});
            tblcoefG2 = mdl.Coefficients;
            coefMroisiginterG2{1,x} = tblcoefG2;
            pvalMroisiginterG2(:,x) = tblcoefG2.pValue(2:end);
            x = x + 1;
            clear mdl tblcoefG2
        end
    end
    
    clear r x n_sig
    
    %%tableaux des résultats%%%%%%%%%%%%%%%%
    if tablegraphmode == 1
        [~,~,~,~,~,~] = Tablegraph_CO(resMroi, pvalMroi, labelfactors, p, dat, lgndroi, labelroiALL{1,R}, MATroimeanG1G2{1,R}, meanroiG1r{1,R}, meanroiG2r{1,R}, roi{1,R}, fileorderconnectogram, savepath, graphmode, coefMroi, coefMroisiginterG1, coefMroisiginterG2, pvalMroisiginterG1, pvalMroisiginterG2);
    end
    
    %%RROI%%%%%%%%%%%%%%%%%%%%%%%
    R = 2;
    dat = lgndroi{1,R};
    x = 1;
    
    coefRroi = [];
    resRroi = [];
    pvalRroi = [];
    for r = 1:(width(tbldataRroi)- 5)
        mdl = fitlm([gr ses],table2array(tblch(:,5+r)),'interactions','CategoricalVars',1,'VarNames',{'Group','SES','Connectivity'});
        tblcoef = mdl.Coefficients;
        coefRroi{1,x} = tblcoef; %[coefRroi, coef];
        tblres = anova(mdl,'components');
        resRroi{1,x} = tblres; %[resRroi, res];
        pvalRroi(:,x) = tblres.pValue(1:end-1);
        x = x + 1;
        clear res coef tblcoef tbres mdl
    end
    
    labelfactors = resRroi{1,1}.Properties.RowNames(1:end-1,1);
    for f = 1:numel(labelfactors)
        n_sig = sum(pvalRroi(f,:) <= p);
        fprintf('%d %s effects are significant p<=%.2f without correction for Rroi\n',n_sig,labelfactors{f,1}, p)
    end
    
    clear r x f n_sig
    
    %%fdr%%%%%%%%%%%%
    if fdrmode == 1
        FDR_CO(pvalRroi, labelfactors, p, dat, savepath)
    end
    
    if any(contains(labelfactors,':')) %si terme interaction
        idx = find(contains(labelfactors,':')); %trouver le num du facteur
        idsiginter = find(pvalRroi(idx,:)<= p); %trouver les interactions sig
        coefRroisiginterG1 = [];
        coefRroisiginterG2 = [];
        pvalRroisiginterG1 = [];
        pvalRroisiginterG2 = [];
        x = 1;
        for r = 1:numel(idsiginter) %pour chaque inter sig
            id = idsiginter(r); %nom de la paire de canaux
            mdl = fitlm(ses(idG1),table2array(tbldataRroi(idG1,5+id)),'VarNames',{'SES','Connectivity'});
            tblcoefG1 = mdl.Coefficients;
            coefRroisiginterG1{1,x} = tblcoefG1;
            pvalRroisiginterG1(:,x) = tblcoefG1.pValue(2:end);
            clear mdl tblcoefG1
            mdl = fitlm(ses(idG2),table2array(tbldataRroi(idG2,5+id)),'VarNames',{'SES','Connectivity'});
            tblcoefG2 = mdl.Coefficients;
            coefRroisiginterG2{1,x} = tblcoefG2;
            pvalRroisiginterG2(:,x) = tblcoefG2.pValue(2:end);
            x = x + 1;
            clear mdl tblcoefG2
        end
    end
    
    clear r x n_sig
    
    %%tableaux des résultats%%%%%%%%%%%%%%%%
    if tablegraphmode == 1
        [~,~,~,~,~,~] = Tablegraph_CO(resRroi, pvalRroi, labelfactors, p, dat, lgndroi, labelroiALL{1,R}, MATroimeanG1G2{1,R}, meanroiG1r{1,R}, meanroiG2r{1,R}, roi{1,R}, fileorderconnectogram, savepath, graphmode, coefRroi, coefRroisiginterG1, coefRroisiginterG2, pvalRroisiginterG1, pvalRroisiginterG2);
    end
    
    %%FROI%%%%%%%%%%%%%%%%%%%%%%%
    R = 3;
    dat = lgndroi{1,R};
    x = 1;
    
    coefFroi = [];
    resFroi = [];
    pvalFroi = [];
    for r = 1:(width(tbldataFroi)- 5)
        mdl = fitlm([gr ses],table2array(tblch(:,5+r)),'interactions','CategoricalVars',1,'VarNames',{'Group','SES','Connectivity'});
        tblcoef = mdl.Coefficients;
        coefFroi{1,x} = tblcoef; %[coefFroi, coef];
        tblres = anova(mdl,'components');
        resFroi{1,x} = tblres; %[resFroi, res];
        pvalFroi(:,x) = tblres.pValue(1:end-1);
        x = x + 1;
        clear res coef tblcoef tbres mdl
    end
    
    labelfactors = resFroi{1,1}.Properties.RowNames(1:end-1,1);
    for f = 1:numel(labelfactors)
        n_sig = sum(pvalFroi(f,:) <= p);
        fprintf('%d %s effects are significant p<=%.2f without correction for Froi\n',n_sig,labelfactors{f,1}, p)
    end
    
    clear r x f n_sig
    
    %%fdr%%%%%%%%%%%%
    if fdrmode == 1
        FDR_CO(pvalFroi, labelfactors, p, dat, savepath)
    end
    
    if any(contains(labelfactors,':')) %si terme interaction
        idx = find(contains(labelfactors,':')); %trouver le num du facteur
        idsiginter = find(pvalFroi(idx,:)<= p); %trouver les interactions sig
        coefFroisiginterG1 = [];
        coefFroisiginterG2 = [];
        pvalFroisiginterG1 = [];
        pvalFroisiginterG2 = [];
        x = 1;
        for r = 1:numel(idsiginter) %pour chaque inter sig
            id = idsiginter(r); %nom de la paire de canaux
            mdl = fitlm(ses(idG1),table2array(tbldataFroi(idG1,5+id)),'VarNames',{'SES','Connectivity'});
            tblcoefG1 = mdl.Coefficients;
            coefFroisiginterG1{1,x} = tblcoefG1;
            pvalFroisiginterG1(:,x) = tblcoefG1.pValue(2:end);
            clear mdl tblcoefG1
            mdl = fitlm(ses(idG2),table2array(tbldataFroi(idG2,5+id)),'VarNames',{'SES','Connectivity'});
            tblcoefG2 = mdl.Coefficients;
            coefFroisiginterG2{1,x} = tblcoefG2;
            pvalFroisiginterG2(:,x) = tblcoefG2.pValue(2:end);
            x = x + 1;
            clear mdl tblcoefG2
        end
    end
    
    clear r x n_sig
    
    %%tableaux des résultats%%%%%%%%%%%%%%%%
    if tablegraphmode == 1
        [~,~,~,~,~,~] = Tablegraph_CO(resFroi, pvalFroi, labelfactors, p, dat, lgndroi, labelroiALL{1,R}, MATroimeanG1G2{1,R}, meanroiG1r{1,R}, meanroiG2r{1,R}, roi{1,R}, fileorderconnectogram, savepath, graphmode, coefFroi, coefFroisiginterG1, coefFroisiginterG2, pvalFroisiginterG1, pvalFroisiginterG2);
    end
    
    %%AROI%%%%%%%%%%%%%%%%%%%%%%%
    R = 4;
    dat = lgndroi{1,R};
    x = 1;
    
    coefAroi = [];
    resAroi = [];
    pvalAroi = [];
    for r = 1:(width(tbldataAroi)- 5)
        mdl = fitlm([gr ses],table2array(tblch(:,5+r)),'interactions','CategoricalVars',1,'VarNames',{'Group','SES','Connectivity'});
        tblcoef = mdl.Coefficients;
        coefAroi{1,x} = tblcoef; %[coefAroi, coef];
        tblres = anova(mdl,'components');
        resAroi{1,x} = tblres; %[resAroi, res];
        pvalAroi(:,x) = tblres.pValue(1:end-1);
        x = x + 1;
        clear res coef tblcoef tbres mdl
    end
    
    labelfactors = resAroi{1,1}.Properties.RowNames(1:end-1,1);
    for f = 1:numel(labelfactors)
        n_sig = sum(pvalAroi(f,:) <= p);
        fprintf('%d %s effects are significant p<=%.2f without correction for Aroi\n',n_sig,labelfactors{f,1}, p)
    end
    
    clear r x f n_sig
    
    %%fdr%%%%%%%%%%%%
    if fdrmode == 1
        FDR_CO(pvalAroi, labelfactors, p, dat, savepath)
    end
    
    if any(contains(labelfactors,':')) %si terme interaction
        idx = find(contains(labelfactors,':')); %trouver le num du facteur
        idsiginter = find(pvalAroi(idx,:)<= p); %trouver les interactions sig
        coefAroisiginterG1 = [];
        coefAroisiginterG2 = [];
        pvalAroisiginterG1 = [];
        pvalAroisiginterG2 = [];
        x = 1;
        for r = 1:numel(idsiginter) %pour chaque inter sig
            id = idsiginter(r); %nom de la paire de canaux
            mdl = fitlm(ses(idG1),table2array(tbldataAroi(idG1,5+id)),'VarNames',{'SES','Connectivity'});
            tblcoefG1 = mdl.Coefficients;
            coefAroisiginterG1{1,x} = tblcoefG1;
            pvalAroisiginterG1(:,x) = tblcoefG1.pValue(2:end);
            clear mdl tblcoefG1
            mdl = fitlm(ses(idG2),table2array(tbldataAroi(idG2,5+id)),'VarNames',{'SES','Connectivity'});
            tblcoefG2 = mdl.Coefficients;
            coefAroisiginterG2{1,x} = tblcoefG2;
            pvalAroisiginterG2(:,x) = tblcoefG2.pValue(2:end);
            x = x + 1;
            clear mdl tblcoefG2
        end
    end
    
    clear r x n_sig
    
    %%tableaux des résultats%%%%%%%%%%%%%%%%
    if tablegraphmode == 1
        [~,~,~,~,~,~] = Tablegraph_CO(resAroi, pvalAroi, labelfactors, p, dat, lgndroi, labelroiALL{1,R}, MATroimeanG1G2{1,R}, meanroiG1r{1,R}, meanroiG2r{1,R}, roi{1,R}, fileorderconnectogram, savepath, graphmode, coefAroi, coefAroisiginterG1, coefAroisiginterG2, pvalAroisiginterG1, pvalAroisiginterG2);
    end
end

 X = ['Results saved in ', savepath];
disp(X)
clear   
end

toc