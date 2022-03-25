%%%%%%%%%%%%%%ANCOVA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

datapath = 'C:\data\Malnutrition\Resting\NIRS\Analyses\Stats\CORR0,01_0,08\PCAPW_nocriteria\fisher\';
load ([datapath 'workspace.mat'])
%load ([datapath 'workspacemat.mat'])

disp(['Computing ANCOVA on ' datapath])

savepathinitial = fullfile(fileparts(datapath), filesep, 'ANCOVA', filesep);
if ~isfolder(savepathinitial)
    mkdir(savepathinitial)
end
channelmode = 1;
importedROImode = 0;
calculatedROImode = 1;

graphmode = 1;
fdrmode = 1;
nointeractionmode = 1;
interactionmode = 1;
p = 0.07;

fileorderconnectogram = {'C:\data\Malnutrition\Resting\NIRS\Analyses\Connectogram_Mixte.txt',...
    'C:\data\Malnutrition\Resting\NIRS\Analyses\Connectogram_Region.txt',...
    'C:\data\Malnutrition\Resting\NIRS\Analyses\Connectogram_Fonction.txt',...
    'C:\data\Malnutrition\Resting\NIRS\Analyses\Connectogram_Aire.txt'};


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
        dat = lgndch{1,:};

        %%ANOVA%%%%%%%%%%%%%%%%%%%
        x = 1;
        resch = [];
        pvalch = [];
        for r = 1:(width(tblch)- 5)
            [pvalch(:,x),res] = anovan(table2array(tblch(:,r+5)),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
            close hidden;
            resc = res(2:end,2:end);
            idx = find(cellfun(@isempty,resc));
            resc(idx) = {[0]};
            tblres = array2table(cell2mat(reshape(resc,4,6)),'RowNames',res(2:end,1),'VariableNames',res(1,2:end));
            resch{1,x} = tblres; %[resch, tblres];
            x = x + 1;
        end
        clear x r res idx resc tblres
        % delete(findall(0));

        %       [~,atab,ctab,stats] = aoctool(ses,table2array(tblch(:,7)),gr,0.05,'SES','Connectivity','Group','on','separate lines');
        %       [c,m] = multcompare(stats,'Estimate','slope');
        %       [c,m] = multcompare(stats,'Estimate','pmm');
        %       [c,m] = multcompare(stats,'Estimate','intercept');

        labelfactors = resch{1,1}.Properties.RowNames(1:end-2,1);
        for f = 1:numel(labelfactors)
            n_sig = sum(pvalch(f,:) <= p);
            fprintf('%d %s effects are significant p<=%.2f without correction for channels\n',n_sig,labelfactors{f,1}, p)
        end
        clear r x res n_sig

        %%fdr%%%%%%%%%%%%
        if fdrmode == 1
            FDR_CO(pvalch, labelfactors, p, dat, savepath)
        end

        %%tableaux des résultats%%%%%%%%%%%%%%%%
        [~,~,~,~,~,~] = Tablegraph_CO(resch, pvalch, labelfactors, p, dat, lgndch, ...
            labelch, MATallG1G2, datameanchG1, datameanchG2, DATA, fileorderconnectogram, ...
            savepath, graphmode, ListChinv);
    end

    if importedROImode
        dat = lgndch{1,:};
        R = 1;
        %%ANOVA%%%%%%%%%%%%%%%%%%%
        x = 1;
        resroi = [];
        for r = 1:(width(tblroi)- 5)
            [pvalroi(:,x),res] = anovan(table2array(tblroi(:,r+5)),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
            close hidden;
            resroi = [resroi, res];
            x = x + 1;
        end
        % delete(findall(0));

        labelfactors = resroi(2:end-2,1);
        for f = 1:numel(labelfactors)
            n_sig = sum(pvalroi(f,:) <= p);
            fprintf('%d %s effects are significant p<=%.2f without correction for ROIs\n',n_sig,labelfactors{f,1}, p)
        end

        clear r x res n_sig

        %%fdr%%%%%%%%%%%%
        if fdrmode == 1
            FDR_CO(pvalroi, labelfactors, p, dat, savepath)
        end

        %%tableaux des résultats%%%%%%%%%%%%%%%%
        if graphmode == 1
            %[~,~,~,~,~,~] = Tablegraph_CO(resroi, pvalroi, p, dat, lgndroi, labelroiALL{1,R}, MATroiG1G2{1,R}, datameanroiG1{1,R}, datameanroiG2{1,R}, roi{1,R}, fileorderconnectogram, savepath);

            tblresroi = array2table(resroi);
            for r = 1:(width(tblroi)- 5)
                tblresroi = mergevars(tblresroi,(r:r+6));
            end
            tblresroi.Properties.VariableNames = labelroi;
            tblresroi.Properties.RowNames = {'Source','Group','SES','Error','Total'};

            tblproi = array2table(pvalroi);
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
                    MATproi(r,rr)= pvalroi(1,idx); %Mettre les  p sous forme de matrice
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
        dat = lgndroi{1,R};
        x = 1;
        resMroi = [];
        pvalMroi = [];
        for r = 1:(width(tbldataMroi)- 5)
            [pvalMroi(:,x),res] = anovan(table2array(tbldataMroi(:,r+5)),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
            close hidden;
            resc = res(2:end,2:end);
            idx = find(cellfun(@isempty,resc));
            resc(idx) = {[0]};
            tblres = array2table(cell2mat(reshape(resc,4,6)),'RowNames',res(2:end,1),'VariableNames',res(1,2:end));
            resMroi{1,x} = tblres; %[resch, tblres];
            x = x + 1;
        end
        clear x r res idx resc tblres
        % delete(findall(0));

        labelfactors = resMroi{1,1}.Properties.RowNames(1:end-2,1);
        for f = 1:numel(labelfactors)
            n_sig = sum(pvalMroi(f,:) <= p);
            fprintf('%d %s effects are significant p<=%.2f without correction for Mroi\n',n_sig,labelfactors{f,1}, p)
        end
        clear r x res n_sig

        %%fdr%%%%%%%%%%%%
        if fdrmode == 1
            FDR_CO(pvalMroi, labelfactors, p, dat, savepath)
        end

        %%tableaux des résultats%%%%%%%%%%%%%%%%
        [~,~,~,~,~,~] = Tablegraph_CO(resMroi, pvalMroi, labelfactors, p, dat, lgndroi,...
            labelroiALL{1,R}, MATroiG1G2{1,R}, datameanroiG1{1,R}, datameanroiG2{1,R},...
            roi{1,R}, fileorderconnectogram, savepath, graphmode);

        %%RROI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%ANOVA%%%%%%%%%%%%%%%%%%%
        R = 2;
        dat = lgndroi{1,R};
        x = 1;
        resRroi = [];
        pvalRroi = [];
        for r = 1:(width(tbldataRroi)- 5)
            [pvalRroi(:,x),res] = anovan(table2array(tbldataRroi(:,r+5)),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
            close hidden;
            resc = res(2:end,2:end);
            idx = find(cellfun(@isempty,resc));
            resc(idx) = {[0]};
            tblres = array2table(cell2mat(reshape(resc,4,6)),'RowNames',res(2:end,1),'VariableNames',res(1,2:end));
            resRroi{1,x} = tblres; %[resch, tblres];
            x = x + 1;
        end
        clear x r res idx resc tblres
        % delete(findall(0));

        labelfactors = resRroi{1,1}.Properties.RowNames(1:end-2,1);
        for f = 1:numel(labelfactors)
            n_sig = sum(pvalRroi(f,:) <= p);
            fprintf('%d %s effects are significant p<=%.2f without correction for Rrois\n',n_sig,labelfactors{f,1}, p)
        end
        clear r x res n_sig

        %%fdr%%%%%%%%%%%%
        if fdrmode == 1
            FDR_CO(pvalRroi, labelfactors, p, dat, savepath)
        end

        %%tableaux des résultats%%%%%%%%%%%%%%%%
        [~,~,~,~,~,~] = Tablegraph_CO(resRroi, pvalRroi, labelfactors, p, dat, ...
            lgndroi, labelroiALL{1,R}, MATroiG1G2{1,R}, datameanroiG1{1,R}, ...
            datameanroiG2{1,R}, roi{1,R}, fileorderconnectogram, savepath, graphmode);

        %%FROI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%ANOVA%%%%%%%%%%%%%%%%%%%
        R = 3;
        dat = lgndroi{1,R};
        x = 1;
        resFroi = [];
        pvalFroi = [];
        for r = 1:(width(tbldataFroi)- 5)
            [pvalFroi(:,x),res] = anovan(table2array(tbldataFroi(:,r+5)),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
            close hidden;
            resc = res(2:end,2:end);
            idx = find(cellfun(@isempty,resc));
            resc(idx) = {[0]};
            tblres = array2table(cell2mat(reshape(resc,4,6)),'RowNames',res(2:end,1),'VariableNames',res(1,2:end));
            resFroi{1,x} = tblres;
            x = x + 1;
        end
        clear x r res idx resc tblres
        % delete(findall(0));

        labelfactors = resFroi{1,1}.Properties.RowNames(1:end-2,1);
        for f = 1:numel(labelfactors)
            n_sig = sum(pvalFroi(f,:) <= p);
            fprintf('%d %s effects are significant p<=%.2f without correction for Froi\n',n_sig,labelfactors{f,1}, p)
        end
        clear r x res n_sig

        %%fdr%%%%%%%%%%%%
        if fdrmode == 1
            FDR_CO(pvalFroi, labelfactors, p, dat, savepath)
        end

        %%tableaux des résultats%%%%%%%%%%%%%%%%
        [~,~,~,~,~,~] = Tablegraph_CO(resFroi, pvalFroi, labelfactors, p, dat, lgndroi, ...
            labelroiALL{1,R}, MATroiG1G2{1,R}, datameanroiG1{1,R}, datameanroiG2{1,R}, ...
            roi{1,R}, fileorderconnectogram, savepath, graphmode);

        %%AROI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%ANOVA%%%%%%%%%%%%%%%%%%%
        R = 4;
        dat = lgndroi{1,R};
        x = 1;
        resAroi = [];
        pvalAroi = [];
        for r = 1:(width(tbldataAroi)- 5)
            [pvalAroi(:,x),res] = anovan(table2array(tbldataAroi(:,r+5)),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
            close hidden;
            resc = res(2:end,2:end);
            idx = find(cellfun(@isempty,resc));
            resc(idx) = {[0]};
            tblres = array2table(cell2mat(reshape(resc,4,6)),'RowNames',res(2:end,1),'VariableNames',res(1,2:end));
            resAroi{1,x} = tblres;
            x = x + 1;
        end
        clear x r res idx resc tblres
        % delete(findall(0));

        labelfactors = resAroi{1,1}.Properties.RowNames(1:end-2,1);
        for f = 1:numel(labelfactors)
            n_sig = sum(pvalAroi(f,:) <= p);
            fprintf('%d %s effects are significant p<=%.2f without correction for Aroi\n',n_sig,labelfactors{f,1}, p)
        end

        clear r x res n_sig

        %%%%fdr%%%%%%%%%%%%%
        if fdrmode == 1
            FDR_CO(pvalAroi, labelfactors, p, dat, savepath)
        end

        %%tableaux des résultats%%%%%%%%%%%%%%%%
        [~,~,~,~,~,~] = Tablegraph_CO(resAroi, pvalAroi, labelfactors, p, dat, lgndroi, ...
            labelroiALL{1,R}, MATroiG1G2{1,R}, datameanroiG1{1,R}, datameanroiG2{1,R}, ...
            roi{1,R}, fileorderconnectogram, savepath, graphmode);

        %         if graphmode == 1
        %             tblgrpsigALL = [tblgrpsigMroi, tblgrpsigRroi, tblgrpsigFroi, tblgrpsigAroi];
        %             writetable(tblgrpsigALL,[savepath date '_grpsigresultsALL.xls'],'WriteRowNames',true);
        %         end
    end
    X = ['Results saved in ', savepath];
    disp(X)
    clear X
end

clearvars -except datapath savepathinitial channelmode importedROImode calculatedROImode graphmode fileorderconnectogram fdrmode nointeractionmode interactionmode p
load ([datapath 'workspace.mat'])
load ([datapath 'workspacemat.mat'])
%savepathinitial='C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Stats\CORR0,01_0,08\PhysioSatRespEKG\ANCOVA\';
if ~isfolder(savepathinitial)
    mkdir(savepathinitial)
end

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
        dat = lgndch{1,:};
        x = 1;
        resch = [];
        pvalch = [];
        for r = 1:(width(tblch)- 5)
            [pvalch(:,x),res] = anovan(table2array(tblch(:,r+5)),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
            close hidden;
            resc = res(2:end,2:end);
            idx = find(cellfun(@isempty,resc));
            resc(idx) = {[0]};
            tblres = array2table(cell2mat(reshape(resc,5,6)),'RowNames',res(2:end,1),'VariableNames',res(1,2:end));
            resch{1,x} = tblres;
            x = x + 1;
        end

        clear x r res idx resc tblres
        % delete(findall(0));

        labelfactors = resch{1,1}.Properties.RowNames(1:end-2,1);
        for f = 1:numel(labelfactors)
            n_sig = sum(pvalch(f,:) <= p);
            fprintf('%d %s effects are significant p<=%.2f without correction for channels\n',n_sig,labelfactors{f,1}, p)
        end

        clear r x res n_sig

        %%fdr%%%%%%%%%%%%
        if fdrmode == 1
            FDR_CO(pvalch, labelfactors, p, dat, savepath)
        end

        %%tableaux des résultats%%%%%%%%%%%%%%%%
        [~,~,~,~,~,~] = Tablegraph_CO(resch, pvalch, labelfactors, p, dat, lgndch, ...
            labelch, MATallG1G2, datameanchG1, datameanchG2, DATA, fileorderconnectogram, ...
            savepath, graphmode, ListChinv);
    end

    if importedROImode
        dat = lgndch{1,:};
        R = 1;
        %%ANOVA%%%%%%%%%%%%%%%%%%%
        x = 1;
        resroi = [];
        for r = 1:(width(tblroi)- 5)
            [pvalroi(:,x),res] = anovan(table2array(tblroi(:,r+5)),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
            close hidden;
            resroi = [resroi, res];
            x = x + 1;
        end
        % delete(findall(0));

        labelfactors = resroi(2:end-2,1);
        for f = 1:numel(labelfactors)
            n_sig = sum(pvalroi(f,:) <= p);
            fprintf('%d %s effects are significant p<=%.2f without correction for ROIs\n',n_sig,labelfactors{f,1}, p)
        end

        clear r x res n_sig

        %%fdr%%%%%%%%%%%%
        if fdrmode == 1
            FDR_CO(pvalroi, labelfactors, p, dat, savepath)
        end

        %%tableaux des résultats%%%%%%%%%%%%%%%%
        if graphmode == 1
            %[~,~,~,~,~,~] = Tablegraph_CO(resroi, pvalroi, p, dat, lgndroi, labelroiALL{1,R}, MATroiG1G2{1,R}, datameanroiG1{1,R}, datameanroiG2{1,R}, roi{1,R}, fileorderconnectogram, savepath);

            tblresroi = array2table(resroi);
            for r = 1:(width(tblroi)- 5)
                tblresroi = mergevars(tblresroi,(r:r+6));
            end
            tblresroi.Properties.VariableNames = labelroiALL;
            tblresroi.Properties.RowNames = {'Source','Group','SES','Group*SES','Error','Total'};

            tblproi = array2table(pvalroi);
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
                    MATproi(r,rr)= pvalroi(1,idx); %Mettre les  p sous forme de matrice
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
        dat = lgndroi{1,R};
        x = 1;
        resMroi = [];
        pvalMroi = [];
        for r = 1:(width(tbldataMroi)- 5)
            [pvalMroi(:,x),res] = anovan(table2array(tbldataMroi(:,r+5)),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
            close hidden;
            resc = res(2:end,2:end);
            idx = find(cellfun(@isempty,resc));
            resc(idx) = {[0]};
            tblres = array2table(cell2mat(reshape(resc,5,6)),'RowNames',res(2:end,1),'VariableNames',res(1,2:end));
            resMroi{1,x} = tblres; %[resMroi, tblres];
            x = x + 1;
        end

        clear x r res idx resc tblres
        % delete(findall(0));

        labelfactors = resMroi{1,1}.Properties.RowNames(1:end-2,1);
        for f = 1:numel(labelfactors)
            n_sig = sum(pvalMroi(f,:) <= p);
            fprintf('%d %s effects are significant p<=%.2f without correction for Mroi\n',n_sig,labelfactors{f,1}, p)
        end

        clear r x res n_sig

        %%fdr%%%%%%%%%%%%
        if fdrmode == 1
            FDR_CO(pvalMroi, labelfactors, p, dat, savepath)
        end

        %%tableaux des résultats%%%%%%%%%%%%%%%%
        [~,~,~,~,~,~] = Tablegraph_CO(resMroi, pvalMroi, labelfactors, p, dat, lgndroi, ...
            labelroiALL{1,R}, MATroiG1G2{1,R}, datameanroiG1{1,R}, datameanroiG2{1,R}, ...
            roi{1,R}, fileorderconnectogram, savepath, graphmode);

        %%RROI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        R = 2;
        dat = lgndroi{1,R};
        x = 1;
        resRroi = [];
        pvalRroi = [];
        for r = 1:(width(tbldataRroi)- 5)
            [pvalRroi(:,x),res] = anovan(table2array(tbldataRroi(:,r+5)),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
            close hidden;
            resc = res(2:end,2:end);
            idx = find(cellfun(@isempty,resc));
            resc(idx) = {[0]};
            tblres = array2table(cell2mat(reshape(resc,5,6)),'RowNames',res(2:end,1),'VariableNames',res(1,2:end));
            resRroi{1,x} = tblres;
            x = x + 1;
        end

        clear x r res idx resc tblres
        % delete(findall(0));

        labelfactors = resRroi{1,1}.Properties.RowNames(1:end-2,1);
        for f = 1:numel(labelfactors)
            n_sig = sum(pvalRroi(f,:) <= p);
            fprintf('%d %s effects are significant p<=%.2f without correction for Rroi\n',n_sig,labelfactors{f,1}, p)
        end

        clear r x res n_sig

        %%fdr%%%%%%%%%%%%
        if fdrmode == 1
            FDR_CO(pvalRroi, labelfactors, p, dat, savepath)
        end

        %%tableaux des résultats%%%%%%%%%%%%%%%%
        [~,~,~,~,~,~] = Tablegraph_CO(resRroi, pvalRroi, labelfactors, p, dat, lgndroi, ...
            labelroiALL{1,R}, MATroiG1G2{1,R}, datameanroiG1{1,R}, datameanroiG2{1,R}, ...
            roi{1,R}, fileorderconnectogram, savepath, graphmode);

        %%FROI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        R = 3;
        dat = lgndroi{1,R};
        x = 1;
        resFroi = [];
        pvalFroi = [];
        for r = 1:(width(tbldataFroi)- 5)
            [pvalFroi(:,x),res] = anovan(table2array(tbldataFroi(:,r+5)),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
            close hidden;
            resc = res(2:end,2:end);
            idx = find(cellfun(@isempty,resc));
            resc(idx) = {[0]};
            tblres = array2table(cell2mat(reshape(resc,5,6)),'RowNames',res(2:end,1),'VariableNames',res(1,2:end));
            resFroi{1,x} = tblres;
            x = x + 1;
        end

        clear x r res idx resc tblres
        % delete(findall(0));

        labelfactors = resFroi{1,1}.Properties.RowNames(1:end-2,1);
        for f = 1:numel(labelfactors)
            n_sig = sum(pvalFroi(f,:) <= p);
            fprintf('%d %s effects are significant p<=%.2f without correction for Froi\n',n_sig,labelfactors{f,1}, p)
        end

        clear r x res n_sig

        %%fdr%%%%%%%%%%%%
        if fdrmode == 1
            FDR_CO(pvalFroi, labelfactors, p, dat, savepath)
        end

        %%tableaux des résultats%%%%%%%%%%%%%%%%
        [~,~,~,~,~,~] = Tablegraph_CO(resFroi, pvalFroi, labelfactors, p, dat, lgndroi, ...
            labelroiALL{1,R}, MATroiG1G2{1,R}, datameanroiG1{1,R}, datameanroiG2{1,R}, ...
            roi{1,R}, fileorderconnectogram, savepath, graphmode);

        %%AROI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        R = 4;
        dat = lgndroi{1,R};
        x = 1;
        resAroi = [];
        pvalAroi = [];
        for r = 1:(width(tbldataAroi)- 5)
            [pvalAroi(:,x),res] = anovan(table2array(tbldataAroi(:,r+5)),{gr,ses},'Continuous',2,'model','interaction','varnames',{'Group','SES'});
            close hidden;
            resc = res(2:end,2:end);
            idx = find(cellfun(@isempty,resc));
            resc(idx) = {[0]};
            tblres = array2table(cell2mat(reshape(resc,5,6)),'RowNames',res(2:end,1),'VariableNames',res(1,2:end));
            resAroi{1,x} = tblres;
            x = x + 1;
        end

        clear x r res idx resc tblres
        % delete(findall(0));

        labelfactors = resAroi{1,1}.Properties.RowNames(1:end-2,1);
        for f = 1:numel(labelfactors)
            n_sig = sum(pvalAroi(f,:) <= p);
            fprintf('%d %s effects are significant p<=%.2f without correction for Aroi\n',n_sig,labelfactors{f,1}, p)
        end

        clear r x res n_sig

        %%%%fdr%%%%%%%%%%%%%
        if fdrmode == 1
            FDR_CO(pvalAroi, labelfactors, p, dat, savepath)
        end

        %%tableaux des résultats%%%%%%%%%%%%%%%%
        [~,~,~,~,~,~] = Tablegraph_CO(resAroi, pvalAroi, labelfactors, p, dat, lgndroi, ...
            labelroiALL{1,R}, MATroiG1G2{1,R}, datameanroiG1{1,R}, datameanroiG2{1,R}, ...
            roi{1,R}, fileorderconnectogram, savepath, graphmode);

        %         if graphmode == 1
        %             tblgrpsigALL = [tblgrpsigMroi, tblgrpsigRroi, tblgrpsigFroi, tblgrpsigAroi];
        %             writetable(tblgrpsigALL,[savepath date '_grpsigresultsALL.xls'],'WriteRowNames',true);
        %         end
    end

    X = ['Results saved in ', savepath];
    disp(X)
    clear
end

toc