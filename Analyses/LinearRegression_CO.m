%%%%%%%%%%%%%%ANCOVA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
disp('Computing Linear Regression')

savepathi = 'C:\data\Malnutrition\Resting\NIRS\Analyses\CORR0,01_0,08\PCAPW\CorrPairC0,1 ExcY\fisher\';
load ([savepathi 'workspace.mat'])

disp(['Computing Linear Regression on ' savepathi])

%% PARAMETERS TO MODIFY %%%%%
channelmode = 1;
importedROImode = 0;
calculatedROImode = 1;
graphmode = 1;
%fileorderconnectogram  = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Mixte.txt';
fdrmode = 1;
%nointeractionmode = 1;
%interactionmode = 1;
savemode = 1;

p = 0.05;
factors = {gr,ses};
labelfactors = {'Group','SES'};
labeldatai = {'Pearson Correlation'};
type = 'component';

%%END MODIFY %%%%%%%%%%%%%
savepath = fullfile(savepathi, 'Regression',filesep);
if ~isfolder(savepath)
    mkdir(savepath)
end

% savepathgraph = fullfile(savepath,'Graphs',filesep);
% if ~isfolder(savepathgraph)
%     mkdir(savepathgraph)
% end

%% CHANNELS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if channelmode ==1
    disp('Computing Linear Regression on channels')
    %dat = lgndch{1,:};
    labeldim = labelch;
    labeldata = [labeldatai; lgndch{1,:}];

    %%REGRESSION%%%%%%%%%%%%%%%%%%%
    [tblpvalresch, tblresch, tblcoefch, tblintersigch] = Regression_job(type, table2array(tblch(:,6:end)),...
        factors, labelfactors, p, labeldata, labeldim, fdrmode, graphmode, savepath, savemode);

    %%tableaux des résultats%%%%%%%%%%%%%%%%
    %[~,~,~,~,~,~] = Tablegraph_CO(resch, pvalch, labelfactors, p, labeldim, lgndch, labelch, MATmeanG1G2, meanG1ch, meanG2ch, DATA, fileorderconnectogram, savepath, graphmode, coefch, coefchsiginterG1, coefchsiginterG2, pvalchsiginterG1, pvalchsiginterG2);
end

if calculatedROImode ==1
    disp('Computing Linear Regression on ROIs')

    %% MROI%%%%%%%%%%%%%%%%%%%%%%%
    R = 1;
    labeldim = labelroiALL{1,R};
    labeldata = [labeldatai; lgndroi{1,R}];

    [tblpvalresMroi, tblresMroi, tblcoefMroi, tblintersigMroi] = Regression_job(type, table2array(tbldataMroi(:,6:end)),...
        factors, labelfactors, p, labeldata, labeldim, fdrmode, graphmode, savepath, savemode);


    %%tableaux des résultats%%%%%%%%%%%%%%%%
    %     if tablegraphmode == 1
    %         [~,~,~,~,~,~] = Tablegraph_CO(resMroi, pvalMroi, labelfactors, p, labeldim, lgndroi, labelroiALL{1,R}, MATroimeanG1G2{1,R}, meanroiG1r{1,R}, meanroiG2r{1,R}, roi{1,R}, fileorderconnectogram, savepath, graphmode, coefMroi, coefMroisiginterG1, coefMroisiginterG2, pvalMroisiginterG1, pvalMroisiginterG2);
    %     end

    %% RROI%%%%%%%%%%%%%%%%%%%%%%%
    R = 2;
    labeldim = labelroiALL{1,R};
    labeldata = [labeldatai; lgndroi{1,R}];

    [tblpvalresRroi, tblresRroi, tblcoefRroi, tblintersigRroi] = Regression_job(type, table2array(tbldataRroi(:,6:end)),...
        factors, labelfactors, p, labeldata, labeldim, fdrmode, graphmode, savepath, savemode);

    %%tableaux des résultats%%%%%%%%%%%%%%%%
    %     if tablegraphmode == 1
    %         [~,~,~,~,~,~] = Tablegraph_CO(resRroi, pvalRroi, labelfactors, p, labeldim, lgndroi, labelroiALL{1,R}, MATroimeanG1G2{1,R}, meanroiG1r{1,R}, meanroiG2r{1,R}, roi{1,R}, fileorderconnectogram, savepath, graphmode, coefRroi, coefRroisiginterG1, coefRroisiginterG2, pvalRroisiginterG1, pvalRroisiginterG2);
    %     end

    %% FROI%%%%%%%%%%%%%%%%%%%%%%%
    R = 3;
    labeldim = labelroiALL{1,R};
    labeldata = [labeldatai; lgndroi{1,R}];

    [tblpvalresFroi, tblresFroi, tblcoefFroi, tblintersigFroi] = Regression_job(type, table2array(tbldataFroi(:,6:end)),...
        factors, labelfactors, p, labeldata, labeldim, fdrmode, graphmode, savepath, savemode);

    %%tableaux des résultats%%%%%%%%%%%%%%%%
    %     if tablegraphmode == 1
    %         [~,~,~,~,~,~] = Tablegraph_CO(resFroi, pvalFroi, labelfactors, p, labeldim, lgndroi, labelroiALL{1,R}, MATroimeanG1G2{1,R}, meanroiG1r{1,R}, meanroiG2r{1,R}, roi{1,R}, fileorderconnectogram, savepath, graphmode, coefFroi, coefFroisiginterG1, coefFroisiginterG2, pvalFroisiginterG1, pvalFroisiginterG2);
    %     end

    %% AROI%%%%%%%%%%%%%%%%%%%%%%%
    R = 4;
    labeldim = labelroiALL{1,R};
    labeldata = [labeldatai; lgndroi{1,R}];

    [tblpvalresAroi, tblresAroi, tblcoefAroi, tblintersigAroi] = Regression_job(type, table2array(tbldataAroi(:,6:end)),...
        factors, labelfactors, p, labeldata, labeldim, fdrmode, graphmode, savepath, savemode);

    %%tableaux des résultats%%%%%%%%%%%%%%%%
    %     if tablegraphmode == 1
    %         [~,~,~,~,~,~] = Tablegraph_CO(resAroi, pvalAroi, labelfactors, p, labeldim, lgndroi, labelroiALL{1,R}, MATroimeanG1G2{1,R}, meanroiG1r{1,R}, meanroiG2r{1,R}, roi{1,R}, fileorderconnectogram, savepath, graphmode, coefAroi, coefAroisiginterG1, coefAroisiginterG2, pvalAroisiginterG1, pvalAroisiginterG2);
    %     end
end

X = ['Results saved in ', savepath];
disp(X)
clear

toc