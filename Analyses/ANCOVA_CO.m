%%%%%%%%%%%%%%ANCOVA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

savepathi = 'C:\data\Malnutrition\Resting\NIRS\Analyses\CORR0,01_0,08\PCAPW\CorrPairC0,1 ExcY\fisher\';
load ([savepathi 'workspace.mat'])
%load ([savepath 'workspacemat.mat'])

disp(['Computing ANCOVA on ' savepathi])

%% PARAMETERS TO MODIFY %%%%%
channelmode = 1;
importedROImode = 0;
calculatedROImode = 1;
graphmode = 1;
fdrmode = 1;
nointeractionmode = 0;
interactionmode = 1;
savemode = 1;

p = 0.05;
factors = {gr,ses};
labelfactors = {'Group','SES'};
labeldatai = {'Pearson Correlation'};

fileorderconnectogram = {'C:\data\Malnutrition\Resting\NIRS\Analyses\Connectogram_Mixte.txt',...
    'C:\data\Malnutrition\Resting\NIRS\Analyses\Connectogram_Region.txt',...
    'C:\data\Malnutrition\Resting\NIRS\Analyses\Connectogram_Fonction.txt',...
    'C:\data\Malnutrition\Resting\NIRS\Analyses\Connectogram_Aire.txt'};
%%END MODIFY %%%%%%%%%%%%%

%% ANCOVA SANS INTERACTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nointeractionmode == 1
    disp('Computing ANCOVA without interaction')

    type = 'main';
    savepath = fullfile(savepathi, ['ANOVA' type],filesep);
    if ~isfolder(savepath)
        mkdir(savepath)
    end

    savepathgraph = fullfile(savepath,'Graphs',filesep);
    if ~isfolder(savepathgraph)
        mkdir(savepathgraph)
    end
    %% CHANNELS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if channelmode ==1

        %%ANOVA%%%%%%%%%%%%%%%%%%%
        labeldim = labelch;
        labeldata = [labeldatai; lgndch{1,:}];
        [tblpvalch, tblresch, slopes] = ANOVA_job(type, table2array(tblch(:,6:end)),...
            factors, labelfactors, p, labeldata, labeldim, fdrmode, graphmode, savepath, savemode);
        %
        %         [~,atab,ctab,stats] = aoctool(ses,table2array(tblch(:,7)),gr,0.05,'SES','Connectivity','Group','on','separate lines');
        %         [c,m] = multcompare(stats,'Estimate','slope');
        %         [c,m] = multcompare(stats,'Estimate','pmm');
        %         [c,m] = multcompare(stats,'Estimate','intercept');

        %%Graph des résultats%%%%%%%%%%%%%%%%
        if graphmode
            ANCOVAGraph_CO(tblpvalch, p, labeldata, MATallG1G2, ...
                DATA, fileorderconnectogram{1,1}, savepathgraph, ListChinv, slopes);
        end
    end

    if importedROImode %%not tested
        %%ANOVA%%%%%%%%%%%%%%%%%%%
        labeldim = labelroi;
        labeldata = [labeldatai; lgndroi{1,:}];
        [tblpvalroi, tblresroi, slopes] = ANOVA_job(type, table2array(tblroi(:,6:end)),...
            factors, labelfactors, p, labeldata, labeldim, fdrmode, graphmode, savepath, savemode);

        %%Graph des résultats%%%%%%%%%%%%%%%%
        if graphmode
            ANCOVAGraph_CO(tblpvalroi, p, labeldata, MATallG1G2, ...
                DATA, fileorderconnectogram{1,1}, savepathgraph, ListChinv, slopes);
        end
    end

    if calculatedROImode
        %% MROI%%%%%%%%%%%%%%%%%%%%%%%
        %ANOVA%%%%%%%%%%%%%%%%%%%
        R = 1;
        labeldim = labelroiALL{1,R};
        labeldata = [labeldatai; lgndroi{1,R}];
        [tblpvalMroi, tblresMroi, slopes] = ANOVA_job(type, table2array(tbldataMroi(:,6:end)),...
            factors, labelfactors, p, labeldata, labeldim, fdrmode, graphmode, savepath, savemode);

        %%Graph des résultats%%%%%%%%%%%%%%%%
        if graphmode
            ANCOVAGraph_CO(tblpvalMroi, p, labeldata, MATroiG1G2{1,R}, ...
                roi{1,R}, fileorderconnectogram{1,R}, savepathgraph, slopes);
        end

        %% RROI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%ANOVA%%%%%%%%%%%%%%%%%%%
        R = 2;
        labeldim = labelroiALL{1,R};
        labeldata = [labeldatai; lgndroi{1,R}];

        [tblpvalRroi, tblresRroi, slopes] = ANOVA_job(type, table2array(tbldataRroi(:,6:end)),...
            factors, labelfactors, p, labeldata, labeldim, fdrmode, graphmode, savepath, savemode);

        %%Graph des résultats%%%%%%%%%%%%%%%%
        if graphmode
            ANCOVAGraph_CO(tblpvalRroi, p, labeldata, MATroiG1G2{1,R}, ...
                roi{1,R}, fileorderconnectogram{1,R}, savepathgraph, slopes);
        end

        %% FROI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%ANOVA%%%%%%%%%%%%%%%%%%%
        R = 3;
        labeldim = labelroiALL{1,R};
        labeldata = [labeldatai; lgndroi{1,R}];

        [tblpvalFroi, tblresFroi, slopes] = ANOVA_job(type, table2array(tbldataFroi(:,6:end)),...
            factors, labelfactors, p, labeldata, labeldim, fdrmode, graphmode, savepath, savemode);

        %%Graph des résultats%%%%%%%%%%%%%%%%
        if graphmode
            ANCOVAGraph_CO(tblpvalFroi, p, labeldata, MATroiG1G2{1,R}, ...
                roi{1,R}, fileorderconnectogram{1,R}, savepathgraph, slopes);
        end

        %% AROI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%ANOVA%%%%%%%%%%%%%%%%%%%
        R = 4;
        labeldim = labelroiALL{1,R};
        labeldata = [labeldatai; lgndroi{1,R}];

        [tblpvalAroi, tblresAroi, slopes] = ANOVA_job(type, table2array(tbldataAroi(:,6:end)),...
            factors, labelfactors, p, labeldata, labeldim, fdrmode, graphmode, savepath, savemode);

        %%Graph des résultats%%%%%%%%%%%%%%%%
        if graphmode
            ANCOVAGraph_CO(tblpvalAroi, p, labeldata, MATroiG1G2{1,R}, ...
                roi{1,R}, fileorderconnectogram{1,R}, savepathgraph, slopes);
        end

        %         if graphmode == 1
        %             tblgrpsigALL = [tblgrpsigMroi, tblgrpsigRroi, tblgrpsigFroi, tblgrpsigAroi];
        %             writetable(tblgrpsigALL,[savepath date '_grpsigresultsALL.xls'],'WriteRowNames',true);
        %         end
    end
    fprintf('Results saved in %s\n', savepath)

    %clearvars -except factors labelfactors labeldatai savepathi channelmode importedROImode calculatedROImode graphmode fileorderconnectogram fdrmode nointeractionmode interactionmode p
    %load ([savepathi 'workspace.mat'])
    %load ([savepathi 'workspacemat.mat'])

end

%% ANCOVA AVEC INTERACTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if interactionmode == 1
    disp('Computing ANCOVA with interaction')

    type = 'inter';
    savepath = fullfile(savepathi, ['ANOVA' type],filesep);
    if ~isfolder(savepath)
        mkdir(savepath)
    end

    savepathgraph = fullfile(savepath,'Graphs',filesep);
    if ~isfolder(savepathgraph)
        mkdir(savepathgraph)
    end
    %% CHANNELS%%%%%%%%%%%%%%%%%%%%%
    if channelmode ==1
        labeldim = labelch; %'1-2','1-3','1-4'...
        labeldata = [labeldatai; lgndch{1,:}]; %'Ch'
        [tblpvalch, tblresch, slopes] = ANOVA_job(type, table2array(tblch(:,6:end)),...
            factors, labelfactors, p, labeldata, labeldim, fdrmode, graphmode, savepath, savemode);

        %%Graph des résultats%%%%%%%%%%%%%%%%
        if graphmode
            ANCOVAGraph_CO(tblpvalch, p, labeldata, MATallG1G2, ...
                DATA, fileorderconnectogram{1,1}, savepathgraph, ListChinv, slopes);
        end
    end

    if importedROImode %%not tested
        %%ANOVA%%%%%%%%%%%%%%%%%%%
        labeldim = labelroi;
        labeldata = [labeldatai; lgndroi{1,:}];
        [tblpvalroi, tblresroi, slopes] = ANOVA_job(type, table2array(tblroi(:,6:end)), factors,...
            labelfactors, p, labeldata, labeldim, fdrmode, graphmode, savepath, savemode);

        %%Graph des résultats%%%%%%%%%%%%%%%%
        if graphmode
            ANCOVAGraph_CO(tblpvalroi, p, labeldata, MATallG1G2, ...
                DATA, fileorderconnectogram{1,1}, savepathgraph, ListChinv, slopes);
        end
    end

    if calculatedROImode

        %% MROI%%%%%%%%%%%%%%%%%%%%%%%
        R = 1;
        labeldim = labelroiALL{1,R};
        labeldata = [labeldatai; lgndroi{1,R}];

        [tblpvalMroi, tblresMroi, slopes] = ANOVA_job(type, table2array(tbldataMroi(:,6:end)),...
            factors, labelfactors, p, labeldata, labeldim, fdrmode, graphmode, savepath, savemode);

        %%Graph des résultats%%%%%%%%%%%%%%%%
        if graphmode
            ANCOVAGraph_CO(tblpvalMroi, p, labeldata, MATroiG1G2{1,R}, ...
                roi{1,R}, fileorderconnectogram{1,R}, savepathgraph, slopes);
        end

        %% RROI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        R = 2;
        labeldim = labelroiALL{1,R};
        labeldata = [labeldatai; lgndroi{1,R}];

        [tblpvalRroi, tblresRroi, slopes] = ANOVA_job(type, table2array(tbldataRroi(:,6:end)),...
            factors, labelfactors, p, labeldata, labeldim, fdrmode, graphmode, savepath, savemode);

        %%Graph des résultats%%%%%%%%%%%%%%%%
        if graphmode
            ANCOVAGraph_CO(tblpvalRroi, p, labeldata, MATroiG1G2{1,R}, ...
                roi{1,R}, fileorderconnectogram{1,R}, savepathgraph, slopes);
        end

        %% FROI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        R = 3;
        labeldim = labelroiALL{1,R};
        labeldata = [labeldatai; lgndroi{1,R}];

        [tblpvalFroi, tblresFroi, slopes] = ANOVA_job(type, table2array(tbldataFroi(:,6:end)),...
            factors, labelfactors, p, labeldata, labeldim, fdrmode, graphmode, savepath, savemode);

        %%Graph des résultats%%%%%%%%%%%%%%%%
        if graphmode
            ANCOVAGraph_CO(tblpvalFroi, p, labeldata, MATroiG1G2{1,R}, ...
                roi{1,R}, fileorderconnectogram{1,R}, savepathgraph, slopes);
        end

        %% AROI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        R = 4;
        labeldim = labelroiALL{1,R};
        labeldata = [labeldatai; lgndroi{1,R}];

        [tblpvalAroi, tblresAroi, slopes] = ANOVA_job(type, table2array(tbldataAroi(:,6:end)),...
            factors, labelfactors, p, labeldata, labeldim, fdrmode, graphmode, savepath, savemode);

        %%Graph des résultats%%%%%%%%%%%%%%%%
        if graphmode
            ANCOVAGraph_CO(tblpvalAroi, p, labeldata, MATroiG1G2{1,R}, ...
                roi{1,R}, fileorderconnectogram{1,R}, savepathgraph, slopes);
        end

        %         if graphmode == 1
        %             tblgrpsigALL = [tblgrpsigMroi, tblgrpsigRroi, tblgrpsigFroi, tblgrpsigAroi];
        %             writetable(tblgrpsigALL,[savepath date '_grpsigresultsALL.xls'],'WriteRowNames',true);
        %         end
    end
    fprintf('Results saved in %s\n', savepath)
end
clear
toc

%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ANCOVAGraph_CO(tblpval, p, labeldata, MATmeandiff, ...
    DATA, fileorderconnectogram, savepath, varargin)

if ~isempty(varargin)
    if numel(varargin) == 2
        slopes = varargin{2};
    elseif numel(varargin) == 1
        slopes = varargin{1};
    end
end

labelrows = tblpval.Row;
sig = tblpval{:,:} <= p;
if any(sig,'all')

    for f = 1:numel(labelrows)
        fac = labelrows{f,1};
        if contains(fac,'*')
            fac = erase(fac,'*');
        elseif contains(fac,':')
            fac = erase(fac,':');
        end

        sig = tblpval{f,:} <= p;
        if any(sig,'all')

            %%%% Connectogramme%%%%%%%%%%%%
            if strcmp(fac,'Group')
                MATpval = zeros(length(MATmeandiff));
                x = 1;
                for c = 1:length(MATmeandiff)
                    for cc =(c + 1):length(MATmeandiff)
                        MATpval(c,cc) = tblpval{f,x}; %Mettre les  p sous forme de matrice
                        x = x + 1;
                    end
                end

                MATpval(length(MATpval),length(MATpval)) = 0; %Mettre la dernière donnée sur la diagonale à 0
                MATpval = MATpval + triu(MATpval,1)'; %rendre la matrice symétrique
                tf = MATpval <=p & MATpval > 0; %trouver les p significatif
                MATsig = MATmeandiff;
                MATsig(~tf) = 0; %extraire juste les connexions significatives
                clear c cc x tf

                sign = {'Neg'; 'Pos'};
                for s = 1:numel(sign)
                    if strcmp('Neg',sign{s})
                        MATneg = MATsig < 0; %negative differences
                        MAT = MATsig.*(MATneg); %DATA{id}.MAT %loader matrice
                    elseif strcmp('Pos',sign{s})
                        MATpos = MATsig > 0; %positive differences
                        MAT = MATsig.*(MATpos);%DATA{id}.MAT %loader matrice
                    end
                    id = 1;
                    if any(MAT,'all') && contains(labeldata{2,1},'Ch')
                        List = strvcat(DATA{id}.ZoneList); %liste des paires SD
                        ML = DATA{id}.zone.ml; %Loader S/D/ROI/Gr
                        plotLst = DATA{id}.zone.plotLst;
                        labelzone =  DATA{id}.zone.label;
                        labelzone = strrep(labelzone, '_', ' ');
                        plotconnectogram(fileorderconnectogram,MAT(varargin{1,1},varargin{1,1}),List,labelzone,plotLst,ML)
                        fig = gcf;
                        %pbaspect([1 1 1])
                        fig.WindowState = 'maximized';
                        colorbar('Position',[0.663 0.334 0.02 0.4],'FontSize',12)
                        savefig([savepath fac 'effect_' labeldata{2,1} '_' sign{s} 'ConnectG1G2']);
                        exportgraphics(gcf,[savepath fac 'effect_' labeldata{2,1} '_' sign{s} 'ConnectG1G2.png'])
                        close
                    elseif any(MAT,'all') && contains(labeldata{2,1},'roi')
                        plotLst = DATA(2,:);
                        labelzone =  DATA(1,:);
                        plotconnectogramroi(fileorderconnectogram,MAT,labelzone,plotLst)
                        fig = gcf;
                        %pbaspect([1 1 1])
                        fig.WindowState = 'maximized';
                        if contains(labeldata{2,1},'roiR')
                            colorbar('Position',[0.75 0.334 0.02 0.4],'FontSize',12)
                        else
                            colorbar('Position',[0.663 0.334 0.02 0.4],'FontSize',12)
                        end
                        savefig([savepath fac 'effect_' labeldata{2,1} '_' sign{s} 'ConnectG1G2']);
                        exportgraphics(gcf,[savepath fac 'effect_' labeldata{2,1} '_' sign{s} 'ConnectG1G2.png'])
                        close
                    else
                    end
                    clear MAT
                end

            elseif strcmp(fac,'SES')
                idx = find(sig);
                Vslopes = zeros(1,length(sig));
                Vslopes(idx) = slopes.(fac);
                MATslopes = zeros(length(MATmeandiff));
                x = 1;
                for c = 1:length(MATmeandiff)
                    for cc =(c + 1):length(MATmeandiff)
                        MATslopes(c,cc) = Vslopes(x); %Mettre les  p sous forme de matrice
                        x = x + 1;
                    end
                end

                MATslopes(length(MATslopes),length(MATslopes)) = 0; %Mettre la dernière donnée sur la diagonale à 0
                MATslopes = MATslopes + triu(MATslopes,1)'; %rendre la matrice symétrique
                clear c cc x

                if any(MATslopes,'all') && contains(labeldata{2,1},'Ch')
                    id = 1;
                    List = strvcat(DATA{id}.ZoneList); %liste des paires SD
                    ML = DATA{id}.zone.ml; %Loader S/D/ROI/Gr
                    plotLst = DATA{id}.zone.plotLst;
                    labelzone =  DATA{id}.zone.label;
                    labelzone = strrep(labelzone, '_', ' ');
                    plotconnectogram(fileorderconnectogram,MATslopes(varargin{1,1},varargin{1,1}),List,labelzone,plotLst,ML)
                    fig = gcf;
                    %pbaspect([1 1 1])
                    fig.WindowState = 'maximized';
                    colorbar('Position',[0.663 0.334 0.02 0.4],'FontSize',12)
                    savefig([savepath fac 'effect_' labeldata{2,1} '_Connect']);
                    exportgraphics(gcf,[savepath fac 'effect_' labeldata{2,1} '_Connect.png'])
                    close
                elseif any(MATslopes,'all') && contains(labeldata{2,1},'roi')
                    plotLst = DATA(2,:);
                    labelzone =  DATA(1,:);
                    plotconnectogramroi(fileorderconnectogram,MATslopes,labelzone,plotLst)
                    fig = gcf;
                    %pbaspect([1 1 1])
                    fig.WindowState = 'maximized';
                    if contains(labeldata{2,1},'roiR')
                        colorbar('Position',[0.75 0.334 0.02 0.4],'FontSize',12)
                    else
                        colorbar('Position',[0.663 0.334 0.02 0.4],'FontSize',12)
                    end
                    savefig([savepath fac 'effect_' labeldata{2,1} '_Connect']);
                    exportgraphics(gcf,[savepath fac 'effect_' labeldata{2,1} '_Connect.png'])
                    close
                else
                end
                clear MATslopes Vslopes

            elseif strcmp(fac,'GroupSES')
                for g = 1:size(slopes.(fac),1)
                    idx = find(sig);
                    Vslopes = zeros(1,length(sig));
                    Vslopes(idx) = slopes.(fac)(g,:);
                    MATslopes = zeros(length(MATmeandiff));
                    x = 1;
                    for c = 1:length(MATmeandiff)
                        for cc =(c + 1):length(MATmeandiff)
                            MATslopes(c,cc) = Vslopes(x); %Mettre les  p sous forme de matrice
                            x = x + 1;
                        end
                    end

                    MATslopes(length(MATslopes),length(MATslopes)) = 0; %Mettre la dernière donnée sur la diagonale à 0
                    MATslopes = MATslopes + triu(MATslopes,1)'; %rendre la matrice symétrique
                    clear c cc x

                    if any(MATslopes,'all') && contains(labeldata{2,1},'Ch')
                        id = 1;
                        List = strvcat(DATA{id}.ZoneList); %liste des paires SD
                        ML = DATA{id}.zone.ml; %Loader S/D/ROI/Gr
                        plotLst = DATA{id}.zone.plotLst;
                        labelzone =  DATA{id}.zone.label;
                        labelzone = strrep(labelzone, '_', ' ');
                        plotconnectogram(fileorderconnectogram,MATslopes(varargin{1,1},varargin{1,1}),List,labelzone,plotLst,ML)
                        fig = gcf;
                        %pbaspect([1 1 1])
                        fig.WindowState = 'maximized';
                        colorbar('Position',[0.663 0.334 0.02 0.4],'FontSize',12)
                        savefig([savepath fac 'effect_' labeldata{2,1} '_ConnectG' num2str(g)]);
                        exportgraphics(gcf,[savepath fac 'effect_' labeldata{2,1} '_ConnectG' num2str(g) '.png'])
                        close
                    elseif any(MATslopes,'all') && contains(labeldata{2,1},'roi')
                        plotLst = DATA(2,:);
                        labelzone =  DATA(1,:);
                        plotconnectogramroi(fileorderconnectogram,MATslopes,labelzone,plotLst)
                        fig = gcf;
                        %pbaspect([1 1 1])
                        fig.WindowState = 'maximized';
                        if contains(labeldata{2,1},'roiR')
                            colorbar('Position',[0.75 0.334 0.02 0.4],'FontSize',12)
                        else
                            colorbar('Position',[0.663 0.334 0.02 0.4],'FontSize',12)
                        end
                        savefig([savepath fac 'effect_' labeldata{2,1} '_ConnectG' num2str(g)]);
                        exportgraphics(gcf,[savepath fac 'effect_' labeldata{2,1} '_ConnectG' num2str(g) '.png'])
                        close
                    else
                    end
                    clear MATslopes Vslopes
                end
            end
        end
    end
end
end
