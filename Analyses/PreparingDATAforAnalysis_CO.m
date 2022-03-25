%%%%%%%%%%%%%%%%%PREPARING CORR DATA FOR ANALYSIS%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
disp(['Computing PreparingDATAforAnalysis'])
savepath='C:\data\Malnutrition\Resting\NIRS\Analyses\Stats\CORR0,01_0,08\PCAPW_nocriteria\';

fisher = 1;
if fisher == 0
    savepath = [savepath 'nofisher\'];
elseif fisher == 1
    savepath = [savepath 'fisher\'];
end
if ~isfolder(savepath)
    mkdir(savepath)
end

importROI = 0; %0=import channels; 1=importROI%
calculateROI = 1; %Calculer la connectivité moyennée par ROI (seulement possible si channels importés)
channelmode = 1; %Faire les analyses sur les canaux
ROImode = 1; %Faire les analyses sur les ROI

%connectivity = 'CORR'; %Modify 'COH' OR 'CORR'
xlslistfile = 'C:\data\Malnutrition\Resting\NIRS\Analyses\CORRmatrice0,01_0,08\Channels\PCAPW_nocriteria\Subjectlist N=54.xlsx'; %Fichier excel avec le dossier des matrices, leur nom et le groupe
exceltable = 'C:\data\Malnutrition\Resting\NIRS\participants list.xlsx'; %%%% Fichier excel avec les données démographiques d'intérêt

%%%%%% from StatMatrices of LIONIRS toolbox%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Loading the data from ' fileparts(xlslistfile)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Load the excel file in array info%%

[filepath,name,ext] = fileparts([xlslistfile]);
if strcmp(ext,'.xlsx') | strcmp(ext,'.xls')
    [raw, txt, info] = xlsread([xlslistfile]);
elseif strcmp(ext,'.txt')
    [num, txt, info] = readtxtfile_asxlsread([xlslistfile]);
end

%Load the matrices
groupeall = [];
for isubject = 2:size(info,1)
    id = isubject-1;

    %Load the matrix and copy onto DATA
    MAT = load(fullfile(info{isubject,1},[info{isubject,2},'.mat']));
    if isfield(MAT, 'meancorr')
        matcorr = MAT.meancorr;
        if fisher == 1 % add fisher transform
            matcorr =  1/2*log((1+matcorr)./(1-matcorr));
            DATA{id}.MAT = matcorr;
        else
            DATA{id}.MAT = matcorr;
        end
    end
    DATA{id}.name = info{isubject,2};
    DATA{id}.GR = info{isubject,4};

    %Import ZoneList (from the CO data)
    DATA{id}.ZoneList = MAT.ZoneList;

    %création liste des noms de sujets et liste du groupe de chacun%
    list_subject{id} = DATA{id}.name;
    groupeall = [groupeall; info{isubject,4}];

    %Load the zone file (from excel sheet) and copy onto DATA
    load(fullfile(info{isubject,1}, info{isubject,3}),'-mat');
    names = fieldnames(zone);
    for iname = 1:numel(names)
        eval(['DATA{id}.zone.',names{iname},' =zone.',names{iname},';']);
    end
end

%création d'une matrice avec les matrices de tous les sujets et d'une liste du groupe de chacun%
idsubject = 1:numel(groupeall); %channel mode
for isubject = 1:numel(groupeall)
    MATall(:,:,isubject) = DATA{isubject}.MAT;
    groupid(isubject) = DATA{idsubject(isubject)}.GR;
end

%sélectionner le nom du fichier zone et la liste des canaux du dernier participant%
ZONEid = [info{end,3}];
ZoneList =  DATA{end}.ZoneList;
zone = DATA{end}.zone;
labelnode = 'c';
if numel(ZoneList)>1
    if contains(ZoneList{1},'D0')
        Devicename = 'NIRx';
    else
        Devicename = 'ISS';
    end
else
    Devicename = 'NIRx';
end

%Extract channel list from zone file (from excel sheet)
idzone =[];
ChList = {};
for izone = 1:numel(zone.label)
    idch = zone.plot{izone};
    idzone = [idzone,izone, zeros(1,size(idch,1)-1)];
    for ichzone = 1:size(idch,1)
        ich = idch(ichzone,:);
        switch  Devicename %find label according to system
            case 'ISS'
                strDet = SDDet2strboxy_ISS(ich(2));
                strSrs = SDPairs2strboxy_ISS(ich(1));
                ChList = [ChList;{[strDet, ' ',strSrs ]}];
            case 'NIRx'
                strDet = SDDet2strboxy(ich(2));
                strSrs = SDPairs2strboxy(ich(1));
                ChList = [ChList;{[strDet, ' ',strSrs ]}];
            otherwise
                strDet = SDDet2strboxy(ich(2));
                strSrs = SDPairs2strboxy(ich(1));
                ChList = [ChList;{[strDet, ' ',strSrs ]}];
        end
    end
end

%compare ZoneList (from CO data) to channel list (from excel file) to reorder the data in correct order
for i = 1:numel(ChList)
    ListCh(i,:) = find(strcmp(ZoneList(:),ChList{i}));
end

for r = 1:numel(ListCh)
    ListChinv(ListCh(r,:),:)  = r;
end

% ZoneListold = ZoneList;
% ZoneList = ChList;
MATallold = MATall;
MATall = MATall(ListCh,ListCh,:);
save([savepath date '_MATall.mat'],'MATall');
%save([savepath date '_DATA.mat'],'MATall');

idG1 = find(groupeall==1);
MATallG1 = MATall(:,:,idG1);
MATallmeanG1 = squeeze(nanmean(MATallG1,3));

idG2 = find(groupeall==2);
MATallG2 = MATall(:,:,idG2);
MATallmeanG2 = squeeze(nanmean(MATallG2,3));

MATallG1G2 = MATallmeanG1 - MATallmeanG2;
%MATallG2G1 = MATallmeanG2 - MATallmeanG1;

% %%%%%%%% Calculate ROI from StatMatrix%%%%%%%%%%%%%%%%
% if calculateROI == 1 %création d'une matrice pour les zones de tous les sujets
%     MATallroi = zeros(numel(DATA{id}.zone.label),numel(DATA{id}.zone.label),numel(DATA));
%     for isubject = 1:numel(groupeall) %Pour chaque sujet
%         ZoneList = DATA{isubject}.ZoneList; %Extraire la ZoneList
%         if numel(ZoneList)>1 %modKR
%             if contains(ZoneList(1),'D0')
%                 DATA{isubject}.System = 'NIRx';
%             else
%                 DATA{isubject}.System = 'ISS';
%             end
%         else
%             DATA{isubject}.System = 'NIRx';
%         end
%         for izone = 1:numel(DATA{isubject}.zone.label) %Pour chaque zone
%             ML = DATA{isubject}.zone.ml; %Extraire le ML (paire d'optodes)
%             DATA{isubject}.zone.plotLst; %Extraire le plotLst (id canal)
%             idlisti = [];
%             idliststr = [];
%             chzone = DATA{isubject}.zone.plotLst{izone}; %Extraire les canaux de la zone
% 
%             for ichzone = 1:numel(chzone); %pour chaque canal dans la zone
%                 ich = chzone(ichzone); %Identifier le canal
%                 if strcmp(DATA{isubject}.System,'ISS')
%                     strDet = SDDet2strboxy_ISS(ML(ich,2)); %extraire le nom du détecteur
%                     strSrs = SDPairs2strboxy_ISS(ML(ich,1)); %extraire le nom de la source
%                     idch = strmatch([strDet, ' ',strSrs ],ZoneList,'exact'); %trouver la paire d'optodes dans la liste
%                 elseif strcmp(DATA{isubject}.System,'NIRx') %modKR
%                     strDet = SDDet2strboxy(ML(ich,2)); %extraire le nom du détecteur
%                     strSrs = SDPairs2strboxy(ML(ich,1)); %extraire le nom de la source
%                     idch = strmatch([strDet, ' ',strSrs ],ZoneList,'exact'); %trouver la paire d'optodes dans la liste
%                 end
%                 idliststr =[idliststr,{[strDet, ' ',strSrs ]}]; %liste des paires d'optodes de la zone
%                 idlisti = [idlisti, idch]; %liste des canaux de la zone
%             end
% 
%             for jzone = 1:numel(DATA{isubject}.zone.label) %pour chaque zone, refaire la même chose
%                 idlistj = [];
%                 chzone = DATA{isubject}.zone.plotLst{jzone};
% 
%                 for ichzone = 1:numel(chzone); %pour chaque canal dans la zone
%                     ich = chzone(ichzone);
%                     if strcmp(DATA{isubject}.System,'ISS')
%                         strDet = SDDet2strboxy_ISS(ML(ich,2));
%                         strSrs = SDPairs2strboxy_ISS(ML(ich,1));
%                         idch = strmatch([strDet, ' ',strSrs ],ZoneList,'exact');
%                     elseif strcmp(DATA{isubject}.System,'NIRx') %modKR
%                         strDet = SDDet2strboxy(ML(ich,2)); %extraire le nom du détecteur
%                         strSrs = SDPairs2strboxy(ML(ich,1)); %extraire le nom de la source
%                         idch = strmatch([strDet, ' ',strSrs ],ZoneList,'exact'); %trouver la paire d'optodes dans la liste
%                     end
%                     idlistj = [idlistj, idch]; %liste des canaux de la zone
%                 end
%                 matROI = DATA{isubject}.MAT(idlisti,idlistj); %extraire la connectivité entre les canaux des paires de zones sélectionnés(matrice de ROI)
%                 id = find(matROI==0); %trouver si des paires de canaux ont une connectivité de 0
%                 if isempty(id) %s'il n'y a rien dans la matrice, la mettre NAN
%                     matROI(id)=nan;
%                 end
%                 MATallroi(izone,jzone,isubject) = nanmean(matROI(:)); %faire la moyenne de la connectivité des canaux entre la paire de ROI
%                 if izone==jzone
%                     matnbnanbyizone(izone,isubject)=numel(find(sum(double(isnan(matROI))) ==size(matROI,1)));
%                     matnbtotchbyizone(izone,isubject) = size(matROI,1);
%                 end
%             end 
%         end
%         groupid(isubject)= DATA{isubject}.GR;
%         labelnode = 'z';
%     end
% 
% %     %Création de la ZoneList
% %     zoneuse=DATA{isubject}.zone;
% %     ZoneList = [];
% %     plottmp=[];
% %     plotLst = [];
% %     for izoneList = 1:size(MATall,2)
% %         MLfake(izoneList,1) = izoneList;%source
% %         MLfake(izoneList,2) = 1; %detecteur
% %         MLfake(izoneList,3) = 1;
% %         MLfake(izoneList,4) = 1;
% %         strDet = SDDet2strboxy_ISS(MLfake(izoneList,2));
% %         strSrs = SDPairs2strboxy_ISS(MLfake(izoneList,1));
% %         ZoneLabel{izoneList,1}=zoneuse.label{izoneList};
% %         ZoneList{izoneList,1} = [strDet,' ', strSrs];
% %         plottmp{izoneList} = [izoneList,1];
% %         plotLst{izoneList} = [izoneList];
% %     end
% %     %save zone list associate
% %     zone.plot = plottmp;
% %     zone.plotLst = plotLst;
% %     zone.label = ZoneLabel;
% %     zone.color = zone.color;
% %     zone.ml = MLfake;
% %     zone.chMAT = plotLst;
% % 
% %     save(fullfile(info{isubject,1},['avg', info{isubject,3}]),'zone','-mat');
% %     ZONEid = ['avg', info{isubject,3}];
% 
% clear matROI isubject ML idlisti idlistj izone jzone idliststr  ichzone chzone ich idch strDet strSrs
% end

clear id iname i ich izone ichzone idch idzone iname i ilabel isubject strDet strSrs
save([savepath 'workspacemat.mat'])

clear ext raw txt names MAT matcorr idsubject groupid info xlslistfile filepath name ZONEid labelnode
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Organizing the data')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%Extraire les données démographiques%%%%%
T = readtable(exceltable);
T.Properties.VariableNames{'GROUP_C_0_EXP_1_'} = 'GROUP'; %Renommer des variables
T.Properties.VariableNames{'SEX_M_0_F_1_'} = 'SEX';

part = [];
for p = 1:numel(list_subject)
    ipart = list_subject{p};
    idx = strfind(ipart,'_');
    part{p,1} = ipart(1:idx(1)-1);
end
part = string(part);
%     A=[];
%     A = erase(list_subject,'_HBO_Pearson'); %Ajouter Partcorr si s'applique
%     strfind(list_subject,'_')
%     part = string(A');
gr = groupeall;
for p = 1:numel(part)
    tf = strcmp(T.ID,part(p));
    idx = find(tf);
    age(p,1) = T.AGE(idx);
    sex(p,1) = T.SEX(idx);
    ses(p,1) = T.SES(idx);
end

if importROI == 0 & channelmode == 1
    %ListCh = 1:size(MATall,1);
    lgndch = {'Ch';'c'};
end

if ROImode
    if importROI == 1
        %vérification du fichier zone
        if size(zone.label,1)~=size(zone.plotLst,1)
            zone.label = zone.label';
        else
        end

        %extraire la zone des données de connectivité importées
        roi = [zone.label; zone.plotLst]; 

        %%Créer une liste des labels de zones 
        for r = 1:numel(roi(1,:))
            Listroi{r} = ['R' num2str(r)];
        end
        roi = [roi; Listroi];
        Listroi = string(Listroi);
    end

    if calculateROI == 1 & importROI == 0
        %vérification du fichier zone
        if size(zone.label,1)~=size(zone.plotLst,1)
            zone.label = zone.label';
        else
        end

        lgndroi = {'roiM','roiR','roiF','roiA'; 'm', 'r', 'f', 'a'}; %légende des différentes zones
        ListRoi = zone.label;
        for i = 1:numel(ListRoi)
            ilabel = ListRoi{1,i};
            ListRoi{1,i} = strrep(ilabel(1:end),'_',' ');
        end
        roiM = [ListRoi; zone.plotLst];
        roiR = {'Cortex prefrontal ant G','Cortex prefrontal dorsolat G','Pars triangularis G',...
            'Pars opercularis G','Cortex premoteur G','Cortex moteur prim G','Cortex sensorimoteur prim G',...
            'Gyrus supramarginal G','Cortex gustatif G','Gyrus temporal sup G','Gyrus fusiforme G',...
            'Gyrus temporal med G','Aire temporopolaire G','Cortex prefrontal ant D','Cortex prefrontal dorsolat D',...
            'Pars triangularis D','Pars opercularis D','Cortex premoteur D','Cortex moteur prim D',...
            'Cortex sensorimoteur prim D','Gyrus supramarginal D','Cortex gustatif D','Gyrus temporal sup D',...
            'Gyrus fusiforme D','Gyrus temporal med D','Aire temporopolaire D'; 
            23, [2 3 4 22], 5, 7, [1 8 11 12], [13 19], 21, [17 20], 10, 15, 18, [14 16], [6 9],...
            46, [25 26 27 45], 28, 30, [24 41 42 43], [34 44], 36, [35, 40], 33, 38, 41, [37 39], [29, 32]}; %création manuelle d'une zone
        roiF = {'executif G','moteur G','somato G','auditif G','memoire G','reguemo G',...
            'executif D','moteur D','somato D','auditif D','memoire D','reguemo D';
            [2 3 4 5 7 22 23], [1 8 11 12 13 19], [10 17 20 21], 15, [14 16 18], [6 9],...
            [25 26 27 28 30 45 46], [24 31 34 42 43 44], [33 35 36 40], 38, [37 39 41], [29 32]};
        roiA = {'prefrontal G','frontal G','parietal G','temporal G','prefrontal D','frontal D','parietal D','temporal D';
            [2 3 4 5 7 22 23], [1 8 11 12 13 19], [10 17 20 21], [6 9 14 15 16 18],...
            [25 26 27 28 30 45 46], [24 31 34 42 43 44], [33 35 36 40], [29 32 37 38 39 41]};

        for R = 1:numel(lgndroi(1,:)) %%Créer une liste des labels de zones pour chaque moyennage différent
            if R == 1
                for r = 1:numel(roiM(1,:))
                    ListMroi{r} = [lgndroi{2,R} num2str(r)];
                end
            roiM = [roiM; ListMroi];
            end
             if R == 2
                for r = 1:numel(roiR(1,:))
                    ListRroi{r} = [lgndroi{2,R} num2str(r)];
                end
            roiR = [roiR; ListRroi];
             end
             if R == 3
                for r = 1:numel(roiF(1,:))
                    ListFroi{r} = [lgndroi{2,R} num2str(r)];
                end
            roiF = [roiF; ListFroi];
             end
             if R == 4
                for r = 1:numel(roiA(1,:))
                    ListAroi{r} = [lgndroi{2,R} num2str(r)];
                end
            roiA = [roiA; ListAroi];
             end
        end
        roi = {roiM, roiR, roiF, roiA}; %tout mettre dans la même variable

        clear roiM roiR roiF roiA R r
    end
end

clear list_subject groupeall T tf A connectivity idx ipart p ilabel i
    
%PAR CHANNELS%%%%%%%%%%%%%%%%%%%%%%%%%%%
%création des variables de connectivité pour chaque paire de canaux%

if importROI == 0 & channelmode == 1
    %création des labels de chaque paire de canaux
    x = 1; 
    for c = 1:size(MATall,1)
        for cc = (c + 1):size(MATall,1)
            %labelch{x} = [num2str(c) '-' num2str(cc)];
            labelch{x} = [num2str(ListCh(c)) '-' num2str(ListCh(cc))];
            x = x + 1; 
        end
    end
    labelch = string(labelch);


    %extraire les données des matrices initiales dans le nouveau format
    x = 1;
    for p = 1:size(MATall,3) 
        for c = 1:size(MATall,1)
            for cc = (c + 1):size(MATall,1)
                datach(p,x) = MATall(c,cc,p);
                x = x + 1;
            end
        end
        x = 1;
    end
    NC = size(MATall,1);
    NCP = size(datach,2);
    tblch = [table(part, gr, age, sex, ses) array2table([datach],'VariableNames', labelch)]; %créer un tableau avec les données

    clear idx tf p c cc x;
end

%PAR ROI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ROImode
    if importROI == 1
        %création des labels de chaque paire de ROI
        x = 1; 
        for r = 1:size(MATall,1)
            for rr = (r + 1):size(MATall,1)
                labelroi{x} = ['R' num2str(r) '-' 'R' num2str(rr)];
                x = x + 1; 
            end
        end
        labelroi = string(labelroi);

        %extraire les données des matrices initiales dans le nouveau format
        x = 1;
        for p = 1:size(MATall,3) 
            for r = 1:size(MATall,1)
                for rr = (r + 1):size(MATall,1)
                    dataroiALL(p,idx) = MATall(r,rr,p);
                    x = x + 1;
                end
            end
            x = 1;
        end

        tblroi = [table(part, gr, age, sex, ses) array2table([dataroiALL],'VariableNames', labelroi)]; %créer un tableau avec les données

        clear idx tf p r rr x;
    end
    
    if calculateROI == 1 & importROI == 0
        %création des labels pour les paires de ROI%%%
        x = 1;
        for R = 1: numel(roi)
            for r = 1:numel(roi{:,R}(1,:))
                for rr = (r + 1):numel(roi{:,R}(1,:))
                    labelroiALL{1,R}{1,x} = [lgndroi{2,R} num2str(r) '-' lgndroi{2,R} num2str(rr)];
                    x = x + 1; 
                end
            end
            x = 1;
        end

        %création de la liste des paires de channels pour chaque paire de ROI%%
        x = 1;
        y = 1;
        for R = 1:numel(roi) %for each ROI montage
            for m = 1:numel(roi{:,R}(1,:)) %for each R1
                for n = (m + 1):numel(roi{:,R}(1,:)) %for each R2
                    for r = 1:numel(roi{:,R}{2,m}) %for each channel in R1
                        for rr = 1:numel(roi{:,R}{2,n}) %for each channel in R2
                            chroi{1,R}{1,x}(1,y) = join([num2str(roi{:,R}{2,m}(1,r)) "-" num2str(roi{:,R}{2,n}(1,rr))],"");
                            y = y + 1;
                        end
                    end
                x = x + 1;
                y = 1;
                end
            end
            x = 1;
        end

        % Extraire les valeurs des données originales et les organiser par ROI
        z = 1;
        for p = 1:size(MATall,3) %for each participant
            for R = 1:numel(roi) %for each ROI montage
                for r = 1:size(MATall,1) %for each channel1
                    for rr = (r + 1):size(MATall,1) %for each channel2
                        x = [num2str(r) '-' num2str(rr)]; %extract channel pair name
                        y = [num2str(rr) '-' num2str(r)];
                        for n = 1:numel(chroi{1,R}) %for each roi pair
                            tf = strcmp(x, chroi{1,R}{1,n}); %find the chpair in the roich list
                            if any(tf)
                                idx = find(tf);
                                %   idx = strfind(labelch, x);
                                roiALL{1,R}{p,n}(1,idx) = MATallold(r,rr,p); %extract its connectivity
                            end
                            tf = strcmp(y, chroi{1,R}{1,n});
                            if any(tf)
                                idx = find(tf);
                                %   idx = strfind(labelch, x);
                                roiALL{1,R}{p,n}(1,idx) = MATallold(r,rr,p);
                            end
                        end
                    end
                end
            end
        end

        %Extraire chaque colonne du cell array en matrice individuelle et faire la
        %moyenne de la connectivité des canaux de chaque roi pour chaque participant

        for R = 1:numel(roi) %for each ROI montage
            for n = 1:numel(chroi{1,R}) %for each roi pair
                A = cell2mat(roiALL{1,R}(:,n)); %extract connectivity values for each roipair
                roimean = nanmean(A,2); %compute the mean
                dataroiALL{1,R}(:,n) = roimean;
            end
        end

        %mettre les données moyennées dans un tableau
        tbldataMroi = [table(part, gr, age, sex, ses) array2table([dataroiALL{1,1}], 'VariableNames', string(labelroiALL{1,1}))]; 
        tbldataRroi = [table(part, gr, age, sex, ses) array2table([dataroiALL{1,2}], 'VariableNames', string(labelroiALL{1,2}))];
        tbldataFroi = [table(part, gr, age, sex, ses) array2table([dataroiALL{1,3}], 'VariableNames', string(labelroiALL{1,3}))];
        tbldataAroi = [table(part, gr, age, sex, ses) array2table([dataroiALL{1,4}], 'VariableNames', string(labelroiALL{1,4}))];

        clear idx tf p n m R r rr x y z A roimean roiALL;
    end
end

%SAUVEGARDES
if importROI == 0 & channelmode == 1 & calculateROI == 0
    save([savepath date '_datamat.mat'],'datach');
    save([savepath date '_tables.mat'],'tblch');
    
elseif importROI == 0 & calculateROI == 1 & channelmode == 0
    save([savepath date '_datamat.mat'],'dataroiALL');
    save([savepath date '_tables.mat'],'tbldataMroi','tbldataRroi','tbldataFroi','tbldataAroi');

elseif importROI == 0 & calculateROI == 1 & channelmode == 1
    save([savepath date '_datamat.mat'],'datach','dataroiALL');
    save([savepath date '_tables.mat'],'tblch','tbldataMroi','tbldataRroi','tbldataFroi','tbldataAroi');
    
elseif importROI == 1
    save([savepath date '_datamat.mat'],'dataroiALL');
    save([savepath date '_tables.mat'],'tblroi');
end

%% Moyennes pour graph
%données de FC séparée par groupe
datachG1 = datach(idG1,:);
datachG2 = datach(idG2,:);
datameanchG1 = mean(datachG1,'omitnan');
datameanchG2 = mean(datachG2,'omitnan');
datachG1G2 = datameanchG1 - datameanchG2;
%datachG2G1 = datameanchG2 - datameanchG1;

if calculateROI == 1
    for R = 1:numel(roi)
        %données de FC séparée par groupe
        dataroiG1{:,R} = dataroiALL{1,R}(idG1,:);
        dataroiG2{:,R} = dataroiALL{1,R}(idG2,:);
        datameanroiG1{:,R} = nanmean(dataroiG1{:,R});
        datameanroiG2{:,R} = nanmean(dataroiG2{:,R});
        dataroiG1G2{:,R} = datameanroiG1{:,R} - datameanroiG2{:,R};
        %dataroiG2G1{:,R} = meanroiG2{:,R} - meanroiG1{:,R};
        
        %convertir les données en matrices
        for p = 1:size(MATall,3)
            MATroiALL{1,R}(:,:,p) = MATvTBL.TBL2MAT(dataroiALL{1,R}(p,:));
        end
        MATroiG1{:,R} = MATroiALL{1,R}(:,:,idG1);
        MATroiG2{:,R} = MATroiALL{1,R}(:,:,idG2);
        MATroimeanG1{:,R} = mean(MATroiG1{:,R}, 3, 'omitnan');
        MATroimeanG2{:,R} = mean(MATroiG2{:,R}, 3, 'omitnan');
        MATroiG1G2{:,R} = MATroimeanG1{:,R} - MATroimeanG2{:,R};
    end
    clear R p
end

if importROI == 1
    %données de FC séparée par groupe
    dataroiG1 = dataroiALL(idG1,:);
    dataroiG2 = dataroiALL(idG2,:);
    datameanroiG1 = nanmean(dataroiG1);
    datameanroiG2 = nanmean(dataroiG2);
    dataroiG1G2 = datameanroiG1 - datameanroiG2;
    dataroiG2G1 = datameanroiG2 - datameanroiG1;
end

save([savepath 'workspace.mat'])

X = ['Results saved in ', savepath];
disp(X)

clear all

toc