%%%%%%%%%%%%%%%%%PREPARING CORR DATA FOR ANALYSIS%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
savepath='C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Stats\CORR0,01_0,08\PhysioSatRespEKG\';
if ~isfolder(savepath)
    mkdir(savepath)
end

importROI = 0; %0=import channels; 1=importROI%
calculateROI = 1; %Calculer la connectivité moyennée par ROI (seulement possible si channels importés)
channelmode = 1; %Faire les analyses sur les canaux
ROImode = 1; %Faire les analyses sur les ROI

connectivity = 'CORR'; %Modify 'COH' OR 'CORR'
xlslistfile = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\CORRmatrice0,01_0,08\Channels\PhysioRegressed_SatRespEKG\Subjectlist N=54.xlsx'; %Fichier excel avec le dossier des matrices, leur nom et le groupe
exceltable = 'C:\data\Malnutrition\Resting\NIRS\Participants list.xlsx'; %%%% Fichier excel avec les données démographiques d'intérêt

%%%%%% from StatMatrices of LIONIRS toolbox%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Loading the data')
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
for isubject=2:size(info,1)
    id = isubject-1;
    MAT = load(fullfile(info{isubject,1},[info{isubject,2},'.mat']));
    DATA{id}.ZoneList = MAT.ZoneList;
    if isfield(MAT, 'meancorr')
        matcorr = MAT.meancorr;
        DATA{id}.MAT = matcorr;
    end
    DATA{id}.name = info{isubject,2};
    DATA{id}.GR = info{isubject,4};

    %Load the zones
    load(fullfile(info{isubject,1}, info{isubject,3}),'-mat');
    names = fieldnames(zone);
    for iname = 1:numel(names)
        eval(['DATA{id}.zone.',names{iname},' =zone.',names{iname},';']);
    end
    
    %création liste des noms de sujets et liste du groupe de chacun%
    list_subject{id} =DATA{id}.name;
    groupeall = [groupeall; info{isubject,4}];
%     clear MAT
end

%création d'une matrice avec les matrices de tous les sujets et d'une liste du groupe de chacun%
idsubject = 1:numel(groupeall); %channel mode
    for isubject = 1:numel(groupeall)
        MATall(:,:,isubject)=DATA{isubject}.MAT;
        groupid(isubject)= DATA{idsubject(isubject)}.GR;
    end
%sélectionner le nom du fichier zone et la liste des canaux du dernier participant%     
ZONEid = [info{end,3}];
ZoneList =  DATA{end}.ZoneList;
labelnode = 'c';


% %%%%%%%% NEWLY ADDED, Calculate ROI from StatMatrix%%%%%%%%%%%%%%%%
% %création d'une matrice pour les zones de tous les sujets
% MATall = zeros(numel(DATA),numel(DATA{id}.zone.label),numel(DATA{id}.zone.label));
% for isubject = 1:numel(groupeall) %Pour chaque sujet
%     List = DATA{isubject}.ZoneList; %Extraire la ZoneList
%     
%     for izone = 1:numel(DATA{isubject}.zone.label) %Pour chaque zone
%         ML = DATA{isubject}.zone.ml; %Extraire le ML (paire d'optodes)
%         DATA{isubject}.zone.plotLst; %Extraire le plotLst (id canal)
%         idlisti = [];
%         idliststr = [];
%         chzone = DATA{isubject}.zone.plotLst{izone}; %Extraire les canaux de la zone
%         
%         for ichzone = 1:numel(chzone); %pour chaque canal dans la zone
%             ich = chzone(ichzone); %Identifier le canal
%             if strcmp(DATA{isubject}.System,'ISS')
%                 strDet = SDDet2strboxy_ISS(ML(ich,2)); %extraire le nom du détecteur
%                 strSrs = SDPairs2strboxy_ISS(ML(ich,1)); %extraire le nom de la source
%                 idch = strmatch([strDet, ' ',strSrs ],List,'exact'); %trouver la paire d,optodes dans la liste
%             end
%             idliststr =[idliststr,{[strDet, ' ',strSrs ]}]; %liste des paires d'optodes de la zone
%             idlisti = [idlisti, idch]; %liste des canaux de la zone
%         end
%         
%         for jzone = 1:numel(DATA{isubject}.zone.label) %pour chaque zone, refaire la même chose
%             idlistj = [];
%             chzone = DATA{isubject}.zone.plotLst{jzone};
%             
%             for ichzone = 1:numel(chzone); %pour chaque canal dans la zone
%                 ich = chzone(ichzone);
%                 if strcmp(DATA{isubject}.System,'ISS')
%                     strDet = SDDet2strboxy_ISS(ML(ich,2));
%                     strSrs = SDPairs2strboxy_ISS(ML(ich,1));
%                     idch = strmatch([strDet, ' ',strSrs ],List,'exact');
%                 end
%                 idlistj = [idlistj, idch]; %liste des canaux de la zone
%             end
%             matROI = DATA{isubject}.MAT(idlisti,idlistj); %extraire la connectivité entre les canaux des paires de zones sélectionnés(matrice de ROI)
%             id = find(matROI==0); %trouver si des paires de canaux ont une connectivité de 0
%             if isempty(id) %s'il n'y a rien dans la matrice, la mettre NAN
%                 matROI(id)=nan;
%             end
%             MATall(isubject,izone,jzone) = nanmean(matROI(:)); %faire la moyenne de la connectivité des canaux entre la paire de ROI
%             if izone==jzone
%                 matnbnanbyizone(isubject,izone)=numel(find(sum(double(isnan(matROI))) ==size(matROI,1)));
%                 matnbtotchbyizone(isubject,izone) = size(matROI,1);
%             end
%         end
%     end
%     groupid(isubject)= DATA{isubject}.GR;
%     labelnode = 'z';
% end

save([savepath date '_MATall.mat'],'MATall');
%save([savepath date '_DATA.mat'],'MATall');

idG1 = find(groupeall==1);
MATallG1 = MATall(:,:,idG1);
MATmeanG1 = squeeze(nanmean(MATallG1,3));
MATstdG1 = std(MATall(:,:,idG1),0,3,'omitnan');
MATvarG1 = var(MATall(:,:,idG1),0,3,'omitnan');
%MATpctstdG1 = (MATstdG1./MATmeanG1)*100;
%pbstd = sum(MATpctstdG1 > 75);
%pbstdchG1 = ch(pbstd > 10);

idG2 = find(groupeall==2);
MATallG2 = MATall(:,:,idG2);
MATmeanG2 = squeeze(nanmean(MATallG2,3));
MATstdG2 = std(MATall(:,:,idG2),0,3,'omitnan');
MATvarG2 = var(MATall(:,:,idG2),0,3,'omitnan');
%MATpctstdG2 = (MATstdG2./MATmeanG2)*100;
%pbstd = sum(MATpctstdG2 > 75);
%pbstdchG2 = ch(pbstd > 10);

clear id iname

save([savepath 'workspacemat.mat'])

clear ext raw txt names MAT matcorr isubject idsubject groupid info xlslistfile filepath name ZONEid ZoneList labelnode MATstdG1 MATstdG1 MATvarG1 MATvarG2
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Organizing the data')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%Extraire les données démographiques%%%%%
T = readtable(exceltable);
T.Properties.VariableNames{'GROUP_C_0_EXP_1_'} = 'GROUP'; %Renommer des variables
T.Properties.VariableNames{'SEX_M_0_F_1_'} = 'SEX';

if isequal(connectivity,'CORR') % créer des variables différente pour chaque variable démographique
    A=[];
    A = erase(list_subject,'_HBO_Pearson');
    part = string(A');
    gr = groupeall;
    for p=1:numel(part)
        tf = strcmp(T.ID,part(p));
        idx = find(tf);
        age(p,1) = T.AGE(idx);
        sex(p,1) = T.SEX(idx);
        ses(p,1) = T.SES(idx);
    end
end
if isequal(connectivity,'COH') % créer des variables différente pour chaque variable démographique
    A=[];
    A = erase(list_subject,'_HBO_COH FFT');
    part = string(A');
    gr = groupeall;
    for p=1:numel(part)
        tf = strcmp(T.ID,part(p));
        idx = find(tf);
        age(p,1) = T.AGE(idx);
        sex(p,1) = T.SEX(idx);
        ses(p,1) = T.SES(idx);
    end
end

if importROI == 0 & channelmode == 1
    Listch = 1:size(MATall,1);
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
        roiM = [zone.label; zone.plotLst];
        roiR = {'Cortex_prefrontal_ant_G','Cortex_prefrontal_dorsolat_G','Pars_triangularis_G',...
            'Pars_opercularis_G','Cortex_premoteur_G','Cortex_moteur_prim_G','Cortex_sensorimoteur_prim_G',...
            'Gyrus_supramarginal_G','Cortex_gustatif_G','Gyrus_temporal_sup_G','Gyrus_fusiforme_G',...
            'Gyrus_temporal_med_G','Aire_temporopolaire_G','Cortex_prefrontal_ant_D','Cortex_prefrontal_dorsolat_D',...
            'Pars_triangularis_D','Pars_opercularis_D','Cortex_premoteur_D','Cortex_moteur_prim_D',...
            'Cortex_sensorimoteur_prim_D','Gyrus_supramarginal_D','Cortex_gustatif_D','Gyrus_temporal_sup_D',...
            'Gyrus_fusiforme_D','Gyrus_temporal_med_D','Aire_temporopolaire_D'; 
            23, [2 3 4 22], 5, 7, [1 8 11 12], [13 19], 21, [17 20], 10, 15, 18, [14 16], [6 9],...
            46, [25 26 27 45], 28, 30, [24 41 42 43], [34 44], 36, [35, 40], 33, 38, 41, [37 39], [29, 32]}; %création manuelle d'une zone
        roiF = {'executif_G','moteur_G','somato_G','auditif_G','memoire_G','reguemo_G',...
            'executif_D','moteur_D','somato_D','auditif_D','memoire_D','reguemo_D';
            [2 3 4 5 7 22 23], [1 8 11 12 13 19], [10 17 20 21], 15, [14 16 18], [6 9],...
            [25 26 27 28 30 45 46], [24 31 34 42 43 44], [33 35 36 40], 38, [37 39 41], [29 32]};
        roiA = {'prefrontal_G','frontal_G','parietal_G','temporal_G','prefrontal_D','frontal_D','parietal_D','temporal_D';
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

clear list_subject groupeall zone T tf A connectivity
    
%PAR CHANNELS%%%%%%%%%%%%%%%%%%%%%%%%%%%
%création des variables de connectivité pour chaque paire de canaux%

if importROI == 0 & channelmode == 1
    %création des labels de chaque paire de canaux
    x = 1; 
    for c = 1:size(MATall,1)
        for cc = (c + 1):size(MATall,1)
            labelch{x} = [num2str(c) '-' num2str(cc)];
            x = x + 1; 
        end
    end
    labelch = string(labelch);

    %extraire les données des matrices initiales dans le nouveau format
    for p=1:size(MATall,3) 
        for c=1:size(MATall,1)
            for cc =(c + 1):size(MATall,1)
                x = [num2str(c) '-' num2str(cc)];
                tf = strcmp(x, labelch);
                idx = find(tf);
    %            idx = strfind(labelch, x);
                datach(p,idx) = MATall(c,cc,p);
            end
        end
    end

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
        for p=1:size(MATall,3) 
            for r=1:size(MATall,1)
                for rr =(r + 1):size(MATall,1)
                    x = ['R' num2str(r) '-' 'R' num2str(rr)];
                    tf = strcmp(x, labelroi);
                    idx = find(tf);
        %            idx = strfind(labelch, x);
                    dataroi(p,idx) = MATall(r,rr,p);
                end
            end
        end

        tblroi = [table(part, gr, age, sex, ses) array2table([dataroi],'VariableNames', labelroi)]; %créer un tableau avec les données

        clear idx tf p r rr x;
    end
    
    if calculateROI == 1 & importROI == 0
        %création des labels pour les paires de ROI%%%
        x = 1;
        for R = 1: numel(roi)
            for r = 1:numel(roi{:,R}(1,:))
                for rr = (r + 1):numel(roi{:,R}(1,:))
                    %labelroi{1,x} = {char(lgndroi(2,R)) num2str(r) '-' char(lgndroi(2,R)) num2str(rr)};
                    labelroiALL{1,R}{1,x} = [lgndroi{2,R} num2str(r) '-' lgndroi{2,R} num2str(rr)];
                    x = x + 1; 
                end
            end
            x = 1;
        end

        %création de la liste des paires de channels pour chaque paire de ROI%%
        x = 1;
        y = 1;
        for R = 1:numel(roi) %1-4
            for m = 1:numel(roi{:,R}(1,:)) %1-14
                for n = (m + 1):numel(roi{:,R}(1,:)) %2-14
                    for r = 1:numel(roi{:,R}{2,m}) %1-5
                        for rr = 1:numel(roi{:,R}{2,n}) %1-2
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
        for p=1:54
            for R = 1:numel(roi)
                for r=1:46
                    for rr =(r + 1):46
                        x = [num2str(r) '-' num2str(rr)];
                        y = [num2str(rr) '-' num2str(r)];
                        for n = 1:numel(chroi{1,R})
                            tf = strcmp(x, chroi{1,R}{1,n});
                            if any(tf)
                                idx = find(tf);
                                %   idx = strfind(labelch, x);
                                roiALL{1,R}{p,n}(1,idx) = MATall(r,rr,p);
                            end
                            tf = strcmp(y, chroi{1,R}{1,n});
                            if any(tf)
                                idx = find(tf);
                                %   idx = strfind(labelch, x);
                                roiALL{1,R}{p,n}(1,idx) = MATall(r,rr,p);
                            end
                        end
                    end
                end
            end
        end

        %Extraire chaque colonne du cell array en matrice individuelle et faire la
        %moyenne de la connectivité des canaux de chaque roi pour chaque participant

        for R = 1:numel(roi)
            for n = 1:numel(chroi{1,R})
                A = cell2mat(roiALL{1,R}(:,n));
                roimean = nanmean(A,2);
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
    save([savepath date '_datamat.mat'],'dataroi');
    save([savepath date '_tables.mat'],'tblroi');
end

save([savepath 'workspace.mat'])

X = ['Results saved in ', savepath];
disp(X)

clear all

toc