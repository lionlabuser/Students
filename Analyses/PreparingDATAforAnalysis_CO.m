%%%%%%%%%%%%%%%%%PREPARING CORR DATA FOR ANALYSIS%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
savepath='C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Stats\CORR0,01_0,08\Test\';
if ~isfolder(savepath)
    mkdir(savepath)
end

zonemode = 0;
connectivity = 'CORR'; %Modify 'COH' OR 'CORR'
xlslistfile = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\CORRmatrice0,01_0,08\Channels\Subjectlist N=54.xlsx'; %Fichier excel avec le dossier des matrices, leur nom et le groupe
exceltable = 'Participants list.xlsx'; %%%% Fichier excel avec les données démographiques d'intérêt

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
ch = 1:46;

save([savepath date '_MATall.mat'],'MATall');
save([savepath date '_DATA.mat'],'MATall');

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

save([savepath 'workspacemat.mat'])

clear ext raw txt iname names MAT matcorr id isubject idsubject groupid info xlslistfile filepath name ZONEid ZoneList labelnode MATstdG1 MATstdG1 MATvarG1 MATvarG2
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

if zonemode
    lgndroi = {'roiM','roiR','roiF','roiA'; 'm', 'r', 'f', 'a'}; %légende des différentes zones
    roiM = [zone.label; zone.plotLst]; %% extraire la zone des données de connectivité importées
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

clear list_subject groupeall zone T tf
    
%PAR CHANNELS%%%%%%%%%%%%%%%%%%%%%%%%%%%
%création des variables de connectivité pour chaque paire de canaux%

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
            chALL(p,idx) = MATall(c,cc,p);
        end
    end
end

tblch = [table(part, gr, age, sex, ses) array2table([chALL],'VariableNames', labelch)]; %créer un tableau avec les données

clear idx tf p c cc x;

%PAR ROI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if zonemode
    %création des labels pour les paires de ROI%%%
    x = 1;
    for R = 1: numel(roi)
        for r = 1:numel(roi{:,R}(1,:))
            for rr = (r + 1):numel(roi{:,R}(1,:))
                %labelroi{1,x} = {char(lgndroi(2,R)) num2str(r) '-' char(lgndroi(2,R)) num2str(rr)};
                labelroi{1,R}{1,x} = [lgndroi{2,R} num2str(r) '-' lgndroi{2,R} num2str(rr)];
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
            roiALLmean{1,R}(:,n) = roimean;
        end
    end

    %mettre les données moyennées dans un tableau
    tblmeanMroi = [table(part, gr, age, sex, ses) array2table([roiALLmean{1,1}], 'VariableNames', string(labelroi{1,1}))]; 
    tblmeanRroi = [table(part, gr, age, sex, ses) array2table([roiALLmean{1,2}], 'VariableNames', string(labelroi{1,2}))];
    tblmeanFroi = [table(part, gr, age, sex, ses) array2table([roiALLmean{1,3}], 'VariableNames', string(labelroi{1,3}))];
    tblmeanAroi = [table(part, gr, age, sex, ses) array2table([roiALLmean{1,4}], 'VariableNames', string(labelroi{1,4}))];

    clear idx tf p n m R r rr x y z A roimean roiALL;

    save([savepath date '_datamat.mat'],'chALL','roiALLmean');
    save([savepath date '_tables.mat'],'tblch','tblmeanMroi','tblmeanRroi','tblmeanFroi','tblmeanAroi');
end
    
%%
%%Stats descriptives%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Descriptive Statistics')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Channels%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%normalité et valeurs extrêmes%
meanG1ch = nanmean(chALL(idG1,:));
meanG2ch = nanmean(chALL(idG2,:));
stdG1ch = nanstd(chALL(idG1,:));
stdG2ch = nanstd(chALL(idG2,:));
sknG1ch = skewness(chALL(idG1,:));
sknG2ch = skewness(chALL(idG2,:));
krtG1ch = kurtosis(chALL(idG1,:));
krtG2ch = kurtosis(chALL(idG2,:));

tblgrmeanch = [array2table([meanG1ch],'VariableNames',labelch); array2table([meanG2ch],'VariableNames',labelch)];
tblgrmeanch.Properties.RowNames = {'G1','G2'};
tbldescrch = array2table([meanG1ch; meanG2ch; stdG1ch; stdG2ch; sknG1ch; sknG2ch; krtG1ch; krtG2ch],'VariableNames',labelch,'RowNames',{'meanG1','meanG2','stdG1','stdG2','sknG1','sknG2','krtG1','krtG2'});

tf = tbldescrch{{'sknG1','sknG2','krtG1','krtG2'},:} >2 | tbldescrch{{'sknG1','sknG2','krtG1','krtG2'},:} <-2;
n_abnormal = sum(any(tf));
p_abnormal = n_abnormal/numel(tbldescrch(1,:))*100;
fprintf('%.1f percent of the channels variables are not respecting the normal distribution\n',p_abnormal);

nanzscore = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));
zchALL = nanzscore(chALL);
tf = zchALL > 3.29 | zchALL <-3.29;
n_extreme = sum(any(tf));
p_extreme = n_extreme/numel(zchALL)*100;
fprintf('%.1f percent of the channels variables have extreme values\n',p_extreme);

%%NAN Values%%%%%%%%%%
npart_nan = sum(isnan(chALL),2); %nb de NAN channels par participant
percentnan = sum(npart_nan)/numel(chALL)*100; %pourcentage NAN channels pour tous participants
fprintf('%.1f percent of the data is a missing value\n',percentnan)

[sortpartnan,indexpartnan] = sort(npart_nan,'descend'); %ordonner les participants + au - NAN
mostnanpart = sortpartnan > 50/100*size(chALL,2); %50% ou plus de données manquantes
pbnanpart = part(indexpartnan(mostnanpart)); %identifier les participants
fprintf('%s has more than 50 percent missing values\n',pbnanpart)

clear npart_nan percentnan sortpartnan indexpartnan mostnanpart pbnanpart

R = 4; %Aroi
x = 1;
for p = 1:size(MATall,3) %pour chaque participant
    tfnanch = mean(isnan(MATall(:,:,p)),1) == 1; %trouver les canaux complètement NAN
    idxnanch = find(tfnanch); %trouver leur position dans la matrice
    nanch{x,1} = idxnanch; %liste des canaux manquants pour chaque participant
    x = x + 1;
end

nb = 0;
for c = 1:size(MATall,1) %pour chaque canal
    for n = 1 : size(nanch,1) %pour chaque participant
        tf = sum(nanch{n,1} == c); %déterminer si le canal est manquant
        nb = nb + tf; % nombre de fois ou le canal est manquant
    end
   nanfreq{c,1} = nb; %liste du nb de fois ou chaque canal est manquant 
   nb = 0;
    if zonemode
       for n = 1:size(roi{1,R},2) %pour chaque roi
           tf = roi{1,R}{2,n} == c; %trouver si le canal fait partie de la roi
           if any(tf) %si oui
               nanfreq{c,2} = roi{1,R}(1,n); %ajouter le label
           end
       end
   end
end

if zonemode
    for n = 1:size(roi{1,R},2) % pour chaque roi
        tf = strcmp(roi{1,R}(1,n),[nanfreq{:,2}]); %trouver chaque canal qui provient de chaque région
        idx = find(tf); %trouver la position de chaque canal
        nantot(n,1) = sum([nanfreq{idx,1}]); %faire la somme des canaux NAN de chaque région
    end
    nantot = num2cell(nantot);
    nantot(:,2) = roi{1,R}(1,:);
end

[sortchnan,indexchnan] = sort(cell2mat(nanfreq(:,1)),'descend');
mostnanch = sortchnan > 30/100*size(MATall,3); %30% ou plus de données manquantes
pbnanch = ch(indexchnan(mostnanch));
fprintf('Channel %.0f has more than 30 percent missing values\n',pbnanch)

figure
X = categorical(1:size(MATall,1));
bar(X,cell2mat(nanfreq(:,1)))
yline(mean(cell2mat(nanfreq(:,1)),1))
ylabel('Total number of missing')
title('Total number of missing channels')

if zonemode
    figure
    X = categorical(roi{1,R}(1,:));
    bar(X,cell2mat(nantot(:,1)))
    yline(mean(cell2mat(nantot(:,1)),1))
    ylabel('Total number of channels')
    title('Total number of channels with missing data in each region')
    
    clear nantot
end
    
clear R x p c n idx tf nanch tfnanch idxnanch nb sortchnan indexchnan mostnanch pbnanch nanfreq

%%%%%Roi%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if zonemode
    %normalité et valeurs extrêmes%
    ALLroi = [];
    for R = 1:4
        ALLroi = [ALLroi, roiALLmean{1,R}];
    end

        meanG1roi = nanmean(ALLroi(idG1,:));
        meanG2roi = nanmean(ALLroi(idG2,:));
        stdG1roi = nanstd(ALLroi(idG1,:));
        stdG2roi = nanstd(ALLroi(idG2,:));
        sknG1roi = skewness(ALLroi(idG1,:));
        sknG2roi = skewness(ALLroi(idG2,:));
        krtG1roi = kurtosis(ALLroi(idG1,:));
        krtG2roi = kurtosis(ALLroi(idG2,:));
    tblgrmeanroi = [array2table([meanG1roi],'VariableNames',[labelroi{1,1:4}]); array2table([meanG2roi],'VariableNames',[labelroi{1,1:4}])];
    tblgrmeanroi.Properties.RowNames = {'G1','G2'};
    tbldescrroi = array2table([meanG1roi; meanG2roi; stdG1roi; stdG2roi; sknG1roi; sknG2roi; krtG1roi; krtG2roi],'VariableNames',[labelroi{1,1:4}],'RowNames',{'meanG1','meanG2','stdG1','stdG2','sknG1','sknG2','krtG1','krtG2'});

    tf = tbldescrroi{{'sknG1','sknG2','krtG1','krtG2'},:} >2 | tbldescrroi{{'sknG1','sknG2','krtG1','krtG2'},:} <-2;
    n_abnormal = sum(any(tf));
    p_abnormal = n_abnormal/numel(tbldescrroi(1,:))*100;
    fprintf('%.1f percent of the ROIs variables are not respecting the normal distribution\n',p_abnormal);

    zALLroi = nanzscore(ALLroi);
    tf = zALLroi > 3.29 | zALLroi <-3.29;
    n_extreme = sum(any(tf));
    p_extreme = n_extreme/numel(zALLroi)*100;
    fprintf('%.1f percent of the ROIs variables have extreme values\n',p_extreme);

    meanG1Mroi = nanmean(roiALLmean{1,1}(idG1,:));
    meanG2Mroi = nanmean(roiALLmean{1,1}(idG2,:));
    meanG1Rroi = nanmean(roiALLmean{1,2}(idG1,:));
    meanG2Rroi = nanmean(roiALLmean{1,2}(idG2,:));
    meanG1Froi = nanmean(roiALLmean{1,3}(idG1,:));
    meanG2Froi = nanmean(roiALLmean{1,3}(idG2,:));
    meanG1Aroi = nanmean(roiALLmean{1,4}(idG1,:));
    meanG2Aroi = nanmean(roiALLmean{1,4}(idG2,:));

    clear X tf R ALLroi n_abnormal n_extreme p_abnormal p_extreme stdG1ch stdG2ch sknG1ch sknG2ch krtG1ch krtG2ch meanG1roi meanG2roi stdG1roi stdG2roi sknG1roi sknG2roi krtG1roi krtG2roi zchALL zroiALLmean zALLroi
end
    
save([savepath 'workspace.mat'])

X = ['Results saved in ', savepath];
disp(X)

%clear all

toc