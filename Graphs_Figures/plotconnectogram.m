function plotconnectogram(fileorderconnectogram,MAT,List,label,plotLst,ML)

%%%from GUI_LookMatrice of LIONIRS toolbox%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% id = 1;
% MAT = MATsigch; %DATA{id}.MAT %loader matrice
% List = strvcat(DATA{id}.ZoneList); %liste des paires SD
% ML = DATA{id}.zone.ml; %Loader S/D/ROI/Gr
% plotLst = DATA{id}.zone.plotLst;
% label =  DATA{id}.zone.label;
% fileorderconnectogram  = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Mixte.txt';

fid = fopen(fileorderconnectogram,'r');
id = 1;
while ~feof(fid)
    tline = fgetl(fid);
    listselected{id} = tline; %extract list from connectogramfile
    id = id+1;
end
fclose(fid);

%POUR IDENTIFICATION DE TOUTE LES ZONES
%verification des listes de canaux
idzone = [];
idlist = [];
idlabelall = [];
for izone = 1:numel(plotLst) % pour chaque zone
    chzone =  plotLst{izone}; %extraire les canaux
    for ichzone = 1:numel(chzone) % pour chaque canal
        ich = chzone(ichzone); % prendre un canal
        strDet = SDDet2strboxy(ML(ich,2)); %trouver détecteur
        strSrs = SDPairs2strboxy(ML(ich,1)); % trouver paire
        idch = strmatch([strDet, ' ',strSrs ],List,'exact'); %trouver position dans liste
        idlist = [idlist, idch]; % importer la position de chaque canal
        if ichzone==1 %lorsque premier canal de la zone sélectionné
            idzone =[idzone, izone]; %creer une liste avec le num de zone
        else
            idzone =[idzone, 0]; % si pas premier canal, mettre un 0
        end
    end
    idlabelall = [idlabelall, {[label{izone}]}]; %creer une liste avec le label de toutes les zones%
end

%verification des labels et des canaux de chaque zone%
for adji = 1:numel(idlabelall) % pour chaque zone
    for adjj = 1:numel(idlabelall) % pour chaque zone
        labelzone = idlabelall{adji}; %extraire le label de chaque zone
        x = strmatch({labelzone} ,idlabelall, 'exact'); %voir la position de ce label dans le fichier avec tous les labels%
        labelzone = idlabelall{adjj};
        y = strmatch({labelzone} ,idlabelall, 'exact');
        if isempty(x)|isempty(y)
            msgbox('problem zone in subject')
        end
        chzone = plotLst{x}; %nom des channels pour chaque zone
        idlisti = [];
        for ichzone = 1:numel(chzone); % nombre de canal dans chaque zone
            ich = chzone(ichzone); %prendre chaque canal
            strDet = SDDet2strboxy(ML(ich,2)); %trouver D ISS
            strSrs = SDPairs2strboxy(ML(ich,1)); %trouver paire ISS
            idch = strmatch([strDet, ' ',strSrs ],List,'exact'); %trouver le canal dans la liste
            idlisti = [idlisti, idch]; %creer une liste de la position de chaque canal dans la liste
        end
        if numel(x)>1
            x = x(1);
            msgbox('Attention the zone are duplicated')
        else
            chMAT{x} = idlisti; %fichier contenant la position de chaque canal dans la liste
        end
    end
end

%%ICI POUR CHANGER MANUELLEMENT LA LISTE%%
%     listok = {'Fp1-F7' 'F7-T3' 'T3-T5' 'Fp1-F3' 'Fp2-F8' 'T4-T6' 'Fp2-F4' 'F4-C4' 'C4-P4' 'P4-O2'}
%     listok = {'T3-T5' 'F3-C3' 'C3-P3' 'P3-O1' 'P4-O2' 'C4-P4' 'Fp2-F4' 'T4-T6''F8-T4' 'Fp2-F8' 'Fp1-F7' 'F7-T3' }
%     %zone.label; %ENTER HERE DIFFERENT ORDER{'FT_R', 'TR','TL', 'FT_L','CL','FL', 'FTL','CR'}
%     listok  = {'Fp1-F7' 'Fp1-F3'  'F7-T3' 'T3-T5'  'F3-C3' 'C3-P3' 'P3-O1' 'P4-O2','C4-P4','F4-C4'   'T4-T6' 'F8-T4' 'Fp2-F4'  'Fp2-F8'}
%     listselected = get(handles.listbox_selectedzone,'string');
try
    listok = listselected;
end
idlist = [];
idlabel=[];
idzone =[];
for ilistzone = 1:numel(listok) % pour chaque zone dans la liste manuelle
    for izone = 1:numel(plotLst) %pour chaque zone dans data
        chzone = plotLst{izone}; % extraire les canaux de chaque zone dans data
        labelzone = label{izone}; %extraire le label de la zone dans data
        x = strmatch({labelzone} , {listok{ilistzone}}, 'exact'); %trouver la position du label des données dans la liste manuelle
        if ~isempty(x)
            idch = chMAT{izone}; %position de chaque canal dans la liste
            idlist = [idlist, idch]; %liste de la position de chaque canal dans la liste
            idzone =[idzone,izone, zeros(1,numel(idch)-1)];
            idlabel = [idlabel, {[label{izone}, sprintf('_%03.0f',ilistzone)]}];
        end
    end
end
idline = [find(idzone)-0.5,numel(idzone)+0.5];
%FIN

colorMap = jet(100);
c = max([max(MAT,[],'all'), abs(min(MAT,[],'all'))]);
cmin = -c + -c/10;
cmax = c + c/10;
cstep = (cmax-cmin)/100;
cf=cmin:cstep:cmax-cstep;
for i=1:size(MAT,1)
    for j = 1:size(MAT,2)
        colorMatrix(i,j) = sum(cf<MAT(i,j));
    end
end

tmp = [find(idzone), numel(idzone)];
tmphald = floor(find(idzone) + (tmp(2:end) - tmp(1:end-1))/2);
idmiddle = zeros(1,numel(idzone));
idmiddle(tmphald ) = idzone(find(idzone));

% Create custom node labels
label = strrep(label,'_',' '); %Modif KR
myLabel = cell(length(x));
colorlistline = zeros(size(MAT,1),3);
%line(numel(listok))
idcolor = 1;
%if get(handles.popup_connectogramlabel,'value')==1 %zone
    for i = 1:length(idzone)
        if idmiddle(i) %idzone(i)
            myLabel{i,1} =  label{idmiddle(i)};
            idcolor = idcolor + 1;
        else
            myLabel{i,1} =  ' ';
        end
        colorhomemade(i,:) =  colorlistline(idcolor,:);
    end

    for i = 1:length(idzone)
        colorhomemade(i,:) =  [0,0,0];
        % idcolor = idcolor + 1
    end

%elseif get(handles.popup_connectogramlabel,'value')==2 %channels
%   for i=1:numel(idzone)
%       myLabel{i,1} = [sprintf('%03.0f', idlist(i)), List(idlist(i),:)]
%   end
%end

figure; hold on
% Create custom colormap
hold on
x = MAT(idlist,idlist);
%         id0 = find(isnan(x))
%         x(id0)=0;
%         id0 = find(isinf(x))
%         x(id0)=0;
%         thresh = str2num(get(handles.edit_threshold,'string'));
%         x(x >  thresh) = 1;
%         x(x <= thresh) = 0;

xcolorMatrix = colorMatrix(idlist,idlist);

%circularGraph(MATsigch,'Label',labelch)
%save([savepath date '_plot.mat'],'x' 'myLabel' 'colorMatrix', 'colorMap','idlist')
myconnectogram(x,myLabel,xcolorMatrix,colorMap,idlist)
colormap('jet')
caxis([cmin cmax])
ax = gca;
ax.XAxis.Color = 'None';
ax.YAxis.Color = 'None';
ax.OuterPosition = [0.25 0.28 0.5 0.5];
colorbar('Position',[0.9 0.18 0.0305 0.7])