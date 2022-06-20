%Fichier excel avec le dossier des matrices, leur nom et le groupe
xlslistfile = 'C:\data\Malnutrition\Resting\NIRS\Analyses\CORRmatrice0,01_0,08\PCAPW\CorrPairC0,25 ExcY\Subjectlist N=54.xlsx';
fprintf('Computing Inspect_CO on %s\n', xlslistfile)

fisher = 1;
graph = 1;

if fisher == 1
    savepath = [fileparts(xlslistfile) '\graphmatrix_fisher\'];
    if ~isfolder(savepath)
        mkdir(savepath)
    end
elseif fisher == 0
    savepath = [fileparts(xlslistfile) '\graphmatrix_nofisher\'];
    if ~isfolder(savepath)
        mkdir(savepath)
    end
end

[~,~,ext] = fileparts([xlslistfile]);
if strcmp(ext,'.xlsx') | strcmp(ext,'.xls')
    [~, ~, xls] = xlsread([xlslistfile]);
elseif strcmp(ext,'.txt')
    [~, ~, xls] = readtxtfile_asxlsread([xlslistfile]);
end

%Load the matrices
groupeall = [];
list_subject = {};

for isubject = 2:size(xls,1)
    id = isubject-1;
    name = xls{isubject,2};
    idx = strfind(name,'_');
    name = name(1:idx(1)-1);
    fprintf('%s\n',name)
    list_subject{id,1} = name;
    groupeall = [groupeall; xls{isubject,4}];

    %Load the matrix
    load(fullfile(xls{isubject,1},[xls{isubject,2},'.mat']));
    
    if isnan(mean(meancorr,'all','omitnan')) || isnan(mean(matcorr,'all','omitnan'))
        fprintf('No data for %s\n',name)
        blocs(id,1) = 0;
        continue
    else
        blocs(id,1) = size(matcorr,3);
    end

    %Load the zone file
    load(fullfile(xls{isubject,1}, xls{isubject,3}),'-mat');

%         if find(strfind(xls{isubject,3},'v2'))
%             chlistfile = 'C:\data\Malnutrition\Resting\NIRS\channellist_zone_v2.txt';
%         else
%             chlistfile = 'C:\data\Malnutrition\Resting\NIRS\channellist_zone_v1.txt';
%         end

    if numel(ZoneList)>1
        if contains(ZoneList{1},'D0')
            Devicename = 'NIRx';
        else
            Devicename = 'ISS';
        end
    else
        Devicename = 'NIRx';
    end
    
    %Extract channel list from zone file
    idzone =[];
    chlist = {};
    for izone = 1:numel(zone.label)
        idch = zone.plot{izone};
        idzone = [idzone,izone, zeros(1,size(idch,1)-1)];
        for ichzone = 1:size(idch,1)
            ich = idch(ichzone,:);
            switch  Devicename %find label according to system
                case 'ISS'
                    strDet = SDDet2strboxy_ISS(ich(2));
                    strSrs = SDPairs2strboxy_ISS(ich(1));
                    chlist = [chlist;{[strDet, ' ',strSrs ]}];
                case 'NIRx'
                    strDet = SDDet2strboxy(ich(2));
                    strSrs = SDPairs2strboxy(ich(1));
                    chlist = [chlist;{[strDet, ' ',strSrs ]}];
                otherwise
                    strDet = SDDet2strboxy(ich(2));
                    strSrs = SDPairs2strboxy(ich(1));
                    chlist = [chlist;{[strDet, ' ',strSrs ]}];
            end
        end
    end

%     [~,~,ext] = fileparts([chlistfile]);
%     if strcmp(ext,'.xlsx') | strcmp(ext,'.xls')
%         [~, txt, ~] = xlsread([chlistfile]);
%     elseif strcmp(ext,'.txt')
%         [~, txt, ~] = readtxtfile_asxlsread([chlistfile]);
%     end
%     chlist_good = append(txt(:,1),' ',txt(:,2));

    for i = 1:numel(chlist)
        goodlist(i,:) = find(strcmp(ZoneList(:),chlist{i}));
    end
    
    matcorrgood = matcorr(goodlist,goodlist,:);
    meancorrgood = meancorr(goodlist,goodlist);
    MATallgood(:,:,id) = meancorrgood;
    %groupid(isubject)= DATA{idsubject(isubject)}.GR;

    if fisher == 1
        meancorrgood =  1/2*log((1+meancorrgood)./(1-meancorrgood));
        MATallgood(:,:,id) = meancorrgood;
        for b = 1:size(matcorrgood, 3)
            matcorrgood(:,:,b) = 1/2*log((1+matcorrgood(:,:,b))./(1-matcorrgood(:,:,b)));
        end
    end

    labelszone = zone.label;
    for i = 1:numel(labelszone)
        ilabel = labelszone{1,i};
        labelszone{1,i} = strrep(ilabel(1:end),'_',' ');
    end

    if graph == 1
        fig = figure;
        cmin = -max(abs(meancorrgood),[],'all');
        cmax = max(abs(meancorrgood),[],'all');
        clims = [cmin cmax];
        imagesc(meancorrgood,clims); %
        colorbar
        colormap(jet)
        set(gca,'xtick',1:size(meancorrgood,1)); %find(idzone)
        set(gca,'xticklabel', chlist); %labelszone
        %xtickangle(90);
        set(gca,'ytick', 1:size(meancorrgood,1)); %find(idzone)
        set(gca,'yticklabel', chlist); %labelszone
        ax = gca;
        ax.FontSize = 16;
        str = sprintf('%s Connectivity matrix',name);
        title(str)
        pbaspect([1 1 1])
        fig.WindowState = 'maximized';
        savefig([savepath name '_ConnMAT']);
        exportgraphics(gcf,[savepath name '_ConnMAT.png'])
        close
    end

    A = sum(meancorrgood,'omitnan');
    A(A==0) = NaN;
    idxnan = isnan(A);
    pnan = sum(idxnan)/numel(idxnan)*100;
    fprintf('%.2f percent of nan channels\n', pnan)

    idxneg = meancorrgood < 0;
    pneg = numel(nonzeros(idxneg))/numel(meancorrgood)*100;
    fprintf('%.2f percent of negative connections\n', pneg)
    idxnegnan = isnan(meancorrgood);
    idxneg = idxneg + idxnegnan/2;
    idxneg(1:(size(idxneg,1)+1):end) = 0.5;
    
    if graph == 1
        fig = figure;
        cmin = 0;
        cmax = 1;
        clims = [cmin cmax];
        imagesc(idxneg,clims);
        colorbar
        colormap(jet)
        set(gca,'xtick', 1:size(idxneg,1)); %find(idzone)
        set(gca,'xticklabel', chlist); %labelszone
        %xtickangle(90);
        set(gca,'ytick', 1:size(idxneg,1)); %find(idzone)
        set(gca,'yticklabel', chlist); %labelszone
        ax = gca;
        ax.FontSize = 16;
        str = sprintf('%s Negative connections',name);
        title(str)
        pbaspect([1 1 1])
        fig.WindowState = 'maximized';
        savefig([savepath name '_ConnNeg']);
        exportgraphics(gcf,[savepath name '_ConnNeg.png'])
        close
    end

    matstd = std(matcorrgood,0, 3,"omitnan");
    matstd(matstd==0) = NaN;
    stdallch = std(matstd,0,'all',"omitnan");
%     list_stdallch(id,:) = stdallch;
    meanallch = mean(matstd,'all',"omitnan");
    list_meanallch(id,:) = meanallch;
    matz = (matstd - meanallch)/stdallch;
    idxvar = abs(matz) >3.29;
    pvar = numel(nonzeros(idxvar))/numel(meancorrgood)*100;
    fprintf('%.2f percent of extreme variation between blocs\n\n', pvar)
    MATallstd(:,:,id) = matstd;
    MATallz(:,:,id) = matz;
    
    if graph == 1
        fig = figure;
        cmin = 0;
        cmax = 0.65; %display 2SD over the mean to isolate 5% more variable part
        clims = [cmin cmax];
        imagesc(matstd,clims);
        colorbar
        colormap(jet)
        set(gca,'xtick', 1:size(matstd,1)); %find(idzone)
        set(gca,'xticklabel', chlist); %labelszone
        %xtickangle(90);
        set(gca,'ytick', 1:size(matstd,1)); %find(idzone)
        set(gca,'yticklabel', chlist); %labelszone
        ax = gca;
        ax.FontSize = 16;
        str = sprintf('%s Variation(STD) of connection strength between blocks',name);
        title(str)
        pbaspect([1 1 1])
        fig.WindowState = 'maximized';
        savefig([savepath name '_ConnVar']);
        exportgraphics(gcf,[savepath name '_ConnVar.png'])
        close

        fig = figure;
        cmin = -3.29;
        cmax = 3.29;
        clims = [cmin cmax];
        imagesc(matz,clims); %,clims
        colorbar
        colormap(jet)
        set(gca,'xtick', 1:size(matz,1)); %find(idzone)
        set(gca,'xticklabel', chlist); %labelszone
        %xtickangle(90);
        set(gca,'ytick', 1:size(matz,1)); %find(idzone)
        set(gca,'yticklabel', chlist); %labelszone
        ax = gca;
        ax.FontSize = 16;
        str = sprintf('%s Distribution of the variation of connection strength between blocks',name);
        title(str)
        pbaspect([1 1 1])
        fig.WindowState = 'maximized';
        savefig([savepath name '_ConnZ']);
        exportgraphics(gcf,[savepath name '_ConnZ.png'])
        close
   end
end

clear fig ax cmin cmax clims A i b isubject id idx idxvar pvar idxnan pnan ...
    idxneg pneg idxnegnan name ext txt matcorr meancorr matcorrgood meancorrgood ...
    matstd matz meanallch stdallch ich ichzone idch ilabel izone str strDet strSrs

MATmeanG1 = mean(MATallgood(:,:,groupeall==1),3,'omitnan');
MATmeanG2 = mean(MATallgood(:,:,groupeall==2),3,'omitnan');
MATmeanG1G2 = MATmeanG1-MATmeanG2;
n_G1 = squeeze(sum(~isnan(MATmeanG1),3));
n_G2 = squeeze(sum(~isnan(MATmeanG2),3));

if graph == 1
    fig = figure;
    cmin = -1.6; %-max(abs(MATmeanG1),[],'all');
    cmax = 1.6; %max(abs(MATmeanG1),[],'all');
    clims = [cmin cmax];
    imagesc(MATmeanG1,clims); %
    colorbar
    colormap(jet)
    set(gca,'xtick', find(idzone)); %1:size(MATmeanG1,1)
    set(gca,'xticklabel', labelszone); %chlist
    %xtickangle(90);
    set(gca,'ytick', find(idzone)); %1:size(MATmeanG1,1)
    set(gca,'yticklabel', labelszone); %chlist
    ax = gca;
    ax.FontSize = 20;
    title('Mean G1 Connectivity Matrix')
    pbaspect([1 1 1])
    fig.WindowState = 'maximized';
    savefig([savepath 'ConnMATmeanG1']);
    exportgraphics(gcf,[savepath 'ConnMATmeanG1.png'])
    close

    fig = figure;
    cmin = -1.6; %-max(abs(MATmeanG2),[],'all');
    cmax = 1.6; %max(abs(MATmeanG2),[],'all');
    clims = [cmin cmax];
    imagesc(MATmeanG2,clims); %
    colorbar
    colormap(jet)
    set(gca,'xtick', find(idzone)); %1:size(MATmeanG2,1)
    set(gca,'xticklabel', labelszone); %chlist
    %xtickangle(90);
    set(gca,'ytick', find(idzone)); %1:size(MATmeanG2,1)
    set(gca,'yticklabel', labelszone); %chlist
    ax = gca;
    ax.FontSize = 20;
    title('Mean G2 Connectivity Matrix')
    pbaspect([1 1 1])
    fig.WindowState = 'maximized';
    savefig([savepath 'ConnMATmeanG2']);
    exportgraphics(gcf,[savepath 'ConnMATmeanG2.png'])
    close

    fig = figure;
    cmin = -0.5; %-max(abs(MATmeanG1G2),[],'all');
    cmax = 0.5; %max(abs(MATmeanG1G2),[],'all');
    clims = [cmin cmax];
    imagesc(MATmeanG1G2,clims); %
    colorbar
    colormap(jet)
    set(gca,'xtick', find(idzone)); %1:size(MATmeanG1,1)
    set(gca,'xticklabel', labelszone); %chlist
    %xtickangle(90);
    set(gca,'ytick', find(idzone)); %1:size(MATmeanG1,1)
    set(gca,'yticklabel', labelszone); %chlist
    ax = gca;
    ax.FontSize = 20;
    title('G1-G2 Connectivity Matrix')
    pbaspect([1 1 1])
    fig.WindowState = 'maximized';
    savefig([savepath 'ConnMATmeanG1G2']);
    exportgraphics(gcf,[savepath 'ConnMATmeanG1G2.png'])
    close

    fig = figure;
    cmin = 0.35;
    cmax = 0.55; %display 2SD over the mean to isolate 5% more variable part
    clims = [cmin cmax];
    imagesc(mean(MATallstd,3,'omitnan'),clims);
    colorbar
    colormap(jet)
    set(gca,'xtick', find(idzone)); %find(idzone)
    set(gca,'xticklabel', labelszone); %labelszone
    %xtickangle(90);
    set(gca,'ytick', find(idzone)); %find(idzone)
    set(gca,'yticklabel', labelszone); %labelszone
    ax = gca;
    ax.FontSize = 20;
    str = sprintf('Mean variation(STD) of connection strength between blocks');
    title(str)
    pbaspect([1 1 1])
    fig.WindowState = 'maximized';
    savefig([savepath 'ConnMeanVar']);
    exportgraphics(gcf,[savepath 'ConnMeanVar.png'])
    close
end

% % %Determine clims for SD graph for the whole sample
%  meanallpart = mean(list_meanallch);
%  stdallpart = mean(list_stdallch);
%  zallpart = (list_meanallch - meanallpart)/stdallpart;

[~,I] = sort(list_meanallch,'descend');
varall(:,1) = list_subject(I);
varall(:,2) = num2cell(list_meanallch(I));

tblblocs = array2table([blocs groupeall], 'RowNames', list_subject, 'VariableNames',{'NbBlocs', 'Group'});
writetable(tblblocs,[fileparts(xlslistfile) '\GoodBlocs.xls'],'WriteRowNames',true);

clear fig ax cmin cmax clims
clear all