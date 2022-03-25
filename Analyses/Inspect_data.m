%%From nirs_run_chcorrMAT of LIONirs toolbox%%
%Fichier excel avec le dossier des matrices, leur nom et le groupe
xlslistfile = 'C:\data\Malnutrition\Resting\NIRS\DocumentInfo\Subjectlist.xlsx';
savepath='C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Inspect\';
if ~isfolder(savepath)
    mkdir(savepath)
end

[~,~,ext] = fileparts([xlslistfile]);
if strcmp(ext,'.xlsx') | strcmp(ext,'.xls')
    [~, ~, xls] = xlsread([xlslistfile]);
elseif strcmp(ext,'.txt')
    [~, ~, xls] = readtxtfile_asxlsread([xlslistfile]);
end

NIRS = [];
for filenb = 2:size(xls,1) %do it one by one for the associate name
    id = filenb-1;
    name = xls{filenb,1};
    fprintf('Computing %s\n',name);
    NIRSfile = [xls{filenb,2} filesep xls{filenb,3}];
    load(NIRSfile);
    ML_new= [NIRS.Cf.H.C.id(2:3,:)',...
        ones(size(NIRS.Cf.H.C.id,2),1),...
        [ones(size(NIRS.Cf.H.C.id,2)/2,1);ones(size(NIRS.Cf.H.C.id,2)./2,1).*2]];
    lst = length(NIRS.Dt.fir.pp);
    rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
    NC = NIRS.Cf.H.C.N;

    for f = 1:size(rDtp,1) %Loop over all files of a NIRS.mat
        d1 = fopen_NIR(rDtp{f,1},NC); %open NIRS data
        dur{id,1}(1,f) = size(d1,2);
        mrk_type_arr = cellstr('bad_step');
        [dir1,fil1,~] = fileparts(rDtp{f});
        vmrk_path = fullfile(dir1,[fil1 '.vmrk']);
        [ind_dur_ch] = read_vmrk_find(vmrk_path,mrk_type_arr); %load the noise marker
        
        %nan and bad channels check
        A = sum(d1,2,'omitnan');
        A(A==0) = NaN;
        x = 1;
        for i = 1:size(d1,1) % for each channel
            if isnan(A(i)) %If all NaNs
                idxbadch(x,1) = i; %list of bad channels that are completely NaN in the data
                x = x + 1;
            else
            end
        end
        if ~exist('idxbadch','var')
            idxbadch = [];
        end
        str = num2str(idxbadch(:).');
        fprintf('Bad channels: %s\n',str)
        idgoodch = NIRS.Cf.H.C.ok(:,f); %list of good channels on the NIRS.mat
        idxnanch(:,1) = find(~idgoodch);
        str = num2str(idxnanch(:).');
        fprintf('NaN channels: %s\n',str)
        idxallbadch{id,1}(:,f) = union(idxbadch,idxnanch);
        pctbadch{id,1}(:,f) = numel(idxallbadch{id,1}(:,f))/NC*100;
        
        if ~isempty(ind_dur_ch)
            %rejected segments check
            idx = [];
            for n = 1:NC
                idx = find(ind_dur_ch(:,3) == n);
                  if ~isempty(idx)
                      timebad{id,1}(n,f) = sum(ind_dur_ch(idx,2));
                  else
                      timebad{id,1}(n,f) = 0;
                  end
                  idx = [];
            end
            %maxtimebad{id,f} = max(timebad{id,1}(:,f));
            meantimebad{id,f} = mean(timebad{id,1}(:,f));
            pctrejch{id,1}(:,f) = timebad{id,1}(:,f)./dur{id,1}(:,f)*100;
            pctrejch{id,1}(idxallbadch{id,1}(:,f),f) = NaN;
        else
            timebad{id,1}(:,f) = 0;
            pctrejch{id,1}(1:NC,f) = 0;
            pctrejch{id,1}(idxallbadch{id,1}(:,f),f) = NaN;
        end

       clearvars d1 mrk_type_arr dir1 fil1 vmrk_path ind_dur_ch A x i idxbadch idgoodch idxnanch idx n
    end
    clearvars id NIRSfile ML_new lst rDTP NC f
end

tbl = table(dur, pctbadch, meantimebad,'RowNames', xls(2:end,1));
writetable(tbl,[savepath 'inspectres.xls'],'WriteRowNames',true);