function out = nirs_run_physio_PCA(job)

NIRS = [];
load(job.NIRSmat{1});
ML_new= [NIRS.Cf.H.C.id(2:3,:)',...
    ones(size(NIRS.Cf.H.C.id,2),1),...
    [ones(size(NIRS.Cf.H.C.id,2)/2,1);ones(size(NIRS.Cf.H.C.id,2)./2,1).*2]];
lst = length(NIRS.Dt.fir.pp);
rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
NC = NIRS.Cf.H.C.N;
fs = NIRS.Cf.dev.fs;
idpart = job.ID;
fprintf('Computing %s filter for %s\n',job.e_NIRSmatdirnewbranch, idpart)

%check if SelectedFactors file exist or create one
[nirsPATH,~,~] = fileparts(job.NIRSmat{1});
if exist(fullfile(nirsPATH,'SelectedFactors.mat'),'file')
    load(fullfile(nirsPATH,'SelectedFactors.mat'))
    if any(contains({PARCOMP.type}, 'PCA')) && any(contains({PARCOMP.label}, ['ExtPhysio_' job.e_NIRSmatdirnewbranch '_Bl01']))
        answer = questdlg('PCA already computed, what do you want to do?','Warning','Overwrite','Skip Participant','Skip Participant');
        switch answer
            case 'Overwrite'
                corroverwrite = 1;
                Prow = find(contains({PARCOMP.label}, ['ExtPhysio_' job.e_NIRSmatdirnewbranch '_Bl01']));
                disp(['The old PCA filter will be overwritten'])
            case 'Skip Participant'
                return
        end
    else
        Prow = length(PARCOMP) + 1; %number of components computed for each participants (1 per block) A REMETTRE
    end
else
    PARCOMP = struct;
    Prow = 1;
end

IDrows=[]; %Id of rows for the current script in the PARCOMP - will be use later to create a graph output regarding the betas

for f = 1:size(rDtp,1) %For each block - Loop over all files of a NIRS.mat
    d = fopen_NIR(rDtp{f,1},NC); %open NIRS data
    dur = NIRS.Dt.fir.sizebloc{f,1};
    t = [0:1/fs:dur/fs]';
    t = t(1:end-1);
    nanch = NIRS.Cf.H.C.ok(:,f);
    mrk_type_arr = cellstr('bad_step');
    [dir1,fil1,~] = fileparts(rDtp{f});
    vmrk_path = fullfile(dir1,[fil1 '.vmrk']);
    [ind_dur_ch] = read_vmrk_find(vmrk_path,mrk_type_arr); %load the noise marker
    MA = zeros(size(d));
    if ~isempty (ind_dur_ch)
        for j = 1:NC
            idx = find(ind_dur_ch(:,3)==j);
            if ~isempty(idx)
                for a = 1:numel(idx)
                    MA(j,ind_dur_ch(idx(a),1):(ind_dur_ch(idx(a),1)+ind_dur_ch(idx(a),2)-1)) = 1;
                end
            end
            idx = [];
        end
    else
        fprintf('No artifact found for %s\n',job.ID)
    end

    %% interpolate bad intervals from Brain Analyzer toolbox%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    d1 = d';
    list = isnan(d1);

    if any(list(:))
        try
            for j = 1:size(d1,2) %for each channel
                if(all(list(:,j))) %if all values of the channel are 1
                    d1(list(:,j),j) = job.ifFailReplaceWith;
                elseif any(list(:,j)) %if some values are 1
                    % interpolation
                    l = list(:,j); %extract nan tp for channel j
                    d1(l,j) = interp1(t(~l), d1(~l,j), t(l),job.type);
                    %data(i).data = d;
                end
            end

        catch
            % just replace with white noise
            d1(list) = job.ifFailReplaceWith;
        end

    end

    % repeat to get the edges using nearest
    list = isnan(d1);

    if any(list(:))
        try
            for j = 1:size(d1,2)
                if(all(list(:,j)))
                    d1(list(:,j),j) = job.ifFailReplaceWith;
                elseif any(list(:,j))
                    % interpolation
                    l = list(:,j);
                    d1(l,j) = interp1(t(~l), d1(~l,j), t(l),'nearest','extrap');
                    %data(i).data = d;
                end
            end

        catch
            % just replace with white noise
            d1(list) = job.ifFailReplaceWith;
        end

    end

    clear list

    %% PCA from Brain Analyzer toolbox%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if(job.keep_components) %ajuster selon nos données
        comp = data(i);
        comp.auxillary('eigenvalues') = [];
        data(i).auxillary('PCAfilter') = comp;
    end

    types = unique(NIRS.Cf.dev.wl)'; %extract wavelengths
    for tI = 1:length(types) %for each wavelength

        if(tI>1 & ~job.splittypes)
            continue;
        end

        if(job.splittypes)
            if(~iscell(types))
                types = num2cell(types);
            end
            list = find(ML_new(:,4) == tI); %find selected wavelength channels
        else
            list = 1:size(d1,2); %take all channels
        end

        % resample data
        d2 = d1(:,list);

        % remove mean
        m = mean(d2,1);
        d2 = bsxfun(@minus, d2, m);

        % svd
        [u, s, v] = svd(d2,'econ');
        s = diag(s);

        if(job.ncomp>0 & job.ncomp<1)
            % if fraction, then remove that % of the noise
            n = max(find(cumsum(s)/sum(s)<=job.ncomp));
            disp(['Removing first ' num2str(n) ' components (' num2str(sum(s(1:n))/sum(s)*100) '%)']);

            if(n==0)
                disp(['     lowest component = ' num2str(sum(s(1))/sum(s)*100) '%)']);
            end
        else
            n = job.ncomp;
        end

        if(job.keep_components) % à ajuster à notre type de données
            s2 = s;
            s2(n+1:end) = 0;
            comp = data(i).auxillary('PCAfilter');
            comp.data(:,list) = u*diag(s2)*v';
            eig = comp.auxillary('eigenvalues');
            eig(:,end+1) = s;
            data(i).auxillary('PCAfilter') = comp;
        end

        % remove n components
        s(1:n) = 0;
        d2 = u*diag(s)*v';

        % add mean back
        d2 = bsxfun(@plus, d2, m);

        % save new data
        d3(:,list) = d2;

        if job.runPW == 0
            d5 = d3;
        end
    end

    %% AR prewhitening from Brain Analyzer toolbox%%
    if job.runPW == 1
    
        medians = median(d3);
        d4 = bsxfun( @minus , d3 , medians );
    
        %goodCH = [];
        for j = 1:size(d4, 2)
            %goodCH = [goodCH j];
            y   = d4(:,j);
    
            if(isstr(job.modelorder))
                p = fs*str2num(job.modelorder(1:strfind(job.modelorder,'x')-1));
            else
                p = job.modelorder;
            end
    
            [~,~,~,ymoco] = nirs.math.robust_ari1_fit(y, round(p), job.tune);
    
            d5(:,j) = ymoco; %data(i).data(:,j) = ymoco;
        end
        d5 = bsxfun(@plus,bsxfun(@minus,d5,median(d5)),medians);
    end

    %% write SELECTED FACTORS new info
    PARCOMP(Prow).file= f;
    PARCOMP(Prow).filestr =  sprintf('Bl%02.0f', f);
    PARCOMP(Prow).label= ['ExtPhysio_' job.e_NIRSmatdirnewbranch '_' PARCOMP(Prow).filestr]; %[cov.labels{~(1:length(cov.labels)==cov.ConstantID)}]
    if job.runPW == 1
        PARCOMP(Prow).type = 'PCAPW';
    else
        PARCOMP(Prow).type = 'PCA';
    end
    PARCOMP(Prow).beta = []; %tmpbeta; %beta
    PARCOMP(Prow).std = []; %tmpErrorVariance;
    PARCOMP(Prow).AUX = []; %données AUX
    PARCOMP(Prow).data = d'; %données NIRS
    PARCOMP(Prow).Xm = d' - d5; %résiduels entre les données originales et corrigées
    PARCOMP(Prow).dataCORR = d5; %données NIRS régressées
    PARCOMP(Prow).indt = [1 :size(d,2)]; %indice de temps.
    PARCOMP(Prow).listgood = []; %goodCH;
    PARCOMP(Prow).module  = lst;
    PARCOMP(Prow).modulestr = NIRS.Dt.fir.pp(lst).pre;
    PARCOMP(Prow).ComponentToKeep = 1;
    PARCOMP(Prow).idreg = 1;
    PARCOMP(Prow).topo = []; %tmpbeta(PARCOMP(Prow).ComponentToKeep,:);

    IDrows = [IDrows Prow]; %nb of blocks
    Prow = Prow + 1;

    disp(['block ' num2str(f) ' done'] );

    %% FIGURES
    yy=[];
    fig = figure; %('units','normalized','outerposition',[0 0 1 1]);
    fig.WindowState = 'maximized';
    if job.runPW
        fig = tiledlayout(4,1,'TileSpacing','Compact','Padding','Compact');
    else
        fig = tiledlayout(3,1,'TileSpacing','Compact','Padding','Compact');
    end
    ylabel(fig, ['Block ' num2str(PARCOMP(IDrows(f)).file)],'fontweight','bold','FontSize',12)
    title(fig, [job.ID ' ' job.e_NIRSmatdirnewbranch ' correction'])

    nexttile;
    plot(d1(:,1:(NC/2)))
    %set(gca,'xtick',t);
    yy = [yy ylim];
    hold on
    plot(mean(d1(:,1:(NC/2)),2,'omitnan'),'Color','k','LineWidth',2);
    title('HBO interpolated data');

    nexttile;
    plot(d3(:,1:(NC/2)));
    %set(gca,'xtick',t);
    %yy = [yy ylim];
    hold on
    plot(mean(d3(:,1:(NC/2)),2,'omitnan'),'Color','k','LineWidth',2);
    title('PCA corrected data');

    if job.runPW
        nexttile;
        plot(d5(:,1:(NC/2)));
        %set(gca,'xtick',t);
        %yy=[yy ylim];
        hold on
        plot(mean(d5(:,1:(NC/2)),2,'omitnan'),'Color','k','LineWidth',2);
        title('Prewhitened data');
    end

    nexttile;
    d5nan = d5;
    d5nan(logical(MA')) = NaN;
    d5nan(:,~nanch) = NaN;
    plot(d5nan(:,1:(NC/2)));
    %set(gca,'xtick',t);
    %yy=[yy ylim];
    hold on
    plot(mean(d5(:,1:(NC/2)),2,'omitnan'),'Color','k','LineWidth',2);
    title('Final data');

    %adjust the Y limits so that all are the same!
    for ir = 1:size(fig.Children,1)
        nexttile(ir);
        xlim([0 size(d1,1)]) %mettre la même limite des axes x et y à toutes les tuiles
        ylim(yy)
    end

    %save fig
    saveas(fig,fullfile(nirsPATH,['PhysioCorr_b' num2str(f) '.png']))
    saveas(fig,fullfile(nirsPATH,['PhysioCorr_b' num2str(f) '.fig']))
    close; clear fig
end

%% save nirs mat and PARCOMP

fprintf('PCA filter done computing for %s\n',job.ID)
save(job.NIRSmat{1},'NIRS');
save(fullfile(nirsPATH,'SelectedFactors.mat'),'PARCOMP');

out.NIRSmat = job.NIRSmat;

end