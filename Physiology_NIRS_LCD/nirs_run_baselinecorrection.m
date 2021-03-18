function out = nirs_run_baselinecorrection(job)
% EXACTLY SAME INPUTS AS NORMALIZATION FUNCTION
% SAME FUNCTION, EXCEPT THAT NO LOG10.
% The corrected data is recorded in a new .nir binary file.
%

%filename prefix
prefix = 'b'; %for "baseline correction"
DelPreviousData  = job.DelPreviousData;

for filenb=1:size(job.NIRSmat,1) %Loop over all subjects
    %Load NIRS.mat information
    NIRS = [];
    load(job.NIRSmat{filenb,1});
    [dir2,~,~] = fileparts(job.NIRSmat{filenb,1});
    %use last step of preprocessing
    lst = length(NIRS.Dt.fir.pp);
    rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
    NC = NIRS.Cf.H.C.N;
    fs = NIRS.Cf.dev.fs;
    ifile = 1; %Utiliser si on remets chaque stim normaliser dans des fichier diff�rent
    for f=1:numel(rDtp) %Loop over all files of a NIRS.mat
        d = fopen_NIR(rDtp{f},NC);
        dnorm = zeros(size(d));
        
        if isfield(job.normtype,'b_choiceglobal') %Global normalization
            if job.normtype.b_choiceglobal.m_choiceNan
                [dir1,fil1,~]=fileparts(rDtp{f});
                vmrk_path = fullfile(dir1,[fil1 '.vmrk']);
                [ind_dur_ch] = read_vmrk_find(vmrk_path,'bad_step');
                nch = size(d,1);
                nsample = size(d,2);
                noise = ind_dur_ch2mat(ind_dur_ch, nsample,nch)';
                idnoise = find(double(noise));
                dnan = d;
                if ~isempty(idnoise)
                    dnan(idnoise)=nan;
                end
                %figure;imagesc(dnan);
                for Idx = 1:NC
                    meansub = nanmean(dnan(Idx,:)); %Normalize whole blocks
                    dnorm(Idx,:) = d(Idx,:)-meansub; %BASELINE SUBSTRACTION
                end
            else
                for Idx = 1:NC
                    meansub = nanmean(d(Idx,:)); %Normalize whole blocks
                    dnorm(Idx,:) = d(Idx,:)-meansub; %BASELINE SUBSTRACTION
                end
            end
        elseif isfield(job.normtype,'b_choiceinternan') %Inter NAN normalization
            [dir1,fil1,~]=fileparts(rDtp{f});
            vmrk_path = fullfile(dir1,[fil1 '.vmrk']);
            [ind_dur_ch] = read_vmrk_find(vmrk_path,'bad_step');
            nch = size(d,1);
            nsample = size(d,2);
            noise = ind_dur_ch2mat(ind_dur_ch, nsample,nch)';
            idnoise = find(double(noise));
            dnan = d;
            if ~isempty(idnoise)
                dnan(idnoise)=nan;
            end
            for Idx = 1:NC
                segnanid = isnan(dnan(Idx,1:end-1))-isnan(dnan(Idx,2:end));
                idstartint = find( segnanid==+1);
                idstopint = find( segnanid==-1);
                if numel(idstartint) == numel( idstopint) %start stop
                    segnan = [[1,find( segnanid==+1) ];[find( segnanid==-1),size(d,2)]];
                elseif numel(idstartint) < numel(idstopint)
                    numel(segnanid);
                    segnan = [[1,idstartint ];[idstopint(1:end-1),size(d,2)]];
                elseif numel(idstartint) > numel(idstopint)
                    numel(segnanid);
                    segnan = [[1,idstartint(1:end-1) ];[idstopint,size(d,2)]];
                end
                for iseg = 1:size(segnan,2)
                    meansub = nanmean(dnan(Idx,segnan(1,iseg):segnan(2,iseg))); %Normalize using mean between nan periode
                    dnorm(Idx,segnan(1,iseg):segnan(2,iseg)) = d(Idx,segnan(1,iseg):segnan(2,iseg))-meansub;
                end
            end
        elseif isfield(job.normtype,'b_choicenormstim') %Stim normalization
            trigger = job.normtype.b_choicenormstim.trigger;
            pretime = round(fs*str2num(job.normtype.b_choicenormstim.pretime))-1;
            posttime = round(fs*str2num(job.normtype.b_choicenormstim.posttime));
            NIRS.Dt.fir.aux5bloc = NIRS.Dt.fir.aux5{f};
            aux5 = NIRS.Dt.fir.aux5{f};
            if trigger == 0 %If 0, consider every trigger
                indstim = aux5(:,2);
                ind = [indstim; size(d,2)];
                itrigger = 1;
            else
                indstim = [];
                itrigger = [];
                for itypestim = 1:numel(trigger)
                    idstim= aux5((aux5(:,1) == trigger(itypestim)),2);
                    indstim  = [indstim  ;idstim ];
                    ind = [1; idstim; size(d,2)];
                    itrigger =[itrigger, ones(1,numel(idstim)).*trigger(itypestim)];
                end
                
            end
            if job.normtype.b_choicenormstim.m_choiceNan
                [dir1,fil1,~]=fileparts(rDtp{f});
                vmrk_path = fullfile(dir1,[fil1 '.vmrk']);
                [ind_dur_ch] = read_vmrk_find(vmrk_path,'bad_step');
                nch = size(d,1);
                nsample = size(d,2);
                noise = ind_dur_ch2mat(ind_dur_ch, nsample,nch)';
                idnoise = find(double(noise));
                dnan = d;
                if ~isempty(idnoise)
                    dnan(idnoise)=nan;
                end
            else
                dnan = d;
            end
            
            for istim = 1:numel(indstim)
                
                if find((indstim(istim)-pretime) <= 0)
                    istart = 1;
                else
                    istart = indstim(istim)-pretime;
                end
                if find((indstim(istim)+posttime) > size(d,2))
                    istop = size(d,2);
                else
                    istop = indstim(istim)+posttime;
                end
                
                if job.normtype.b_choicenormstim.m_NormType==0      %I/Io Io = pretime
                    for Idx = 1:NC
                        meansub= mean(dnan(Idx,istart:indstim(istim)),2,'omitnan');
                        dnorm(Idx, istart:istop) = d(Idx,istart:istop)-meansub;
                    end
                elseif job.normtype.b_choicenormstim.m_NormType==1  %I/Io Io = pretime to posttime
                    for Idx = 1:NC
                        meansub= nanmean(dnan(Idx,istart:istop),2);
                        dnorm(Idx, istart:istop) = d(Idx,istart:istop)-meansub;
                    end
                end
            end
        end
        
        [dir1,fil1,ext1] = fileparts(rDtp{f});
        infilevmrk = fullfile(dir1,[fil1 '.vmrk']);
        infilevhdr = fullfile(dir1,[fil1 '.vhdr']);
        if ~exist(dir2,'dir'), mkdir(dir2); end
        outfile = fullfile(dir2,[prefix fil1 ext1]);
        outfilevmrk = fullfile(dir2,[prefix fil1  '.vmrk']);
        outfilevhdr = fullfile(dir2,[prefix fil1  '.vhdr']);
        fwrite_NIR(outfile,dnorm);
        copyfile(infilevmrk,outfilevmrk);
        
        
        try
            ChannelLabels = ConvertmlIDsrs2label(NIRS);
            SamplingInterval =floor(1000000/NIRS.Cf.dev.fs);
            nirs_boxy_write_vhdr(outfilevhdr,... %Output file
                outfile,... %DataFile
                outfilevmrk,... %MarkerFile,...
                'nirs_run_baselinecorrection',... %Function that created the header
                '',... %Channel Resolution
                '',... %Channel Units
                ChannelLabels,... %names given as a column of cells
                SamplingInterval,...
                size(d,2)); %SamplingInterval in microseconds
        catch
        end
        
        NIRS.Dt.fir.pp(lst+1).p{f,1} = outfile;
        fprintf('%s\n',outfile);
        
        
        if DelPreviousData
            delete(rDtp{f,1});
            delete(infilevmrk);
            delete(infilevhdr);
            try
                infileAC = fullfile(dir1,[fil1 'AC' '.nir']);
                delete(infileAC)
                infilePH = fullfile(dir1,[fil1 'PH' '.nir']);
                delete(infilePH)
            catch
            end
        end
        NIRS.Dt.fir.pp(lst+1).pre = 'BaselineCorrection';
        NIRS.Dt.fir.pp(lst+1).job = job;
    end
    
    save(fullfile(dir2,'NIRS.mat'),'NIRS');
    job.NIRSmat{1} =fullfile(dir2,'NIRS.mat');
    
end


out.NIRSmat = job.NIRSmat;


