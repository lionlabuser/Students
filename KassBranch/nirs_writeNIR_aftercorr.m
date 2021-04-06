function out = nirs_writeNIR_aftercorr(job)
DelPreviousData  = job.DelPreviousData;
prefix = 'c'; %for "corrected"

corr=job.globalmethod;

load(job.NIRSmat{1});

%use last step of preprocessing
lst = length(NIRS.Dt.fir.pp);
rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
NC = NIRS.Cf.H.C.N;
fs = NIRS.Cf.dev.fs;

%Verify job.NIRS.mat location, dir2 will be the new location
[dir2,~,~] = fileparts(job.NIRSmat{1});

% Load Selected factors
load(fullfile(dir2,'SelectedFactors.mat'),'PARCOMP');

Grows=find(contains({PARCOMP.label},corr)); %rows where the selected global regression is localized
validblocks=[PARCOMP(Grows).file];


for f=1:size(rDtp,1) %Loop over all files of a NIRS.mat
    
    [dir1,fil1,ext1] = fileparts(rDtp{f});
    infilevmrk = fullfile(dir1,[fil1 '.vmrk']);
    infilevhdr = fullfile(dir1,[fil1 '.vhdr']);
    
    outfile = fullfile(dir2,[prefix fil1 ext1]);
    outfilevmrk = fullfile(dir2,[prefix fil1 '.vmrk']);
    outfilevhdr = fullfile(dir2,[prefix fil1 '.vhdr']);
    
    d = fopen_NIR(rDtp{f,1},NC);
    
    if any(validblocks==f) %if the block is identified as good
        newDATA=PARCOMP(Grows(find(validblocks==f))).dataCORR;
        if NC==size(newDATA,1)
            %do nothing
        elseif NC==size(newDATA,2)
            newDATA=permute(newDATA,[2 1]);
        else
            error('Corrected data doesnt have the same number of channels as the original data. Please check.')
        end
        
    else %all channels of the block have been previously marked as bad
        newDATA=nan(size(d));
        
    end
    
    fwrite_NIR(outfile,newDATA);
    
    fprintf('%s\n',outfile);
    
    %write new .vmrk file
    try
        copyfile(infilevmrk,outfilevmrk);
    catch
        disp('problem with writing VMRK')
    end
    try
        copyfile(infilevhdr,outfilevhdr);
        ChannelLabels = ConvertmlIDsrs2label(NIRS);
        SamplingInterval =floor(1000000/NIRS.Cf.dev.fs);
        nirs_boxy_write_vhdr(outfilevhdr,... %Output file
            outfile,... %DataFile
            outfilevmrk,... %MarkerFile,...
            'nirs_writeNIR_aftercorr',... %Function that created the header
            '',... %Channel Resolution
            '',... %Channel Units
            ChannelLabels,... %names given as a column of cells
            SamplingInterval,...
            size(newDATA,2)); %SamplingInterval in microseconds
    catch
        disp('problem with writing VHDR')
    end
    if DelPreviousData
        try  delete(rDtp{f,1}); catch;   end
        delete(infilevmrk)
        delete(infilevhdr)
    end
    %add outfile name to NIRS
    if f == 1
        NIRS.Dt.fir.pp(lst+1).pre = ['Corrected with ' corr];
        NIRS.Dt.fir.pp(lst+1).job = job;
    end
    NIRS.Dt.fir.pp(lst+1).p{f,1} = outfile;
    clear newDATA
end
save(fullfile(dir2,'NIRS.mat'),'NIRS');
out =fullfile(dir2,'NIRS.mat');

end
