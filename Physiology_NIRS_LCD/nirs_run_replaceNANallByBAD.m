function out=nirs_run_replaceNANallByBAD(job)
% job structure needs:
% - job.NIRSmat={'C:....\NIRS.mat'};
% - job.DelPreviousData= 0 or 1 : 0 if KEEP the previous
% WHAT IT DOES: it will search for the nirs file BEFORE correction (from
% AUX or SpatialPCA) and BEFORE the Nullify Bad interval Step (therefore
% the file before the use of the prefix <cnull>). It will copy paste the
% bad yellow identification (where there were NaN values for all channels)
% + replace the NaN values that are present for at least 80% of channels by
% the data point value before to nan.

DelPreviousData  = job.DelPreviousData;
prefix ='b';
NIRS = [];
load(job.NIRSmat{1});
%use last step of operation
lst = length(NIRS.Dt.fir.pp);
rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
NC = NIRS.Cf.H.C.N;
fs = NIRS.Cf.dev.fs;

%Verify job.NIRS.mat location, dir2 will be the new location
[dir2,~,~] = fileparts(job.NIRSmat{1});

for f=1:numel(rDtp) %Loop over all files of a NIRS.mat
    
    [dir1,fil1,ext1] = fileparts(rDtp{f});
    
    dirt=split(dir1,filesep);
    findprefix='cnull'; %5caracters
    if exist(fullfile(dir1,[fil1(6:end) ext1]))
        
        infilevmrk = fullfile(dir1,[fil1(6:end) '.vmrk']);
        infilevhdr = fullfile(dir1,[fil1(6:end) '.vhdr']);
    elseif exist(fullfile(dirt{1:end-1},[fil1(6:end) ext1]))
        infilevmrk = fullfile(dirt{1:end-1},[fil1(6:end) '.vmrk']);
        infilevhdr = fullfile(dirt{1:end-1},[fil1(6:end) '.vhdr']);
    elseif exist(fullfile(dirt{1:end-2},[fil1(6:end) ext1]))
        infilevmrk = fullfile(dirt{1:end-2},[fil1(6:end) '.vmrk']);
        infilevhdr = fullfile(dirt{1:end-2},[fil1(6:end) '.vhdr']);
    else
        error('There are no files without %s: %s\nIn one of those directories: \n- %s\n- %s\n- %s\n',...
            findprefix,[fil1(6:end) ext1],dir1,fullfile(dirt{1:end-1}),fullfile(dirt{1:end-2}));
    end
    outfile = fullfile(dir2,[prefix fil1 ext1]);
    outfilevmrk = fullfile(dir2,[prefix fil1 '.vmrk']);
    outfilevhdr = fullfile(dir2,[prefix fil1 '.vhdr']);
    

    
    d = fopen_NIR(rDtp{f},NC);
 if size(d,1)==NC
        %do nothing
        repermute=0;
    elseif size(d,2)==NC
        d=permute(d,[2 1]); %organize the data according to Dim1=channels and Dim2 = time
        repermute=1;
    else
        error ('The size of your data does not match the number of channels...')
    end
    
        badch = NIRS.Cf.H.C.ok(:,f);
        if badch==0 %if all channels are bad, switch to next block!
            continue
        end
    
    nanid=isnan(d); %identifies with a 1 all data positions that ARE  a NAN value
    nansum=find(sum(nanid)>(NC*.8)); %for each time point, 
                     % check if you have more valid channels than 
                     % the cut-off percentage. If so, return a 1 value (if
                     % not, return a 0 value)
for nn=nansum
    if nansum(nn)==1
        d(:,nansum(nn))=mean(d,2,'omitnan');
    else
    d(:,nansum(nn))=d(:,nansum(nn)-1);
    end
end
npoints=size(d,2);
if repermute==1
    d=permute(d,[2 1]);
end
    fwrite_NIR(outfile,d);
    
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
            'replaceNANallByBAD',... %Function that created the header
            '',... %Channel Resolution
            '',... %Channel Units
            ChannelLabels,... %names given as a column of cells
            SamplingInterval,...
            npoints); %SamplingInterval in microseconds
    catch
        disp('problem with writing VHDR')
    end
    
    if DelPreviousData
        try  delete(rDtp{f}); catch;   end
        try  delete(fullfile(dir1,[fil1 '.vmrk']));catch;   end
        try  delete(fullfile(dir1,[fil1 '.vhdr']));catch;   end
    end
    %add outfile name to NIRS
    if f == 1
        NIRS.Dt.fir.pp(lst+1).pre = 'Replace NAN by badID';
        NIRS.Dt.fir.pp(lst+1).job = job;
    end
    NIRS.Dt.fir.pp(lst+1).p{f,1} = outfile;
    clear d
end
save(fullfile(dir2,'NIRS.mat'),'NIRS');
out =fullfile(dir2,'NIRS.mat');
    

end