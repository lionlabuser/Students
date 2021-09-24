

clc;diskspace;
question=input('Do you want to delete matlabbatch before running the script?\n[yes = 1 / no = 0]:');
if question==1
    clear
end
generalpath='C:\Users\laura\OneDrive - Universite de Montreal\BACKup_Analyses\MART0m\Summary_badChannels\';
newfilterfolder='CardiacCrossCorr';
highpass='1'; %if none = 'none'
lowpass='3';

%% VARIABLES TO CREATE
part={ ...
    'BB002' 'BB003' 'BB004' 'BB008' 'BB009' 'BB010' ...'BB001'
    'BB011' 'BB012'         'BB014' 'BB015'   'BB016'   'BB017' 'BB018'...
    'BB019' 'BB020' ...  'BB013'
    'BB021' 'BB022' 'BB023' 'BB024' 'BB025' 'BB026' ...    'BB031''BB036'
    'BB027' 'BB028'         'BB030'  'BB032' 'BB033' 'BB034' 'BB035' 'BB037'  'BB040' ...  'BB038'
    'BB041' 'BB042' 'BB043' 'BB044' 'BB045' 'BB046' 'BB047' 'BB048' 'BB049' 'BB050' ...
    'BB051' 'BB052' 'BB053' 'BB055' 'BB056' 'BB057' 'BB058' 'BB059' 'BB060'  'BB062' 'BB063' 'BB065'   ...
    'BB066'   ...'BB067'
    };

for p=1:length(part)
    directory{p}=['C:\Users\laura\AnalysesNIRS\MART0m\' part{p} '\Segment\NormV2\'];
end

matlabbatch{1}.spm.tools.nirsHSJ.M_Utility.E_NIRSmatcreatenewbranch.NIRSmat = '<UNDEFINED>';
matlabbatch{1}.spm.tools.nirsHSJ.M_Utility.E_NIRSmatcreatenewbranch.e_NIRSmatdirnewbranch = newfilterfolder; %%HERE CHANGE
matlabbatch{1}.spm.tools.nirsHSJ.M_Utility.E_NIRSmatcreatenewbranch.m_newbranchcomponent = 1;

matlabbatch{2}.spm.tools.nirsHSJ.M_preprocessing.bpfilt.NIRSmat(1) = cfg_dep('dCONC: NIRS.mat', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','NIRSmat'));
matlabbatch{2}.spm.tools.nirsHSJ.M_preprocessing.bpfilt.DelPreviousData = 0;
matlabbatch{2}.spm.tools.nirsHSJ.M_preprocessing.bpfilt.lowcutfreq = lowpass;
matlabbatch{2}.spm.tools.nirsHSJ.M_preprocessing.bpfilt.highcutfreq = highpass;
matlabbatch{2}.spm.tools.nirsHSJ.M_preprocessing.bpfilt.filterorder = 4;
matlabbatch{2}.spm.tools.nirsHSJ.M_preprocessing.bpfilt.paddingsymfilter = 1;
matlabbatch{2}.spm.tools.nirsHSJ.M_preprocessing.bpfilt.interpolatebadfilter = 1;

%% run loop
for p= 1:length(part)
    
    %new folder + filter around 1-3hz
    matlabbatch{1}.spm.tools.nirsHSJ.M_Utility.E_NIRSmatcreatenewbranch.NIRSmat = {[directory{p} 'NIRS.mat']};
    jobid = cfg_util('initjob',matlabbatch);
    cfg_util('run', jobid)
    fprintf('New done: %s \n', part{p});
    
    %call DataQuality Check
    DataQualityCheck_CardiacCrossCorrelation(part{p},[directory{p} filesep newfilterfolder filesep],directory{p});
    
    %move .txt (summary file) and .mat (with correlation values) to a generic folder
    Qfiles=dir([directory{p} filesep newfilterfolder filesep part{p} '*']);
    for q=1:length(Qfiles)
        movefile(fullfile(Qfiles(q).folder,Qfiles(q).name),fullfile(generalpath,Qfiles(q).name));
    end
end