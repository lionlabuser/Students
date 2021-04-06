%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% parameter to adjust %%%%%%%%%%%%%%%%%%
fileXLS = 'C:\data\Malnutrition\Resting\NIRS\DocumentInfo\ReadDataRestMal.xlsx';

FolderAdjustment = 0;
AUXPHYSIOLOGY=1; %if you dont want, put 0
GLOBALPHYSIOLOGY=0; %if you dont want, put 0
WHOLEsPCA=0; %if you dont want, put 0
VIEWplotGLOBALphys=0; %if you dont want, put 0

%paths={'C:\data\Malnutrition\Resting\NIRS\G10115\DetectAuto\DetectManual\Filter0,01_0,08_dCONC\'};
%multimodalDirectory= {'C:\data\Malnutrition\Resting\NIRS\Multimodal\'};
%zoneFile ={'C:\data\Malnutrition\Resting\NIRS\Zones\Regressor HbO.zone'};

%You can automatically run multiple participants one after the other, using
%this format:
% paths={'C:\Data\ELAN\Martine_0m\BB016\Segment\NormV2\Test2021\'...
%        'C:\Data\ELAN\Martine_0m\BB017\Segment\NormV2\Test2021\'};
% multimodalDirectory={'C:\Data\ELAN\Martine_0m\Multimodal\BB016\'...
%        'C:\Data\ELAN\Martine_0m\Multimodal\BB017\'};
% zoneFile ={'C:\Data\ELAN\Martine_0m\BB016\GlobalZoneHBO.zone'...
%     'C:\Data\ELAN\Martine_0m\BB017\GlobalZoneHBO.zone'};


% ------------ parameters to adjust for AUXPHYSIOLOGY
jobA.outAUXfolder='filAUX'; %for nirs_run_filterAUX
jobA.copynirs=1; %for nirs_run_filterAUX
jobA.covariables='Sat,Resp'; %for nirs_run_GLM_regressAUX
jobA.e_NIRSmatdirnewbranch='SatResp'; %name of the new branch to create

%------------- parameters for GLOBALPHYSIOLOGY
jobG.trig = [0]; %for nirs_run_GlobalPhysio % 0 for resting OR trig number for task
jobG.globalavg = 0; %for nirs_run_GlobalPhysio
jobG.globalpca = 0; %for nirs_run_GlobalPhysio
jobG.spatialpca = 1; %for nirs_run_GlobalPhysio
jobG.e_NIRSmatdirnewbranch='SpatialPCA';%name of the new branch to create

%------------- parameters for WHOLEsPCA
jobW.goodPercent = 0.5;
jobW.e_NIRSmatdirnewbranch='SpatialPCA';%name of the new branch to create

%------------- parameters for VIEWplotGLOBALphys
PhysioLabels4FIGURES={'SatResp' 'SpatialPCA'}; %need to be in the SelectedFactors.mat

%%%%%%%%%%%%%%%%%%%%%%%%%  end  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%start of the automatic script
%for p=1:length(paths) %for each dataset. multiple paths can be enter in the paths variable.
%tempdirectory=paths{p};
[num,txt,raw] = xlsread(fileXLS);
for icol=1:size(raw,2)
    if strcmp(deblank(raw{1,icol}),'AnalysisFolder')
        ipath = icol;
    end
    if strcmp(deblank(raw{1,icol}),'AUX')
        iaux = icol;
    end
    if strcmp(deblank(raw{1,icol}),'ZoneAUX')
        izoneaux = icol;
    end
    if strcmp(deblank(raw{1,icol}),'Batch_Physio Fait (0) / À faire (1)')
        idoit  = icol;
    end
end

for isubject=2:size(raw,1)
    if isnan(raw{isubject,idoit})
        break
    end
    if raw{isubject,idoit} %afaire
        %try
        if FolderAdjustment==1
            
            jobF.NIRSmat = {fullfile([raw{isubject,ipath},'\DetectAuto\DetectManual\Filter0,01_0,08_dCONC\'],'NIRS.mat')}; %modify
            jobF.c_MultimodalPath.b_MultimodalPath_yes.e_MultimodalPath = cellstr(fileparts(raw{isubject,iaux})); %raw{isubject,iaux}
            
            nirs_run_NIRSmatdiradjust(jobF) %run script
        end
        
        if AUXPHYSIOLOGY==1
            %% filter AUX and extract regressed data
            jobA.NIRSmat = {fullfile([raw{isubject,ipath},'\DetectAuto\DetectManual\Filter0,01_0,08_dCONC\'],'NIRS.mat')}; %modify
            
            nirs_run_filterAUX(jobA); %run script
            nirs_run_GLM_regressAUX(jobA); %run script
            
            %% new branch + save corrected data
            
            jobA.m_newbranchcomponent=1; %1 if copy paste SelectedFactors.mat and CorrectionApply.mat
            jobA2=nirs_run_NIRSmatcreatenewbranch(jobA); %run script
            
            jobA2.globalmethod=jobA.e_NIRSmatdirnewbranch; %label that will be search into the SelectedFactors.mat PARCOMP variable
            jobA2.DelPreviousData=0;
            nirs_writeNIR_aftercorr(jobA2);%run script
            
            %suppress SelectedFactors.mat file in the new directory
            directory=fileparts(jobA2.NIRSmat{1});
            delete([directory filesep 'SelectedFactors.mat'])
        end
        
        if GLOBALPHYSIOLOGY==1
            
            %% run global physio + extract regressed data
            
            jobG.physzone = {raw{isubject,izoneaux}};
            jobG.NIRSmat = {fullfile([raw{isubject,ipath},'\DetectAuto\DetectManual\Filter0,01_0,08_dCONC\'],'NIRS.mat')}; %modify
            jobG.multimodalPATH = {[fileparts(raw{isubject,iaux}) '\AUXglobal']};
            
            nirs_run_GlobalPhysio(jobG) %run script
            
            %% new branch + save corrected data
            %create new branch
            jobG.m_newbranchcomponent=1;
            jobG2=nirs_run_NIRSmatcreatenewbranch(jobG);%run script
            
            %overwrite nirs data with the corrected data
            jobG2.globalmethod=jobA.e_NIRSmatdirnewbranch; %label that will be search into the SelectedFactors.mat PARCOMP variable
            jobG2.DelPreviousData=0;
            nirs_writeNIR_aftercorr(jobG2);
            
            %suppress SelectedFactors.mat file in the new directory
            directory=fileparts(jobG2.NIRSmat{1});
            delete([directory filesep 'SelectedFactors.mat'])
            clear job*
        end
        
        if WHOLEsPCA==1
            
            %% run SpatialPCA based on NAN segmentation + extract regressed data
            
            jobW.physzone = {raw{isubject,izoneaux}};
            jobW.NIRSmat = {fullfile([raw{isubject,ipath},'\DetectAuto\DetectManual\Filter0,01_0,08_dCONC\'],'NIRS.mat')}; %modify
            
            nirs_run_NANsegmentSpatialPCA(jobW) %run script
            
            %% new branch + save corrected data
            %create new branch
            jobW.m_newbranchcomponent=1;
            jobW2=nirs_run_NIRSmatcreatenewbranch(jobW);%run script
            
            %overwrite nirs data with the corrected data
            jobW2.globalmethod=jobW.e_NIRSmatdirnewbranch; %label that will be search into the SelectedFactors.mat PARCOMP variable
            jobW2.DelPreviousData=0;
            nirs_writeNIR_aftercorr(jobW2);
            
            %suppress SelectedFactors.mat file in the new directory
            directory=fileparts(jobW2.NIRSmat{1});
            delete([directory filesep 'SelectedFactors.mat'])
            %clear job*
        end
        
        if VIEWplotGLOBALphys==1
            
            physiolabels=PhysioLabels4FIGURES;
            
            %here is the script to visualize the Original data - Global
            %component - corrected data. nothing else to adapt
            load([tempdirectory 'SelectedFactors.mat'])
            for gc=1:length(physiolabels)
                idrow{gc}=find(contains({PARCOMP.label},physiolabels{gc}));
            end
            
            NC=size(PARCOMP(idrow{1}(1)).data,2)/2; %number of HBO channels
            NT=size(PARCOMP(idrow{1}(1)).data,1);
            for q=1:length(idrow{1})
                yy=[];
                figg=figure('units','normalized','outerposition',[0 0 1 1]);
                figg=tiledlayout(length(physiolabels),3,'TileSpacing','Compact','Padding','Compact');
                
                for gc=1:length(physiolabels)
                    tr=idrow{gc}(q);
                    nexttile;
                    plot(PARCOMP(tr).data(:,1:NC));
                    ylabel(physiolabels{gc},'fontweight','bold','FontSize',14);
                    if gc==1
                        title('HBO initial data');
                    end
                    yy=[yy ylim];
                    
                    nexttile;
                    plot(PARCOMP(tr).Xm(:,1:NC));
                    hold on
                    plot(mean(PARCOMP(tr).Xm(:,1:NC),2,'omitnan'),'Color','k','LineWidth',2);
                    yy=[yy ylim];
                    if gc==1
                        title('Global component (GC)');
                    end
                    
                    nexttile;
                    plot(PARCOMP(tr).dataCORR(:,1:NC));
                    if gc==1
                        title('Corrected data (GC subtracted)');
                    end
                    yy=[yy ylim];
                end
                
                %adjust the Y limits so that all are the same!
                for ir=1:(length(physiolabels)*3)
                    nexttile(ir);
                    xlim([0 NT])
                    ylim([min(yy) max(yy)])
                end
                
                if ~exist([tempdirectory 'GC_figure\'],'dir')
                    mkdir([tempdirectory 'GC_figure\']);
                end
                title(figg,['Block ' num2str(PARCOMP(tr).file) ])
                saveas(figg,[tempdirectory 'GC_figure\HBO_B' num2str(PARCOMP(tr).file) '.fig'])
                saveas(figg,[tempdirectory 'GC_figure\HBO_B'  num2str(PARCOMP(tr).file) '.png'])
                
                clear figg yy
                
            end
        end
        
        raw{isubject,idoit} = 0;
        %catch
        %disp(['Error subject ' raw{isubject,1},' fail'])
    disp([raw{isubject,1} ' done'])    
    end

end

try
    xlswrite(fileXLS,raw)
catch
    [filepath,name,ext] = fileparts(fileXLS);
    xlswrite(fullfile(filepath,['ERROR_OPENXLS',name,ext]),raw);
    disp(['ERROR XLS IS OPEN LOOK AT ', fullfile(filepath,['ERROR_OPENXLS',name,ext])])
end
