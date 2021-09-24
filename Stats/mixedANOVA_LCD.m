%% final anova script
outputdirectory='C:\Users\laura\OneDrive - Universite de Montreal\CYCLES SUP\THÈSE\ARTICLE1_MARTINE-0M\FinalResults\';
figdirectory='C:\Users\laura\OneDrive - Universite de Montreal\CYCLES SUP\THÈSE\ARTICLE1_MARTINE-0M\FinalResults\Figures\';

 matfile='C:\Users\laura\AnalysesNIRS\MART0m\0FINAL\2021-08-21_DATAmat_v3detrend_AUX5z3.29.mat';  tag='AUX5z3.29';
load(matfile)
disp(['Matfile: ' matfile])

cond={'FR' 'AL' 'HE'};
longcond={'French' 'German' 'Hebrew'};
hemis={'LH' 'RH'};

partID={DATA.part}'; %participant ID
groupID=[DATA.group]+1; %group identification (here as 1-2-3)
group={'CTL' 'XPH' 'XPA'}; 
groupfig=renamecats(categorical(groupID),{'Control' 'XP_H_e_b_r_e_w' 'XP_G_e_r_m_a_n'});


%% create database
clc

% montage specifications
roinames={'Fro'  'Temp' 'TPJ' 'Post'};
roinameslong={'Frontal'  'Temporal' 'Temporo-Parietal' 'Posterior' };
vROI='v1';
roi{1}=[ 1 41 ;2 42 ;12 54 ; 13 53 ; 6 49 ; 4 50 ]; 
roi{2}=[8 46 ; 16 40  ;18 37 ; 19 36];    ...
roi{3}=[17 38; 20 35 ;21 34 ; 22 33];
roi{4}=[26 29 ; 25 30 ; 24 31 ; 23 32]; %channel ID

% create A MATRICE #participants X #conditions X #roi X #hemispheres
clear d t* within*
for p=1:length(DATA) %PARTICIPANTS
    for c=1:3 %CONDITIONS 
        for r=1:length(roi) %ROI
            for h=1:2 %HEMISPHERE
                d.hbo.activ(p,r,h,c)=['ENTER THE "WAY TO ACCESS YOUR DATA HERE'];  % TO CHANGE
                %E.G. d.hbo.activ(p,r,h,c)=DATA(p).(cond{c})(r,h)
            end
        end
    end
end
 
fprintf('Number of included participants:\n...HBO= %s \n\n',num2str(sum(sum(sum(sum(isnan(d.hbo.activ),4),3),2)==0)))

  d.hbo.activ(sum(sum(sum(isnan(d.hbo.activ),4),3),2)>0,:,:,:)=NaN; %if one 
  %participant have a missing data on one of its repeated measures, it will 
  %not be included in an repeated measure anova. Therefore I change all its
  %values to NAN, so that graphics reflect really the data we used in the
  %anova

  %% construct the table for the ANOVA (as a SPSS format)
thbo=table(partID,'VariableNames',{'part'}); %create table with a PARTICIPANT column
thbo.group=categorical(groupID); %inter-subject factors FOR CATEGORICAL VARIABLE (EX. GROUP)
thbo.group=renamecats(thbo.group, group);
thbo.group=reordercats(thbo.group,group);
thbo.age=partAGE; %inter-subject factors FOR CONTINUOUS VARIABLE (EX. AGE)
[~,numISF]=size(thbo); %number of rows associated to inter-subject factors (ISF)
x=1;
for c=1:length(cond)
    for r=1:length(roi)
        for h=1:length(hemis)
            thbo.([cond{c} hemis{h} roinames{r}])=d.hbo.activ(:,r,h,c); % one column = data from one repeated measure
            within_table{x,1}=cond{c}; %additionnal table to save what levels of the repeated measure for each column
            within_table{x,2}=hemis{h};
            within_table{x,3}=roinames{r};
            x=x+1;
        end
    end
end

within_table = cell2table(within_table,'VariableNames', {'Cond' 'Hemis' 'ROI'});
within_table.Cond=categorical(within_table.Cond);
within_table.Hemis=categorical(within_table.Hemis);
within_table.ROI=categorical(within_table.ROI);
vars=thbo.Properties.VariableNames;

%% ANOVA GROUP X COND X HEMIS X ROI
rm1a = fitrm(thbo,...
    [ vars{numISF+1} '-' vars{end}  '~group'],... it means: from the first column of repeated measures to the last column = dependent variables. ~ separate DEP var and INDEP var. group is my INter-subject var of interest
    'WithinDesign', within_table);

tbl.roi = ranova(rm1a,'WithinModel', 'Cond*Hemis*ROI');
tbl.roi =addsig(tbl.roi );
disp(tbl.roi)

%% Figures
%roi x group x cond+hemis
 plotgraph(d.hbo.activ,groupfig,{'ROIs' roinameslong},{'Hemispheres' hemis},{'Conditions' longcond},figdirectory,vROI);
% group x cond+hemis
 plotgraphmain(d.hbo.activ,groupfig,{'ROIs' roinameslong},{'Hemispheres' hemis},{'Conditions' longcond},figdirectory,vROI);


%% main anova by group
tmp=thbo(groupID==1 | groupID==3,:); 
tmp.group=removecats(tmp.group);
rmG1 = fitrm(tmp,...
    [ vars{numISF+1} '-' vars{end}  '~group'],...
    'WithinDesign', within_table);

tbl.roi = ranova(rmG1,'WithinModel', 'Cond*Hemis*ROI');
tbl.roi =addsig(tbl.roi );
disp(tbl.roi)
tmp=multcompare(rmG1,'ROI'); newtmp=addsig(tmp);  disp(newtmp(newtmp.Difference>0,:));


%% gr cond hemis effect
tripleinter(rm1a )

%% ROI differences   
tmp=multcompare(rm1a,'ROI');
newtmp=addsig(tmp);  disp(newtmp(newtmp.Difference>0,:));

%% ROI x GROUP
tmp=multcompare(rm1a,'ROI','by','group'); newtmp=addsig(tmp);  disp(newtmp(newtmp.Difference>0,:));
tmp=multcompare(rm1a,'group','by','ROI'); newtmp=addsig(tmp); disp(newtmp(newtmp.Difference>0,:)); % disp(newtmp(newtmp.pValue<.1 & newtmp.Difference>0,:));

%% ROI x HEMIS
tmp=multcompare(rm1a,'ROI','by','Hemis'); newtmp=addsig(tmp); disp(newtmp(newtmp.Difference>0,:));
tmp=multcompare(rm1a,'Hemis','by','ROI'); newtmp=addsig(tmp); disp(newtmp(newtmp.Difference>0,:));% disp(newtmp(newtmp.pValue<.1 & newtmp.Difference>0,:));

%% cond x GROUP
tmp=multcompare(rm1a,'Cond','by','group'); newtmp=addsig(tmp);  disp(newtmp(newtmp.Difference>0,:));
tmp=multcompare(rm1a,'group','by','Cond'); newtmp=addsig(tmp);  disp(newtmp(newtmp.Difference>0,:));% disp(newtmp(newtmp.pValue<.1 & newtmp.Difference>0,:));



%% ANOVA GROUP X COND X HEMIS for each ROI
for r=1:length(roi)
    coll=strcmp(vars,'part') | strcmp(vars,'group') | contains(vars,roinames{r});

    thbo2.(roinames{r})=thbo(:,coll);
    tmp_within=within_table(coll(3:end),:);
    tmpvars=thbo2.(roinames{r}).Properties.VariableNames;
    
    rm2.(roinames{r}) = fitrm(thbo2.(roinames{r}), [tmpvars{3} '-' tmpvars{end} '~1+group'], 'WithinDesign', tmp_within);
    
    tmptable = ranova(rm2.(roinames{r}),'WithinModel', 'Cond*Hemis');
    tbl2.(roinames{r})=addsig(tmptable);
    disp(['********************** roi ' (roinames{r})])
    disp(tbl2.(roinames{r}))
end

%% HEMIS
for r=[1];
    disp(['********************** roi ' (roinames{r})])
tmp=multcompare(rm2.(roinames{r}),'Hemis'); newtmp=addsig(tmp); disp(newtmp(newtmp.Difference>0,:));
end

%% HEMIS X GROUP
 for r=[1];
     disp(['********************** roi ' (roinames{r})])
 tmp=multcompare(rm2.(roinames{r}),'Hemis','by','group'); newtmp=addsig(tmp); disp(newtmp(newtmp.Difference>0,:));
 tmp=multcompare(rm2.(roinames{r}),'group','by','Hemis'); newtmp=addsig(tmp); disp(newtmp( newtmp.Difference>0,:));end

%% HEMIS X COND
for r=[ 3];
    disp(['********************** roi ' (roinames{r})])
tmp=multcompare(rm2.(roinames{r}),'Hemis','by','Cond'); newtmp=addsig(tmp); disp(newtmp(newtmp.Difference>0,:));
tmp=multcompare(rm2.(roinames{r}),'Cond','by','Hemis'); newtmp=addsig(tmp); disp(newtmp(newtmp.Difference>0,:));end

%% GROUP X COND
for r=[3 ];
    disp(['********************** roi ' (roinames{r})])
tmp=multcompare(rm2.(roinames{r}),'Cond','by','group'); newtmp=addsig(tmp); disp(newtmp(newtmp.Difference>0,:));
tmp=multcompare(rm2.(roinames{r}),'group','by','Cond'); newtmp=addsig(tmp); disp(newtmp(newtmp.Difference>0,:));end

%% WHEN GROUP X COND X HEMIS IS SIGNIFICANT
for r=[3 2]
    disp(['********************** roi ' (roinames{r})])
    
    for hi={'LH' 'RH'}
        disp([hi{1} ' HEMIS'])
        tmprows=[true ;true ;categorical(rm2.(roinames{r}).WithinDesign.Hemis)==hi{1}];
        tmprm=fitrm(...
            rm2.(roinames{r}).BetweenDesign(:,tmprows), ...
            [rm2.(roinames{r}).BetweenDesign(:,tmprows).Properties.VariableNames{3}  '-' ...
            rm2.(roinames{r}).BetweenDesign(:,tmprows).Properties.VariableNames{end}  '~1+group'],...
            'WithinDesign', rm2.(roinames{r}).WithinDesign(tmprows(3:end),:));
        tmpanova=ranova(tmprm,'WithinModel','Cond');tmpanova=addsig(tmpanova);disp(tmpanova)
        tmp=multcompare(tmprm,'Cond'); newtmp=addsig(tmp); disp(newtmp(newtmp.Difference>0,:));
        tmp=multcompare(tmprm,'group'); newtmp=addsig(tmp); disp(newtmp(newtmp.Difference>0,:));
        tmp=multcompare(tmprm,'Cond','by','group'); newtmp=addsig(tmp); disp(newtmp(newtmp.Difference>0,:));
        tmp=multcompare(tmprm,'group','by','Cond'); newtmp=addsig(tmp); disp(newtmp(newtmp.Difference>0,:));
    end
    
    for ci=1:3
        disp([cond{ci} ' COND'])
        tmprows=[true ;true ;categorical(rm2.(roinames{r}).WithinDesign.Cond)==cond{ci}];
        tmprm=fitrm(...
            rm2.(roinames{r}).BetweenDesign(:,tmprows), ...
            [rm2.(roinames{r}).BetweenDesign(:,tmprows).Properties.VariableNames{3}  '-' ...
            rm2.(roinames{r}).BetweenDesign(:,tmprows).Properties.VariableNames{end}  '~1+group'],...
            'WithinDesign', rm2.(roinames{r}).WithinDesign(tmprows(3:end),:));
        tmpanova=ranova(tmprm,'WithinModel','Hemis');tmpanova=addsig(tmpanova);disp(tmpanova)
        tmp=multcompare(tmprm,'Hemis'); newtmp=addsig(tmp); disp(newtmp(newtmp.Difference>0,:));
        tmp=multcompare(tmprm,'group'); newtmp=addsig(tmp); disp(newtmp(newtmp.Difference>0,:));
        tmp=multcompare(tmprm,'Hemis','by','group'); newtmp=addsig(tmp); disp(newtmp(newtmp.Difference>0,:));
        tmp=multcompare(tmprm,'group','by','Hemis'); newtmp=addsig(tmp); disp(newtmp(newtmp.Difference>0,:));
    end
    
    for gg=1:3
        disp([group{gg} ' GROUP'])
        tmprows=rm2.(roinames{r}).BetweenDesign.group==group{gg};
        tmprm=fitrm(...
            rm2.(roinames{r}).BetweenDesign(tmprows,:), ...
            [rm2.(roinames{r}).BetweenDesign(tmprows,:).Properties.VariableNames{3}  '-' ...
            rm2.(roinames{r}).BetweenDesign(tmprows,:).Properties.VariableNames{end}  '~1'],...
            'WithinDesign', rm2.(roinames{r}).WithinDesign);
        tmpanova=ranova(tmprm,'WithinModel','Cond*Hemis');tmpanova=addsig(tmpanova);disp(tmpanova)
        tmp=multcompare(tmprm,'Hemis'); newtmp=addsig(tmp); disp(newtmp(newtmp.Difference>0,:));
        tmp=multcompare(tmprm,'Cond'); newtmp=addsig(tmp); disp(newtmp(newtmp.Difference>0,:));
        tmp=multcompare(tmprm,'Hemis','by','Cond'); newtmp=addsig(tmp); disp(newtmp(newtmp.Difference>0,:));
        tmp=multcompare(tmprm,'Cond','by','Hemis'); newtmp=addsig(tmp); disp(newtmp(newtmp.Difference>0,:));
    end
end


%% save data
% save([outputdirectory datestr(now,'yyyy-mm-dd HH:MM') '-anova_grcond.mat'],'d','RES*','saveres');

%% functions

function tab_out=addsig(tab_in)
tab_in.sig=repmat('***',[height(tab_in) 1]);
tab_in.sig(tab_in.pValue>=.001,:)=repmat('** ',[sum(tab_in.pValue>=.001) 1]);
tab_in.sig(tab_in.pValue>=.01,:)=repmat('*  ',[sum(tab_in.pValue>=.01) 1]);
tab_in.sig(tab_in.pValue>=.05,:)=repmat('.  ',[sum(tab_in.pValue>=.05) 1]);
tab_in.sig(tab_in.pValue>=.1,:)=repmat('   ',[sum(tab_in.pValue>=.1) 1]);

tab_out = movevars(tab_in, 'sig', 'After', 'pValue');
end


function tripleinter(rm)

 for hi={'LH' 'RH'}
        disp([hi{1} ' HEMIS'])
        tmprows=[true;true;categorical(rm.WithinDesign.Hemis)==hi{1}];
        tmprm=fitrm(...
            rm.BetweenDesign(:,tmprows), ...
            [rm.BetweenDesign(:,tmprows).Properties.VariableNames{3}  '-' ...
            rm.BetweenDesign(:,tmprows).Properties.VariableNames{end}  '~1+group'],...
            'WithinDesign', rm.WithinDesign(tmprows(3:end),:));
        tmpanova=ranova(tmprm,'WithinModel','Cond');tmpanova=addsig(tmpanova);disp(tmpanova)
        tmp=multcompare(tmprm,'Cond'); newtmp=addsig(tmp); disp(newtmp(newtmp.Difference>0,:));
        tmp=multcompare(tmprm,'group'); newtmp=addsig(tmp); disp(newtmp( newtmp.Difference>0,:));
        tmp=multcompare(tmprm,'Cond','by','group'); newtmp=addsig(tmp); disp(newtmp( newtmp.Difference>0,:));
        tmp=multcompare(tmprm,'group','by','Cond'); newtmp=addsig(tmp); disp(newtmp(newtmp.Difference>0,:));
    end
    
    for ci={'FR' 'AL' 'HE'}
        disp([ci{1} ' COND'])
        tmprows=[true;true;categorical(rm.WithinDesign.Cond)==ci{1}];
        tmprm=fitrm(...
            rm.BetweenDesign(:,tmprows), ...
            [rm.BetweenDesign(:,tmprows).Properties.VariableNames{3}  '-' ...
            rm.BetweenDesign(:,tmprows).Properties.VariableNames{end}  '~1+group'],...
            'WithinDesign', rm.WithinDesign(tmprows(3:end),:));
        tmpanova=ranova(tmprm,'WithinModel','Hemis');tmpanova=addsig(tmpanova);disp(tmpanova)
        tmp=multcompare(tmprm,'Hemis'); newtmp=addsig(tmp); disp(newtmp( newtmp.Difference>0,:));
        tmp=multcompare(tmprm,'group'); newtmp=addsig(tmp);disp(newtmp( newtmp.Difference>0,:));
        tmp=multcompare(tmprm,'Hemis','by','group'); newtmp=addsig(tmp); disp(newtmp(newtmp.Difference>0,:));
        tmp=multcompare(tmprm,'group','by','Hemis'); newtmp=addsig(tmp); disp(newtmp(newtmp.Difference>0,:));
    end
    
    for gg={'CTL' 'XPH' 'XPA'}
        disp([gg{1} ' group'])
        tmprows=rm.BetweenDesign.group==gg{1};
        tmprm=fitrm(...
            rm.BetweenDesign(tmprows,:), ...
            [rm.BetweenDesign(tmprows,:).Properties.VariableNames{3}  '-' ...
            rm.BetweenDesign(tmprows,:).Properties.VariableNames{end}  '~1'],...
            'WithinDesign', rm.WithinDesign);
        tmpanova=ranova(tmprm,'WithinModel','Cond*Hemis');tmpanova=addsig(tmpanova);disp(tmpanova)
        tmp=multcompare(tmprm,'Hemis'); newtmp=addsig(tmp); disp(newtmp(newtmp.Difference>0,:));
        tmp=multcompare(tmprm,'Cond'); newtmp=addsig(tmp); disp(newtmp(newtmp.Difference>0,:));
        tmp=multcompare(tmprm,'Hemis','by','Cond'); newtmp=addsig(tmp); disp(newtmp(newtmp.Difference>0,:));
        tmp=multcompare(tmprm,'Cond','by','Hemis'); newtmp=addsig(tmp); disp(newtmp(newtmp.Difference>0,:));
    end
end



function plotgraph(data,group,dim2,dim3,dim4,figdirectory,vROI)
cat1=categories(group);
cat2=dim2{2};
cat3=dim3{2};
cat4=dim4{2};
xlab=categorical(cat1);xlab=reordercats(xlab,cat1);

legcat={};x=1;
for c=cat4
    for h=cat3
        legcat{x}=[c{1} '-' h{1}];
        x=x+1;
    end
end

clear y*
for g=1:length(cat1)
    idG=group==cat1{g};
    ymean(g,:,:,:)=squeeze(mean(data(idG,:,:,:),'omitnan'));
    yci(g,:,:,:)=(squeeze(std(data(idG,:,:,:),0,'omitnan'))/sqrt(sum(idG))); %SEM
end
ymean=reshape(ymean,size(ymean,1),size(ymean,2),[]);
yci=reshape(yci,size(ymean,1),size(ymean,2),[]);
clear idG g

ccolor={[.16 .18 .16] [.6 .63 .6] [0.8    0.05    0.2] [0.95    0.6    0.6] [0.2    0.3    0.75] [ 0.46    0.62    0.95]}; 

if length(cat2)==1
    figg=figure('Units','normalized','Position',[0.592708333333333,0.10462962962963,0.45,0.4],'Color','white');
    figg=tiledlayout(1,1,'TileSpacing','compact','Padding','tight');
elseif length(cat2)==2
    figg=figure('Units','normalized','Position',[0.592708333333333,0.10462962962963,0.45,0.773148148148148],'Color','white');
    figg=tiledlayout(2,1,'TileSpacing','compact','Padding','tight');
elseif length(cat2)==3
    figg=figure('Units','normalized','Position',[0.592708333333333,0.037962962962963,0.45,0.883333333333333],'Color','white');
    figg=tiledlayout(3,1,'TileSpacing','compact','Padding','tight');
elseif length(cat2)==4
    figg=figure('Units','normalized','Position',[0.3671875,0.10462962962963,0.66,0.72037037037037],'Color','white');
    figg=tiledlayout(2,2,'TileSpacing','compact','Padding','tight');
elseif length(cat2)==5 || length(cat2)==6
    figg=figure('Units','normalized','Position',[0.416145833333333,0.037962962962963,0.66,0.883333333333333],'Color','white');
    figg=tiledlayout(3,2,'TileSpacing','compact','Padding','tight');
end
yyaxis=[];
for r=1:length(cat2)
    nexttile;
    y=squeeze(ymean(:,r,:));
    err=squeeze(yci(:,r,:));
    b.(['roi' num2str(r)])=bar(xlab,y,'BarWidth',.85);
    for a=1:length(b.(['roi' num2str(r)]))
        b.(['roi' num2str(r)])(a).FaceColor=ccolor{a};
    end
    ax = gca;
    ax.FontSize=12;
    
    ngroups = size(ymean, 1);
    nbars = size(ymean, 3);
    % Calculating the width for each bar group
    groupwidth = min(0.85, nbars/(nbars + 1.5));
    for i = 1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        hold on; errorbar(x, y(:,i), err(:,i), '.k','LineWidth',1);
    end
    clear i x y err ax a
    yyaxis=[yyaxis ylim];
    ylabel('Mean \DeltaHbO (\muM)','FontSize',14)
    title([cat2{r} ' ROI'])
    
end
for r=1:length(cat2)
    nexttile(r);
    ylim([min(yyaxis) max(yyaxis)])
end
xlabel('Groups','FontSize',14)
lgd=legend(legcat,'FontSize',13);
legend('boxoff')
lgd.Layout.Tile = 'east';
title(lgd,'Conditions')
if ~exist(figdirectory,'dir')
    mkdir(figdirectory);
end

tmpname=[figdirectory vROI 'ROIinteraction_v1.fig'];
tmpn=2;
while exist(tmpname,'file')
    tmpname=[tmpname(1:end-5) num2str(tmpn) '.fig'];
    tmpn=tmpn+1;
end
saveas(figg,tmpname)
saveas(figg,[tmpname(1:end-4) '.png'])
end

function plotgraphmain(data,group,dim2,dim3,dim4,figdirectory,vROI)
cat1=categories(group);
cat2=dim2{2};
cat3=dim3{2};
cat4=dim4{2};
xlab=categorical(cat1);xlab=reordercats(xlab,cat1);

legcat={};x=1;
for c=cat4
    for h=cat3
        legcat{x}=[c{1} '-' h{1}];
        x=x+1;
    end
end

clear y*
for g=1:length(cat1)
    idG=group==cat1{g};
    ymean(g,:,:)=squeeze(mean(mean(data(idG,:,:,:),2,'omitnan'),'omitnan'));
    yci(g,:,:)=(squeeze(std(mean(data(idG,:,:,:),2,'omitnan'),0,'omitnan'))/sqrt(sum(idG))); %SEM
end
ymean=reshape(ymean,size(ymean,1),[]);
yci=reshape(yci,size(ymean,1),[]);
clear idG g

ccolor={[.16 .18 .16] [.6 .63 .6] [0.8    0.05    0.2] [0.95    0.6    0.6] [0.2    0.3    0.75] [ 0.46    0.62    0.95]}; 


    figg=figure('Units','normalized','Position',[0.592708333333333,0.10462962962963,0.45,0.4],'Color','white');
    figg=tiledlayout(1,1,'TileSpacing','compact','Padding','tight');

    nexttile;
    y=ymean;
    err=yci;
    b=bar(xlab,y,'BarWidth',.85);
    for a=1:length(b)
        b(a).FaceColor=ccolor{a};
    end
    ax = gca;
    ax.FontSize=12;
    
    ngroups = size(ymean, 1);
    nbars = size(ymean, 2);
    % Calculating the width for each bar group
    groupwidth = min(0.85, nbars/(nbars + 1.5));
    for i = 1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        hold on; errorbar(x, y(:,i), err(:,i), '.k','LineWidth',1);
    end
    clear i x y err ax a
    
    ylabel('Mean \DeltaHbO (\muM)','FontSize',14)
    title(['Average activation across all ROIs'])
    

xlabel('Groups','FontSize',14)
lgd=legend(legcat,'FontSize',13);
legend('boxoff')
lgd.Layout.Tile = 'east';
title(lgd,'Conditions')
if ~exist(figdirectory,'dir')
    mkdir(figdirectory);
end

tmpname=[figdirectory vROI 'AVGROI_interactionGCH_v1.fig'];
tmpn=2;
while exist(tmpname,'file')
    tmpname=[tmpname(1:end-5) num2str(tmpn) '.fig'];
    tmpn=tmpn+1;
end
saveas(figg,tmpname)
saveas(figg,[tmpname(1:end-4) '.png'])
end
