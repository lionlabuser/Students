%%%%%%%%%%%%%%%%%%%%%%%%%% Graph Theory%%%%%%%%%%%%%%%%%%%%
datapath = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Stats\CORR0,01_0,08\';
%load ([datapath 'workspace.mat'])
load ([datapath 'workspacemat.mat'])

savepath='C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Stats\CORR0,01_0,08\GraphTheory\';
if ~isfolder(savepath)
    mkdir(savepath)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%BCT 

%BCT threshold 
%Weigthed = garder les valeurs de corr, Binary = remplacer tout par 0 ou 1.%

itr = 0.1:0.01:1.0;% Seuil variable BCT soustraction 
fixetr = 0.2; %pour les figures

%Global efficiency, une valeur par threshold par participant
disp('Global Efficiency, running')
GE = zeros(size(MATall,1), numel(itr));
%idtr = 1;
for isubject = 1:size(MATall,3) % pour chaque sujet
    for idtr = 1:numel(itr) % pour chaque threshold
        itrv = itr(idtr); %extraire le threshold actuel
        A = squeeze(MATall(:,:,isubject)); %prendre la matrice de chaque participant%
        A = threshold_absolute(A, itrv); %Rendre la matrice binaire en fonction de chaque threshold%
        A(isnan(A)) = 0; %% Mettre les NAN à 0
        GE(isubject, idtr)=efficiency_bin(A,0); %calculer le GE pour chaque matrice binarisée de chaque participant%              
    end   
end
GEG1 = nanmean(GE(idG1,:),1); %moyenne du G1 pour chaque threshold
GEG2 = nanmean(GE(idG2,:),1); % moyenne du G2 pour chaque threshold


disp('Global Efficiency, done')
save([savepath date '_GE.mat'],'GE','GEG1','GEG2');


figure
subplot(2,4,1);
hold on;
plot(itr,mean(GE(idG1,:),1),'r','displayname','mal')
plot(itr,mean(GE(idG2,:),1),'b','displayname','ctl')
xlabel('Threshold','fontsize',12)
ylabel('Global Efficiency','fontsize',12)

id = sum(itr< fixetr); %trouver la position du threshold sélectionné
d1 = GE(idG1,id);
d2 = GE(idG2,id);

subplot(2,4,5)
hold on
bar(1,mean(d1),'r')
bar(2,mean(d2),'b')
plot(1,d1,'x')
plot(2,d2,'x')
title(['Threshold <', num2str(fixetr)])
xlabel('Groups')
ylabel('Global efficiency')
xlim([0 , 3])

%Local Efficiency Input Distance matrix, une valeur par canal, par threshold, par participant%
disp('Local Efficiency, running')
idtr = 1;
LE = zeros(size(MATall,1),size(MATall,3), numel(itr));
for isubject = 1:size(MATall,3)
    for idtr = 1:numel(itr)
        itrv = itr(idtr);
        A = squeeze(MATall(:,:,isubject)); %prendre la matrice de chaque participant%
        A = threshold_absolute(A, itrv); %Rendre la matrice binaire en fonction de chaque threshold%
        A(isnan(A)) = 0;  %% Mettre les NAN à 0
        LE(:,isubject, idtr)= efficiency_bin(A,1); %calculer le LE pour chaque canal, de chaque matrice binarisée de chaque participant%
    end   
end
%faire la moyenne du LE de tous les canaux, pour chaque matrice de chaque participant%
LEmean=squeeze(mean(LE,1));
LEmeanG1 = nanmean(LEmean(idG1,:),1); %extraire la moyenne de G1
LEmeanG2 = nanmean(LEmean(idG2,:),1); % extraire la moyenne de G2

disp('Local Efficiency, done')
save([savepath date '_LE.mat'],'LE','LEmean','LEmeanG1','LEmeanG2');

subplot(2,4,2)
hold on
plot(itr,mean(LEmean(idG1,:),1),'r','displayname','mal')
plot(itr,mean(LEmean(idG2,:),1),'b','displayname','ctl')
xlabel('Threshold','fontsize',12)
ylabel('Local Efficiency','fontsize',12)

id = sum(itr< fixetr);
clear d1 d2
d1 = LEmean(idG1,id);
d2 = LEmean(idG2,id);

subplot(2,4,6)
hold on
bar(1,mean(d1),'r')
bar(2,mean(d2),'b')
plot(1,d1,'x')
plot(2,d2,'x')
title(['Threshold <', num2str(fixetr)])
xlabel('Groups')
ylabel('Local efficiency')
xlim([0 , 3])


% Characteristic pathlenght charpath(D)(BU, BD, WU, WD networks), une valeur par threshold par participant.
disp('Characteristic path length, running')
idtr = 1;
LL = zeros(size(MATall,3), numel(itr));
for isubject = 1:size(MATall,3)
    for idtr = 1:numel(itr)
        itrv = itr(idtr);
        A = squeeze(MATall(:,:,isubject)); %prendre la matrice de chaque participant%
        A = threshold_absolute(A, itrv); %Rendre la matrice binaire en fonction de chaque threshold%
        A(isnan(A)) = 0; %Mettre les NAN à 0
        D = distance_bin(A); %calculer la matrice de distance%
        [lambda,efficiency,ecc,radius,diameter] = charpath(D,1,0); 
        LL(isubject, idtr) = lambda; %calculer le LL pour chaque matrice de distance binarisée de chaque participant%              
    end   
end
LLG1 = nanmean(LL(idG1,:),1);
LLG2 = nanmean(LL(idG2,:),1);

disp('Characteristic path length, done')
save([savepath date '_LL.mat'],'LL','LLG1','LLG2');

subplot(2,4,3)
hold on
plot(itr,mean(LL(idG1,:),1),'r','displayname','mal')
plot(itr,mean(LL(idG2,:),1),'b','displayname','ctl')
xlabel('Threshold','fontsize',12)
ylabel('Characteristic path lenght','fontsize',12)

id = sum(itr< fixetr);
clear d1 d2
d1 = LL(idG1,id);
d2 = LL(idG2,id);

subplot(2,4,7);
hold on
bar(1,mean(d1),'r')
bar(2,mean(d2),'b')
plot(1,d1,'x')
plot(2,d2,'x')
xlabel('Groups')
ylabel('Characteristic path lenght')
title(['Threshold <', num2str(fixetr)])


%Clustering Coefficient clustering_coef_bu(A)
disp('Clustering coefficient, running')
idtr = 1;
CC = zeros(size(MATall,1),size(MATall,3), numel(itr));
for isubject = 1:size(MATall,3)
    for idtr = 1:numel(itr)
        itrv = itr(idtr);
        A = squeeze(MATall(:,:,isubject)); %prendre la matrice de chaque participant%
        A = threshold_absolute(A, itrv); %Rendre la matrice binaire en fonction de chaque threshold%
        A(isnan(A)) = 0;
        CC(:,isubject, idtr) = clustering_coef_bu(A); %calculer le CC pour chaque canal, de chaque matrice de chaque participant%
    end
end
CCmean = squeeze(mean(CC,1)); %faire la moyenne du CC de tous les canaux, pour chaque matrice binarisée de chaque participant%
CCmeanG1 = nanmean(CCmean(idG1,:),1);
CCmeanG2 = nanmean(CCmean(idG2,:),1);


disp('Clustering coefficient, done')
save([savepath date '_CC.mat'],'CC','CCmean','CCmeanG1','CCmeanG2');

subplot(2,4,4)
hold on
plot(itr,mean(CCmean(idG1,:),1),'r','displayname','mal')
plot(itr,mean(CCmean(idG2,:),1),'b','displayname','ctl')
xlabel('Threshold','fontsize',12)
ylabel('Clustering Coefficient','fontsize',12)

id = sum(itr< fixetr);
clear d1 d2
d1 = CCmean(idG1,id);
d2 = CCmean(idG2,id);

subplot(2,4,8)
hold on
bar(1,nanmean(d1),'r')
bar(2,nanmean(d2),'b')
plot(1,d1,'x')
plot(2,d2,'x')
xlabel('Groups')
ylabel('Clustering Coefficient')
title(['Threshold <', num2str(fixetr)])

savefig([savepath date 'Metrics.fig'])
f=gcf;
exportgraphics(f,[savepath date 'Metrics.png'])


%small worldness index%
disp('Small Worldness index, running')
SW = CCmean./LL;
SWG1 = nanmean(SW(idG1,:),1);
SWG2 = nanmean(SW(idG2,:),1);

disp('Small Worldness index, done')
save([savepath date '_SW.mat'],'SW','SWG1','SWG2');

figure
subplot(1,2,1)
hold on
plot(itr,squeeze(nanmean(CCmean(idG1,:),1)./mean(LL(idG1,:),1)),'r','displayname','mal')
plot(itr,squeeze(nanmean(CCmean(idG2,:),1)./mean(LL(idG2,:),1)),'b','displayname','ctl')
xlabel('Threshold','fontsize',16)
ylabel('Clustering Coefficient/Characteristic Path length','fontsize',16)
title('Small World Index')

subplot(1,2,2)
hold on
plot(itr,squeeze(nanmean(CCmean(idG1,:),1)),'r','displayname','CCmal')
plot(itr,squeeze(mean(LL(idG1,:),1)),'r--','displayname','LLmal')
plot(itr,squeeze(nanmean(CCmean(idG2,:),1)),'b','displayname','CCctl')
plot(itr,squeeze(mean(LL(idG2,:),1)),'b--','displayname','LLctl')
xlabel('Threshold','fontsize',16)
ylabel('Metric value','fontsize',16)
title('Small World Index')
legend

savefig([savepath date,'Comparison.fig'])
f=gcf;
exportgraphics(f,[savepath date 'Comparison.png'])

%%%%% ANOVA %%%%%%%%%%%

% pch(:,x),res] = anovan(table2array(tblch(:,r+5)),{gr,ses},'Continuous',2,'varnames',{'Group','SES'});
%         
% n_sig = sum(pch(1,:) <= p);
% fprintf('%d tests are significant p<=%.2f without correction for channels\n',n_sig, p)
    