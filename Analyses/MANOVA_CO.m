%%%%%%%%%%%%%%MANOVA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
disp('Computing MANOVA')

data = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\CORR\';
load ([data 'workspace.mat'])

savepath='C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\CORR\MANOVA\';


%[dch, pch, statsch] = manova1(table2array(tblch(:,6:end)),gr); %singular nobs < npart
%[dMroi, pMroi, statsMroi] = manova1(table2array(tblmeanMroi(:,6:end)),gr); %singular nobs < npart
%[dRroi, pRroi, statsRroi] = manova1(table2array(tblmeanRroi(:,6:end)),gr); %singular nobs < npart
%[dFroi, pFroi, statsFroi] = manova1(table2array(tblmeanFroi(:,6:end)),gr); %singular nobs < npart

%Problème de singularité, trop de données

[dAroi, pAroi, statsAroi] = manova1(table2array(tblmeanAroi(:,6:end)),gr);

if pAroi <= 0.07;
   fprintf('p= %.3f, MANOVA for pAroi is significant p<=0.07\n',pAroi)
else
   fprintf('p= %.3f, MANOVA for pAroi is not significant p<=0.07\n',pAroi)   
end

save([savepath date '_resultsAroi.mat'],'dAroi','pAroi','statsAroi');


X = ['Results saved in ', savepath];
disp(X)
clear X

toc