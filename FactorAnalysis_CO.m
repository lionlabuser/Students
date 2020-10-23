%%%%%%%%%%%%%%FACTOR ANALYSIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
disp('Computing Factor Analysis')

data = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Stats\CORR\';
load ([data 'workspace.mat'])

savepath='C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Stats\CORR\Factor\';
channelmode = 1;
tablegraphmode = 1;

%%CHANNELS%%%%%%%%%%%%%%%%%%%%%
if channelmode ==1

    X = table2array(tblch(:,6:end));
%     X = X(~isnan(X));
    
   [lambda,psi,T] = factoran(X,8);
end

%%Valeurs manquantes pas acceptées, utiliser des données imputées%%

%     %%tableaux des résultats%%%%%%%%%%%%%%%%
%     if tablegraphmode == 1
%         tblpch = array2table(pch);
%         tblpch.Properties.VariableNames = labelch;
%         tblpch.Properties.RowNames = {'Group'};
% 
%         sig = tblpch.Variables <= 0.07;
%         tblpsigch = tblpch(:,sig);
% 
%         disp(tblpsigch)
% 
%         %writetable(tblsigch,[savepath date '_sigresultsch.xls']); trop gros
%         writetable(tblpsigch,[savepath date '_psigresultsch.xls'],'WriteRowNames',true);
%         save([savepath date '_resultsch.mat'],'tblpch','tblpsigch');
% 
%         %%%% graphique%%%%%%%%%%%%%%%%%%%%
%         figure
%         A = meanG1ch(:,sig);
%         B = meanG2ch(:,sig);
%         C = [A' B'];
%         X = categorical(tblpsigch.Properties.VariableNames);
%         p1 = bar(X,C);
%         p1(1).FaceColor = 'r';
%         p1(2).FaceColor = 'b';
%         ylabel('Pearson correlation');
%         savefig([savepath date '_SigCh']);
% 
%         % figure
%         % p1 = bar(meanG1ch(:,grsig),'r');
%         % hold on
%         % p2 = bar(meanG2ch(:,grsig),'b');
%         % xlabel(tblgrpsigch.Properties.VariableNames);
%         % ylabel('Pearson correlation');
%         % hold off
%     end
% end
% 
% clear r x res sig A B C X p1
% 