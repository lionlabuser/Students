
function openplot(path)
% GUI Plot from Julie's toolbox LIONirs
% directory = directory of NIRS.mat
%             eg. 'C:\data\Data_NIRS\ELAN\ANALYSED\Martine_0m\BB058'
 job.NIRSmat = {[path '\NIRS.mat']};
 job.DelPreviousData = 0;
 plot_sessions_GUI(job.NIRSmat,1,job);
end