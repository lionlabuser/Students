function [tblres, tblpval, tblpvalsig, tblressig, tblfacpvalsig, tblfacressig] = Tablegraph_CO(res, pval, labelfactors, p, dat, lgnd, label, MATmeandiff, meanG1, meanG2, DATA, fileorderconnectogram, savepath, graphmode, varargin)
    
    nbfactors = length(labelfactors);
    R = find(contains(lgnd(1,:),dat));

    tblres = array2table(res,'VariableName',label);
    tblpval = array2table(pval,'VariableName',label,'RowNames',labelfactors');
    
    sig = tblpval.Variables <= p;
    tblpvalsig = tblpval(:,any(sig));
    tblressig = tblres(:,tblpvalsig.Properties.VariableNames);
    
    save([savepath date '_' dat 'Results.mat'],'tblres','tblpval','tblpvalsig','tblressig');
    %writetable(tblsigch,[savepath date '_sigresultsch.xls']); trop gros
    writetable(tblpvalsig,[savepath date '_' dat 'PvalSigresults.xls'],'WriteRowNames',true);
    %writetable(tblgrsigch,[savepath date '_grsigresultsch.xls']);
    
    if find(sig)
        for f = 1:numel(labelfactors)
            fac = labelfactors{f,1};
            if contains(fac,'*')
                fac = erase(fac,'*');
            elseif contains(fac,':')
                fac = erase(fac,':');
            end
            sig = tblpval{f,:} <= p;
            if find(sig)
                tblfacpvalsig.(fac) = tblpval(f,sig);
                tblfacressig.(fac) = tblres(:,sig);
                disp(tblfacpvalsig.(fac))

                if graphmode == 1
                    %%graphique%%%
                    figure
                    A = meanG1(:,sig);
                    B = meanG2(:,sig);
                    C = [A' B'];
                    X = categorical(tblfacpvalsig.(fac).Properties.VariableNames);
                    p1 = bar(X,C);
                    p1(1).FaceColor = 'r';
                    p1(2).FaceColor = 'b';
                    ylabel('Pearson correlation');
                    str = sprintf('Significant %s effect',fac);
                    title(str)
                    savefig([savepath date '_' dat 'Sig' fac]);
                    exportgraphics(gcf,[savepath date '_' dat 'Sig' fac '.png'])

                    % figure
                    % p1 = bar(meanG1ch(:,grsig),'r');
                    % hold on
                    % p2 = bar(meanG2ch(:,grsig),'b');
                    % xlabel(tblgrpsigch.Properties.VariableNames);
                    % ylabel('Pearson correlation');
                    % hold off

                    clear r sig str A B C X p1

                    %%%% Connectogramme%%%%%%%%%%%%
                    %    %MATmeanG2G1 = MATmeanG2-MATmeanG1; %Calculer matrice G2-G1
                    %    MATmeanG1G2 = MATmeanG1 - MATmeanG2;

                    MATpval = zeros(length(MATmeandiff));
                    if contains(dat,'Ch')
                       for c = 1:length(MATmeandiff)
                            for cc =(c + 1):length(MATmeandiff)
                                x = [num2str(c) '-' num2str(cc)];
                                tf = strcmp(x, label);
                                idx = find(tf);
                                MATpval(c,cc)= pval(f,idx); %Mettre les  p sous forme de matrice
                            end
                        end 
                    elseif contains(dat,'roi')
                        for c = 1:length(MATmeandiff)
                            for cc =(c + 1):length(MATmeandiff)
                                x = [lgnd{2,R} num2str(c) '-' lgnd{2,R} num2str(cc)];
                                tf = strcmp(x, label);
                                idx = find(tf);
                                MATpval(c,cc)= pval(f,idx); %Mettre les  p sous forme de matrice
                            end
                        end
                    else
                        error('Problem with data type')
                    end

                    MATpval(length(MATpval),length(MATpval)) = 0;
                    MATpval = MATpval + triu(MATpval,1)';
                    tf = MATpval <=p & MATpval > 0;
                    %         idx = find(tf);
                    %MATsigG2G1ch = MATmeanG2G1;
                    MATsig = MATmeandiff;
                    %MATsigG2G1ch(~tf) = 0;
                    MATsig(~tf) = 0;

                    clear c cc x tf idx

                    id = 1;
                    %         MAT = MATsigch; %DATA{id}.MAT %loader matrice
                    %         List = strvcat(DATA{id}.ZoneList); %liste des paires SD
                    %         ML = DATA{id}.zone.ml; %Loader S/D/ROI/Gr
                    %         plotLst = DATA{id}.zone.plotLst;
                    %         label =  DATA{id}.zone.label;
                    %         fileorderconnectogram  = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Connectogram_Mixte.txt';
                    %         plotconnectogram(fileorderconnectogram,MAT,List,label,plotLst,ML)
                    %         savefig([savepath date '_ConnectSigCh']);

                    MATneg = MATsig < 0;
                    MAT = MATsig.*(MATneg); %DATA{id}.MAT %loader matrice
                    if any(MAT,'all') && contains(dat,'Ch')
                        List = strvcat(DATA{id}.ZoneList); %liste des paires SD
                        ML = DATA{id}.zone.ml; %Loader S/D/ROI/Gr
                        plotLst = DATA{id}.zone.plotLst;
                        labelzone =  DATA{id}.zone.label;
                        plotconnectogram(fileorderconnectogram{1,R},MAT,List,labelzone,plotLst,ML)
                        savefig([savepath date '_' dat 'NegConnectSigG1G2' fac]);
                    elseif any(MAT,'all') && contains(dat,'roi')
                        plotLst = DATA(2,:);
                        labelzone =  DATA(1,:);
                        plotconnectogramroi(fileorderconnectogram{1,R},MAT,labelzone,plotLst)
                        savefig([savepath date '_' dat 'NegConnectSigG1G2' fac]);
                    else
                    end

                    MATpos = MATsig > 0;
                    MAT = MATsig.*(MATpos);%DATA{id}.MAT %loader matrice
                    if any(MAT,'all') && contains(dat,'Ch')
                        List = strvcat(DATA{id}.ZoneList); %liste des paires SD
                        ML = DATA{id}.zone.ml; %Loader S/D/ROI/Gr
                        plotLst = DATA{id}.zone.plotLst;
                        labelzone =  DATA{id}.zone.label;
                        plotconnectogram(fileorderconnectogram{1,R},MAT,List,labelzone,plotLst,ML)
                        savefig([savepath date '_' dat 'PosConnectSigG1G2' fac]);
                    elseif any(MAT,'all') && contains(dat,'roi')
                        plotLst = DATA(2,:);
                        labelzone =  DATA(1,:);
                        plotconnectogramroi(fileorderconnectogram{1,R},MAT,labelzone,plotLst)
                        savefig([savepath date '_' dat 'PosConnectSigG1G2' fac]);
                    else
                    end
                end
                save([savepath date '_' dat 'Results.mat'],'tblfacpvalsig','tblfacressig','-append');
            else
            end
        end
    end
    
    if nargin == 19
        idx = find(contains(labelfactors,':')); %trouver le num du facteur
        idsiginter = find(pval(idx,:)<= p); %trouver les interactions sig
        
        tblcoef = array2table(varargin{1},'VariableName',label);
        tblcoefsiginterG1 = array2table(varargin{2},'VariableNames',label(idsiginter));
        tblcoefsiginterG2 = array2table(varargin{3},'VariableNames',label(idsiginter));
        tblpvalsiginterG1 = array2table(varargin{4},'VariableNames',label(idsiginter));
        tblpvalsiginterG2 = array2table(varargin{5},'VariableNames',label(idsiginter));
        save([savepath date '_' dat 'Results.mat'],'tblcoef', 'tblcoefsiginterG1', 'tblcoefsiginterG2', 'tblpvalsiginterG1', 'tblpvalsiginterG2','-append');
    end
        
end