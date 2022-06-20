function [tblres, tblpval, tblpvalsig, tblressig, tblfacpvalsig, tblfacressig] = ...
    Tablegraph_CO(res, pval, labelfactors, p, dat, lgnd, label, MATmeandiff, meanG1,...
    meanG2, DATA, fileorderconnectogram, savepath, graphmode, varargin)
    
    nbfactors = length(labelfactors);
    R = find(contains(lgnd(1,:),dat));

    tblres = array2table(res,'VariableName',label);
    tblpval = array2table(pval,'VariableName',label,'RowNames',labelfactors');
    
    sig = tblpval.Variables <= p;
    tblpvalsig = tblpval(:,any(sig));
    tblressig = tblres(:,tblpvalsig.Properties.VariableNames);
    
    save([savepath date '_' dat 'Results.mat'],'tblres','tblpval','tblpvalsig','tblressig');
    %writetable(tblsigch,[savepath date '_sigresultsch.xls']); trop gros
try
    writetable(tblpvalsig,[savepath date '_' dat 'PvalSigresults.xls'],'WriteRowNames',true);
catch
end
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
                    fig = figure;
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
                    fig.WindowState = 'maximized';
                    savefig([savepath date '_' dat 'Sig' fac]);
                    exportgraphics(gcf,[savepath dat 'Sig' fac '.png'])
                    close

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
                        x = 1;
                       for c = 1:length(MATmeandiff)
                            for cc = (c + 1):length(MATmeandiff)
%                                 x = [num2str(c) '-' num2str(cc)];
%                                 tf = strcmp(x, label);
%                                 idx = find(tf);
                                MATpval(c,cc)= pval(f,x); %Mettre les  p sous forme de matrice
                                x = x + 1;
                            end
                        end 
                    elseif contains(dat,'roi')
                        x = 1;
                        for c = 1:length(MATmeandiff)
                            for cc =(c + 1):length(MATmeandiff)
%                                 x = [lgnd{2,R} num2str(c) '-' lgnd{2,R} num2str(cc)];
%                                 tf = strcmp(x, label);
%                                 idx = find(tf);
                                MATpval(c,cc) = pval(f,x); %Mettre les  p sous forme de matrice
                                x = x + 1;
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
                    MATneg = MATsig < 0;
                    MAT = MATsig.*(MATneg); %DATA{id}.MAT %loader matrice
                    if any(MAT,'all') && contains(dat,'Ch')
                        List = strvcat(DATA{id}.ZoneList); %liste des paires SD
                        ML = DATA{id}.zone.ml; %Loader S/D/ROI/Gr
                        plotLst = DATA{id}.zone.plotLst;
                        labelzone =  DATA{id}.zone.label;
                        labelzone = strrep(labelzone, '_', ' ');
                        plotconnectogram(fileorderconnectogram{1,R},MAT(varargin{1,1},varargin{1,1}),List,labelzone,plotLst,ML)
                        fig = gcf;
                        %pbaspect([1 1 1])
                        fig.WindowState = 'maximized';
                        colorbar('Position',[0.663 0.334 0.02 0.4],'FontSize',12)
                        savefig([savepath date '_' dat 'NegConnectSigG1G2' fac]);
                        exportgraphics(gcf,[savepath dat 'NegConnectSigG1G2' fac '.png'])
                        close
                    elseif any(MAT,'all') && contains(dat,'roi')
                        plotLst = DATA(2,:);
                        labelzone =  DATA(1,:);
                        plotconnectogramroi(fileorderconnectogram{1,R},MAT,labelzone,plotLst)
                        fig = gcf;
                        %pbaspect([1 1 1])
                        fig.WindowState = 'maximized';
                        if contains(dat,'roiR')
                            colorbar('Position',[0.75 0.334 0.02 0.4],'FontSize',12)
                        else
                            colorbar('Position',[0.663 0.334 0.02 0.4],'FontSize',12)
                        end
                        savefig([savepath date '_' dat 'NegConnectSigG1G2' fac]);
                        exportgraphics(gcf,[savepath dat 'NegConnectSigG1G2' fac '.png'])
                        close
                    else
                    end

                    MATpos = MATsig > 0;
                    MAT = MATsig.*(MATpos);%DATA{id}.MAT %loader matrice
                    if any(MAT,'all') && contains(dat,'Ch')
                        List = strvcat(DATA{id}.ZoneList); %liste des paires SD
                        ML = DATA{id}.zone.ml; %Loader S/D/ROI/Gr
                        plotLst = DATA{id}.zone.plotLst;
                        labelzone =  DATA{id}.zone.label;
                        labelzone = strrep(labelzone, '_', ' ');
                        plotconnectogram(fileorderconnectogram{1,R},MAT(varargin{1,1},varargin{1,1}),List,labelzone,plotLst,ML)
                        fig = gcf;
                        %pbaspect([1 1 1])
                        fig.WindowState = 'maximized';
                        colorbar('Position',[0.663 0.334 0.02 0.4],'FontSize',12)
                        savefig([savepath date '_' dat 'PosConnectSigG1G2' fac]);
                        exportgraphics(gcf,[savepath dat 'PosConnectSigG1G2' fac '.png'])
                        close
                    elseif any(MAT,'all') && contains(dat,'roi')
                        plotLst = DATA(2,:);
                        labelzone =  DATA(1,:);
                        plotconnectogramroi(fileorderconnectogram{1,R},MAT,labelzone,plotLst)
                        fig = gcf;
                        %pbaspect([1 1 1])
                        fig.WindowState = 'maximized';
                        if contains(dat,'roiR')
                            colorbar('Position',[0.75 0.334 0.02 0.4],'FontSize',12)
                        else
                            colorbar('Position',[0.663 0.334 0.02 0.4],'FontSize',12)
                        end
                        savefig([savepath date '_' dat 'PosConnectSigG1G2' fac]);
                        exportgraphics(gcf,[savepath dat 'PosConnectSigG1G2' fac '.png'])
                        close
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