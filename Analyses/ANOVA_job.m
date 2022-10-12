function [tblpval, tblres, varargout] = ANOVA_job(type, data,...
    factors, labelfactors, p, labeldata, labeldim, fdrmode, graphmode, savepath, savemode)
try
    
%savepath = fullfile(savepath,['ANOVA' type],filesep);

    if graphmode
        savepathgraph = fullfile(savepath,'Graphs',filesep);
        if ~isfolder(savepathgraph)
            mkdir(savepathgraph)
        end
    end
    
    nfactors = numel(factors);
    switch type
        case 'inter'
    x = 1;
    for n = 1:nfactors
        for nn = (n+1):nfactors
            labelinter{n} = [labelfactors{n} labelfactors{nn}];
            lgndinter{x,:} = [n,nn];
            ninter = (x-1) + 1;
        end
        comb(n) = nchoosek(nfactors,n);
    end
    ninter = sum(comb) - nfactors;
    labeleffects = [labelfactors, labelinter]';
    lgndinter = [num2cell((1:nfactors)');lgndinter];
    clear x n nn
        case 'main'
            labeleffects = labelfactors;
    end

    sz = size(data);

    if numel(sz) == 2 && sz(2) == 1
        switch type
            case 'inter'
                [pval, resdata] = anovan(real(data),factors,'Continuous',2,'model','full','varnames',labelfactors,'display','off');
            case 'main'
                [pval, resdata] = anovan(real(data),factors,'Continuous',2,'varnames',labelfactors,'display','off');
        end

        labelrows = {resdata{2:end-2,1,1}}';
        for f = 1:numel(labelrows) %for each effect tested
            n_sig(f,:) = sum(pval(f,:) <= p);
            fprintf('\t \t %d %s effects significant p<=%.2f without correction for %s\n',n_sig(f,:), labelrows{f,1}, p, labeldata{2})
        end

        tblres = array2table(resdata(2:end,2:end),'RowNames',resdata(2:end,1),'VariableNames',resdata(1,2:end));
        tblpval = array2table(pval,'VariableName',{'Pval'},'RowNames',labelrows');

        if savemode
            save([savepath date '_' labeldata{2} 'ResultsANOVA' type '.mat'],'tblres','tblpval');
        end

        if graphmode
            for f = 1:numel(labeleffects) %for each effect tested (simple effects and two sided interaction)
                fac = labeleffects{f};
                if n_sig(f,:) > 0 && f<=nfactors %principal effect
                    if min(mod(factors{f},1) == 0) == 0 %continuous
                        figure
                        idx = ~or(isnan(factors{f}), isnan(data));
                        tfit = fit(factors{f}(idx), double(data(idx)),'poly1');
                        plot(tfit, factors{f}, data, 'o')
                        R = corrcoef(factors{f}, data);
                        str = sprintf('y = %.3f x + %.3f\n R2 = %.3f', tfit.p1, tfit.p2, R(2)^2);
                        dim = [.15 .6 .3 .3];
                        annotation('textbox',dim,'String',str,'FitBoxToText','on');
                        xlabel(labelfactors{f},'fontsize',14)
                        ylabel(strrep(labeldata{1},'_',' '),'fontsize',14)
                        title(['Mean ' strrep(labeldata{1},'_',' ') ' of ' strrep(labeldata{2},'_',' ') ' by ' labeleffects{f}],'fontsize',14)
                        fig = gcf;
                        savefig([savepathgraph labeleffects{f} 'effect_' labeldata{2} '.fig'])
                        exportgraphics(fig,[savepathgraph labeleffects{f} 'effect_' labeldata{2} '.png'])
                        close
                        clear R str dim tfit idx
                    elseif  min(mod(factors{f},1) == 0) == 1 %categorical
                        figure
                        meang1 = mean(data(factors{f} == 1),'omitnan');
                        meang2 = mean(data(factors{f} == 2),'omitnan');
                        hold on
                        bar(1,meang1,'r'); %b1
                        bar(2,meang2,'b'); %b2
                        xlabel(labelfactors{f},'fontsize',14)
                        ylabel(strrep(labeldata{1},'_',' '),'fontsize',14)
                        xlim([0 , 3])
                        title(['Mean ' strrep(labeldata{1},'_',' ') ' of ' strrep(labeldata{2},'_',' ') ' by ' labelfactors{f}],'fontsize',14)
                        hold off
                        box on
                        legend Mal Ctl
                        fig = gcf;
                        savefig([savepathgraph labeleffects{f} 'effect_' labeldata{2} '.fig'])
                        exportgraphics(fig,[savepathgraph labeleffects{f} 'effect_' labeldata{2} '.png'])
                        close
                    end
                elseif n_sig(f,:) > 0 && f>nfactors %interaction effect %only for two sided interactions
                    for a = 1:nfactors
                        for b = (a+1):nfactors
                            %factor1 = factors{a}; %Group factor1
                            %factor2 = factors{b}; %SES factor2
                            %min(mod(factors{f},1) == 0) == 1 %cat
                            %min(mod(factors{f},1) == 0) == 0 %continuous
                            if min(mod(factors{a},1) == 0) == 0 && min(mod(factors{b},1) == 0) == 1 %fact1cont %fact2cat
                                figure %SES*Sex factor 2 and factor 3
                                sesg1 = factors{a}(factors{b} == 1);
                                datag1 = data(factors{b} == 1);
                                idx = ~or(isnan(sesg1),isnan(datag1));
                                tfitg1 = fit(sesg1(idx), double(datag1(idx)),'poly1');
                                sesg2 = factors{a}(factors{b} == 2);
                                datag2 = data(factors{b} == 2);
                                idx = ~or(isnan(sesg2),isnan(datag2));
                                tfitg2 = fit(sesg2(idx), double(datag2(idx)),'poly1');
                                hold on
                                plot(tfitg1, '-r',sesg1, datag1, 'or')
                                plot(tfitg2, '-b',sesg2, datag2, 'ob')
                                Rg1 = corrcoef(sesg1, datag1);
                                Rg2 = corrcoef(sesg2, datag2);
                                str = sprintf('y1 = %.3f x + %.3f\n R2 = %.3f\n y2 = %.3f x + %.3f\n R2 = %.3f\n', tfitg1.p1, tfitg1.p2, Rg1(2)^2, tfitg2.p1, tfitg2.p2, Rg2(2)^2);
                                dim = [.15 .6 .3 .3];
                                annotation('textbox',dim,'String',str,'FitBoxToText','on');
                                legend('Data MAL','Trend MAL','Data CTL','Trend CTL')
                                xlabel(labelfactors{a},'fontsize',14)
                                ylabel(strrep(labeldata{1},'_',' '),'fontsize',14)
                                title(['Mean ' strrep(labeldata{1},'_',' ') ' of ' strrep(labeldata{2},'_',' ') ' by ' labelfactors{a}],'fontsize',14)
                                fig = gcf;
                                savefig([savepathgraph labeleffects{f} 'effect_' labeldata{2} '.fig'])
                                exportgraphics(fig,[savepathgraph labeleffects{f} 'effect_' labeldata{2} '.png'])
                                close
                                clear Rg1 Rg2 str dim tfitg1 tfitg2 sesg1 sesg2 datag1 datag2 fig

                            elseif  min(mod(factors{a},1) == 0) == 0 && min(mod(factors{b},1) == 0) == 0 %fact1cont %fact2cont

                            elseif min(mod(factors{a},1) == 0) == 1 && min(mod(factors{b},0) == 1) == 1 %fact1cat %fact2cat
                                figure %Group*Sex factor 1 and factor 3 %à débogguer
                                meang1s1 = mean(data(factors{a} == 1 & factors{b} == 1),'omitnan');
                                meang1s2 = mean(data(factors{a} == 1 & factors{b} == 2),'omitnan');
                                meang2s1 = mean(data(factors{a} == 2 & factors{b} == 1),'omitnan');
                                meang2s2 = mean(data(factors{a} == 2 & factors{b} == 2),'omitnan');
                                X = categorical({'CTL','CHD'});
                                hold on
                                bar(X(1),[meang1s1 meang1s2]); %b1
                                bar(X(2), [meang2s1 meang2s2]); %b2
                                xlabel(labelfactors{a},'fontsize',14)
                                ylabel(strrep(labeldata{1},'_',' '),'fontsize',14)
                                %xlim([0 , 3])
                                title(['Mean ' strrep(labeldata{1},'_',' ') ' of ' strrep(labeldata{2},'_',' ') ' by ' labelfactors{a}],'fontsize',14)
                                hold off
                                box on
                                %legend Ctl CHD
                                fig = gcf;
                                savefig([savepathgraph labeleffects{f} 'effect_' labeldata{2} '.fig'])
                                exportgraphics(fig,[savepathgraph labeleffects{f} 'effect_' labeldata{2} '.png'])
                                close
                                clear Rg1 Rg2 str dim tfitg1 tfitg2 fig
                            
                            elseif  min(mod(factors{a},1) == 0) == 1 && min(mod(factors{b},1) == 0) == 1 %fact1cat %fact2cont
                                figure
                                sesg1 = factors{b}(factors{a} == 1);
                                datag1 = data(factors{a} == 1);
                                idx = ~or(isnan(sesg1), isnan(datag1));
                                tfitg1 = fit(sesg1(idx), double(datag1(idx)),'poly1');
                                sesg2 = factors{b}(factors{a} == 2);
                                datag2 = data(factors{a} == 2);
                                idx = ~or(isnan(sesg2), isnan(datag2));
                                tfitg2 = fit(sesg2(idx), double(datag2(idx)),'poly1');
                                hold on
                                plot(tfitg1, '-r',sesg1, datag1, 'or')
                                plot(tfitg2, '-b',sesg2, datag2, 'ob')
                                Rg1 = corrcoef(sesg1, datag1);
                                Rg2 = corrcoef(sesg2, datag2);
                                str = sprintf('y1 = %.3f x + %.3f\n R2 = %.3f\n y2 = %.3f x + %.3f\n R2 = %.3f\n', tfitg1.p1, tfitg1.p2, Rg1(2)^2, tfitg2.p1, tfitg2.p2, Rg2(2)^2);
                                dim = [.15 .6 .3 .3];
                                annotation('textbox',dim,'String',str,'FitBoxToText','on');
                                legend('Data MAL','Trend MAL','Data CTL','Trend CTL')
                                xlabel(labelfactors{b},'fontsize',14)
                                ylabel(strrep(labeldata{1},'_',' '),'fontsize',14)
                                title(['Mean ' strrep(labeldata{1},'_',' ') ' of ' strrep(labeldata{2},'_',' ') ' by ' labelfactors{b}],'fontsize',14)
                                fig = gcf;
                                savefig([savepathgraph labeleffects{f} 'effect_' labeldata{2} '.fig'])
                                exportgraphics(fig,[savepathgraph labeleffects{f} 'effect_' labeldata{2} '.png'])
                                close
                                clear Rg1 Rg2 str dim tfitg1 tfitg2 fig
                            end
                        end
                    end
                end
            end
        end
        clear n_sig
    end


    if numel(sz) == 2 && sz(2) ~= 1
        for i = 1:sz(2)
            switch type
                case 'inter'
            [pval(:,i),res] = anovan(real(data(:,i)),factors,'Continuous',2,'model','full','varnames',labelfactors,'display','off');
                case 'main'
            [pval(:,i),res] = anovan(real(data(:,i)),factors,'Continuous',2,'varnames',labelfactors,'display','off');
            end
            resi = res(2:end,2:end);
            resi(cellfun(@isempty,resi)) = {[0]};
            resi = cellfun(@double,resi);
            resdata{1,i} = array2table(resi,'RowNames',res(2:end,1),'VariableNames',res(1,2:end));
            clear i res resi
        end

        labelrows = resdata{1,1}.Properties.RowNames(1:end-2,1); %{resdata{2:end-2,1,1}}';
        for f = 1:numel(labelrows)
            isig(f,:) = pval(f,:) <= p;
            idsig{f,1} = find(isig(f,:));
            n_sig(f,:) = sum(isig(f,:));
            fprintf('\t \t %d %s effects are significant p<=%.2f without correction for %s\n', n_sig(f,:), labelrows{f,1}, p, labeldata{2})
        end

        if fdrmode == 1
            FDR_CO(pval, labelrows, labeldim, labeldata{2}, p, savepath)
        end

        tblres = array2table(resdata,'VariableName',labeldim);
        tblpval = array2table(pval,'VariableName',labeldim,'RowNames',labelrows');

        if savemode
            save([savepath date '_' labeldata{2} 'ResultsANOVA' type '.mat'],'tblres','tblpval');

            for f = 1:numel(labelrows)
                fac = labelrows{f,1};
                if contains(fac,'*')
                    fac = erase(fac,'*');
                elseif contains(fac,':')
                    fac = erase(fac,':');
                end
                sig = tblpval{f,:} <= p;
                if find(sig)
                    tblpvalsig.(fac) = tblpval(f,sig);
                    tblressig.(fac) = tblres(:,sig);
                    disp(tblpvalsig.(fac))
                    save([savepath date '_' labeldata{2} 'ResultsANOVA' type '.mat'],'tblpvalsig','tblressig','-append');
                else
                end
            end
        end

        if graphmode
            for f = 1:numel(labeleffects)
                fac = labeleffects{f};
                if n_sig(f,:) > 0 && f<=nfactors %principal effect
                    if min(mod(factors{f},1) == 0) == 0 %continuous
                        figure
                        hold on
                        for i = 1:n_sig(f,:)
                            idx = ~or(isnan(factors{f}), isnan(data(:,idsig{f,1}(i))));
                            tfit = fit(factors{f}(idx), double(data(idx,idsig{f,1}(i))),'poly1');
                            slopes.(fac)(i) = tfit.p1;
                            if tfit.p1 > 0
                                plot(tfit,'r')
                            elseif tfit.p1 < 0
                                plot(tfit,'b')
                            end
                        end
                        xlabel(labelfactors{f},'fontsize',14)
                        ylabel(strrep(labeldata{1},'_',' '),'fontsize',14)
                        title(['Mean ' strrep(labeldata{1},'_',' ') ' of ' strrep(labeldata{2},'_',' ') ' by ' labelfactors{f}],'fontsize',14)
                        legend off
                        hold off
                        %box on
                        %legend(labeldim(isig))
                        fig = gcf;
                        savefig([savepathgraph labeleffects{f} 'effect_' labeldata{2} '.fig'])
                        exportgraphics(fig,[savepathgraph labeleffects{f} 'effect_' labeldata{2} '.png'])
                        close
                        clear tfit i idx fig
                    elseif  min(mod(factors{f},1) == 0) == 1 %categorical
                        figure
                        hold on
                        meang1 = mean(data(factors{f} == 1,:),'omitnan');
                        meang2 = mean(data(factors{f} == 2,:),'omitnan');
                        D = [meang1(isig(f,:))' meang2(isig(f,:))'];
                        X = categorical(labeldim(isig(f,:)));
                        b = bar(X,D);
                        b(1).FaceColor = 'r';
                        b(2).FaceColor = 'b';
                        xlabel(labelfactors{f},'fontsize',14)
                        ylabel(strrep(labeldata{1},'_',' '),'fontsize',14)
                        title(['Mean ' strrep(labeldata{1},'_',' ') ' of ' strrep(labeldata{2},'_',' ') ' by ' labelfactors{f}],'fontsize',14)
                        hold off
                        box on
                        legend Mal Ctl
                        fig = gcf;
                        fig.WindowState = 'maximized';
                        savefig([savepathgraph labeleffects{f} 'effect_' labeldata{2} '.fig'])
                        exportgraphics(fig,[savepathgraph labeleffects{f} 'effect_' labeldata{2} '.png'])
                        close
                        clear D X b fig
                    end

                elseif n_sig(f,:) > 0 && f>nfactors %interaction effect
                    a = lgndinter{f}(1); b = lgndinter{f}(2);

                    if min(mod(factors{a},1) == 0) == 0 && min(mod(factors{b},1) == 0) == 1 %fact1cont %fact2cat
                        b = lgndinter{f}(1); a = lgndinter{f}(2);
                    elseif  min(mod(factors{a},1) == 0) == 1 && min(mod(factors{b},1) == 0) == 0 %fact1cat %fact2cont
                        figure %Group*SES %factor 1 and factor2
                        hold on
                        sesg1 = factors{b}(factors{a} == 1);
                        datag1 = data(factors{a} == 1,:);
                        for i = 1:n_sig(f,:)
                            idx = ~(isnan(sesg1) | isnan(datag1(:,idsig{f,1}(i))));
                            tfitg1 = fit(sesg1(idx), double(datag1(idx,idsig{f,1}(i))),'poly1');
                            slopes.(fac)(1,i) = tfitg1.p1;
                            if slopes.(fac)(1,i)>0
                                p1 = plot(tfitg1, '-r');
                            else
                                p2 = plot(tfitg1, '--r');
                            end
                        end
                        sesg2 = factors{b}(factors{a} == 2);
                        datag2 = data(factors{a} == 2,:);
                        for i = 1:n_sig(f,:)
                            idx = ~(isnan(sesg2) | isnan(datag2(:,idsig{f,1}(i))));
                            tfitg2 = fit(sesg2(idx), double(datag2(idx,idsig{f,1}(i))),'poly1');
                            slopes.(fac)(2,i) = tfitg2.p1;
                            if slopes.(fac)(2,i)>0
                                p3 = plot(tfitg2, '-b');
                            else
                                p4 = plot(tfitg2, '--b');
                            end
                        end
                        h = [];
                        idx = [];
                        if exist('p1','var')
                            h = [h ; p1];
                            idx = [idx; 1];
                        end
                        if exist('p2','var')
                            h = [h ; p2];
                            idx = [idx; 2];
                        end
                        if exist('p3','var')
                            h = [h ; p3];
                            idx = [idx; 3];
                        end
                        if exist('p4','var')
                            h = [h ; p4];
                            idx = [idx; 4];
                        end
                        str = {'Trend MAL positive';'Trend MAL negative'; 'Trend CTL positive'; 'Trend CTL negative'};
                        legend(h, str(idx))
                        xlabel(labelfactors{b},'fontsize',14)
                        ylabel(strrep(labeldata{1},'_',' '),'fontsize',14)
                        title([strrep(labeldata{1},'_',' ') ' of ' strrep(labeldata{2},'_',' ') ' by ' labelfactors{b}],'fontsize',14)
                        fig = gcf;
                        savefig([savepathgraph labeleffects{f} 'effect_' labeldata{2} '.fig'])
                        exportgraphics(fig,[savepathgraph labeleffects{f} 'effect_' labeldata{2} '.png'])
                        close
                        clear tfitg1 tfitg2 sesg1 sesg2 datag1 datag2 idx ifig

                    elseif  min(mod(factors{a},1) == 0) == 0 && min(mod(factors{b},1) == 0) == 0 %fact1cont %fact2cont

                    elseif min(mod(factors{a},1) == 0) == 1 && min(mod(factors{b},1) == 0) == 1 %fact1cat %fact2cat
                        %Group*Sex factor 1 and factor 3 %à débogguer
                    end
                end
            end
            clear f n_sig
            if exist('slopes','var')
                varargout{1} = slopes;
                if savemode
                    save([savepath date '_' labeldata{2} 'ResultsANOVA' type '.mat'],'slopes','-append');
                end
            else
                slopes = [];
                varargout{1} = slopes;
            end
        end
    end

    if numel(sz) == 3
        for j = 1:sz(1)
            for i = 1:sz(3)
                switch type
                    case 'inter'
                [pval(:,i,j),res] = anovan(real(data(j,:,i)),factors,'Continuous',2,'model','full','varnames',labelfactors,'display','off');
                    case 'main'
                [pval(:,i,j),res] = anovan(real(data(j,:,i)),factors,'Continuous',2,'varnames',labelfactors,'display','off');
                end

                resi = res(2:end,2:end);
                resi(cellfun(@isempty,resi)) = {single([0])};
                resi = cellfun(@double,resi);
                resdatai{1,i} = array2table(resi,'RowNames',res(2:end,1),'VariableNames',res(1,2:end));
                clear res resi
            end
            resdata(j,:) = resdatai;
            clear resdatai
        end

        labelrows = resdata{1,1}.Properties.RowNames(1:end-2,1); %{resdata{2:end-2,1,1}}';
        for f = 1:numel(labelrows)
            n_sig(f,:) = sum(squeeze(pval(f,:,:) <= p));
            fprintf('\t \t Mean of %.1f significant p<=%.2f %s effects per channel without correction for %s\n',mean(n_sig(f,:),2), p, labelrows{f,1}, labeldata{2})
        end

        if graphmode
            fig = figure;
            hold on
            for f = 1:numel(labelrows)
                plot(n_sig(f,:));
            end
            xlabel('Channel number')
            ylabel('Nb sig effect')
            legend(labelrows)
            savefig([savepathgraph 'resultsANOVA' type '_' labeldata{2} '.fig'])
            exportgraphics(fig,[savepathgraph 'resultsANOVA' type '_' labeldata{2} '.png'])
            close
        end
        clear n_sig

        if fdrmode == 1
            FDR_CO(pval, labelrows, labeldim, labeldata{2}, p, savepath)
        end

        tblres = array2table(resdata,'VariableName',labeldim{2}, 'RowNames', labeldim{1});
        %tblpval = table(pval,'RowNames',labelrows');
        for n = 1:length(labelrows)
            for i = 1:sz(1)
                tblpval{n,i} = squeeze(pval(n,:,i));
            end
        end
        tblpval = array2table(tblpval,'RowNames',labelrows','VariableNames',labeldim{1});

        if savemode
            save([savepath date '_' labeldata{2} 'ResultsANOVA' type '.mat'],'tblres','tblpval');
        end
    end

    % if numel(labelrows) > nfactors %nargin == 19
    %     idx = find(contains(labelrows,':')); %trouver le num du facteur
    %     idsiginter = find(pval(idx,:)<= p); %trouver les interactions sig
    %
    %     tblcoef = array2table(varargin{1},'VariableName',label);
    %     tblcoefsiginterG1 = array2table(varargin{2},'VariableNames',label(idsiginter));
    %     tblcoefsiginterG2 = array2table(varargin{3},'VariableNames',label(idsiginter));
    %     tblpvalsiginterG1 = array2table(varargin{4},'VariableNames',label(idsiginter));
    %     tblpvalsiginterG2 = array2table(varargin{5},'VariableNames',label(idsiginter));
    %     save([savepath date '_' labeldata 'Results.mat'],'tblcoef', 'tblcoefsiginterG1', 'tblcoefsiginterG2', 'tblpvalsiginterG1', 'tblpvalsiginterG2','-append');
    % end

catch
    fprintf('Error during analyses for %s\n', labeldata{2})
    tblpval = [];
    tblres = [];
end
end