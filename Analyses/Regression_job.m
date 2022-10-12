function [tblpvalres, tblres, tblcoef, tblintersig] = Regression_job(type, data, factors, labelfactors,...
    p, labeldata, labeldim, fdrmode, graphmode, savepath, savemode)

%type = 'component' OR 'summary'
%savepath = fullfile(savepath,'Regression',filesep);

%     if graphmode
%         savepathgraph = fullfile(savepath,'Graphs',filesep);
%         if ~isfolder(savepathgraph)
%             mkdir(savepathgraph)
%         end
%     end

nfactors = numel(factors);
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

sz = size(data);

if numel(sz) == 2 && sz(2) ~= 1
    for i = 1:sz(2)
        mdl = fitlm([factors{:}],data(:,i),'interactions','CategoricalVars',1,'VarNames',{'Group','SES','Connectivity'});
        coef = mdl.Coefficients;
        coefdata{1,i} = coef;

        res = anova(mdl,type);
        resdata{1,i} = res;
        pvalres(:,i) = res.pValue(1:end-1);
        clear mdl res coef
    end

    labelrows = resdata{1,1}.Properties.RowNames(1:end-1,1);
    for f = 1:numel(labelrows)
        isig(f,:) = pvalres(f,:) <= p;
        idsig{f,1} = find(isig(f,:));
        n_sig(f,:) = sum(pvalres(f,:) <= p);
        fprintf('\t \t %d %s effects are significant p<=%.2f without correction for %s\n', n_sig(f,:), labelrows{f,1}, p, labeldata{2})
    end

    clear i f %n_sig

    if fdrmode == 1
        FDR_CO(pvalres, labelrows, labeldim, labeldata{2}, p, savepath)
    end

    tblres = array2table(resdata,'VariableName',labeldim);
    tblcoef = array2table(coefdata,'VariableName',labeldim);
    tblpvalres = array2table(pvalres,'VariableName',labeldim,'RowNames',labelrows');

    if savemode
        save([savepath date '_' labeldata{2} 'ResultsReg' type '.mat'],'tblres','tblpvalres','tblcoef');

        for f = 1:numel(labelrows)
            fac = labelrows{f,1};
            if contains(fac,'*')
                fac = erase(fac,'*');
            elseif contains(fac,':')
                fac = erase(fac,':');
            end
            sig = tblpvalres{f,:} <= p;
            if find(sig)
                tblpvalsig.(fac) = tblpvalres(f,sig);
                tblressig.(fac) = tblres(:,sig);
                disp(tblpvalsig.(fac))
                save([savepath date '_' labeldata{2} 'ResultsReg' type '.mat'],'tblpvalsig','tblressig','-append');
            else
            end
        end
    end

    for f = 1:numel(labeleffects) %for each effect tested (simple effects and two sided interaction)
        fac = labeleffects{f}; %find effect name
        if n_sig(f,:) > 0 && f>nfactors %if significant and interaction effect
            a = lgndinter{f}(1); b = lgndinter{f}(2);
            if min(mod(factors{a},1) == 0) == 0 && min(mod(factors{b},1) == 0) == 1 %fact1cont %fact2cat
                b = lgndinter{f}(1); a = lgndinter{f}(2);
            
            elseif  min(mod(factors{a},1) == 0) == 1 && min(mod(factors{b},1) == 0) == 0 %fact1cat %fact2cont

                %%Slopes differences between level of SES
                p25 = prctile(factors{b},25);
                high = factors{b} - p25;
                p50 = prctile(factors{b},50);
                cr = factors{b} - p50;
                p75 = prctile(factors{b},75);
                low = factors{b} - p75;

                for j = 1:numel(idsig{f,:})
                    i = idsig{f}(j);
                    mdl = fitlm([factors{a} high],data(:,i),'interactions','CategoricalVars',1,'VarNames',{'Group','SEShigh','Connectivity'});
                    coef = mdl.Coefficients;
                    coefdatahigh{1,j} = coef;
                    pvalcoefhigh(:,j) = coef.pValue;
                    res = anova(mdl,type);
                    resdatahigh{1,j} = res;
                    pvalreshigh(:,j) = res.pValue(1:end-1);
                    clear mdl res coef

                    mdl = fitlm([factors{a} cr],data(:,i),'interactions','CategoricalVars',1,'VarNames',{'Group','SEScr','Connectivity'});
                    coef = mdl.Coefficients;
                    coefdatacr{1,j} = coef;
                    pvalcoefcr(:,j) = coef.pValue;
                    res = anova(mdl,type);
                    resdatacr{1,j} = res;
                    pvalrescr(:,j) = res.pValue(1:end-1);
                    clear mdl res coef

                    mdl = fitlm([factors{a} low],data(:,i),'interactions','CategoricalVars',1,'VarNames',{'Group','SESlow','Connectivity'});
                    coef = mdl.Coefficients;
                    coefdatalow{1,j} = coef;
                    pvalcoeflow(:,j) = coef.pValue;
                    res = anova(mdl,type);
                    resdatalow{1,j} = res;
                    pvalreslow(:,j) = res.pValue(1:end-1);
                    clear mdl res coef
                end
                clear i j
                
                tblintersig.pvalgrcoef = array2table([pvalcoeflow(a+1,:); pvalcoefcr(a+1,:);  pvalcoefhigh(a+1,:)], 'VariableNames', labeldim(idsig{f,:}), 'RowNames',{'PvalLow','PvalCr','PvalHigh'});
                tblintersig.res = array2table([resdatalow; resdatacr; resdatahigh],'VariableNames', labeldim(idsig{f,:}), 'RowNames', {'ResLow','ResCr','ResHigh'});
                tblintersig.coef = array2table([coefdatalow; coefdatacr; coefdatahigh],'VariableNames', labeldim(idsig{f,:}), 'RowNames', {'CoefLow','CoefCr','CoefHigh'});
                tblintersig.pvalcoef.Low = array2table(pvalcoeflow, 'VariableNames', labeldim(idsig{f,:}), 'RowNames',coefdatalow{1}.Row);
                tblintersig.pvalcoef.Cr = array2table(pvalcoefcr, 'VariableNames', labeldim(idsig{f,:}), 'RowNames',coefdatacr{1}.Row);
                tblintersig.pvalcoef.High = array2table(pvalcoefhigh, 'VariableNames', labeldim(idsig{f,:}), 'RowNames',coefdatahigh{1}.Row);
                tblintersig.pvalres.Low = array2table(pvalreslow, 'VariableNames', labeldim(idsig{f,:}), 'RowNames',resdatalow{1}.Row(1:end-1));
                tblintersig.pvalres.Cr = array2table(pvalrescr, 'VariableNames', labeldim(idsig{f,:}), 'RowNames',resdatacr{1}.Row(1:end-1));
                tblintersig.pvalres.High = array2table(pvalreshigh, 'VariableNames', labeldim(idsig{f,:}), 'RowNames',resdatahigh{1}.Row(1:end-1));

                %%Slopes differences between groups
                for j = 1:numel(idsig{f,:}) %pour chaque inter sig
                    i = idsig{f}(j);
                    mdl = fitlm(factors{b}(factors{a}==1),data(factors{a}==1,i),'VarNames',{'SES','Connectivity'});
                    coefG1 = mdl.Coefficients;
                    coefdataG1{1,j} = coefG1;
                    pvalcoefG1(:,j) = coefG1.pValue;
                    clear mdl coefG1
                    mdl = fitlm(factors{b}(factors{a}==2),data(factors{a}==2,i),'VarNames',{'SES','Connectivity'});
                    coefG2 = mdl.Coefficients;
                    coefdataG2{1,j} = coefG2;
                    pvalcoefG2(:,j) = coefG2.pValue;
                    clear mdl coefG2
                end
                clear i j

                tblintersig.coefG = array2table([coefdataG1; coefdataG2],'VariableNames', labeldim(idsig{f,:}), 'RowNames', {'CoefG1','CoefG2'});
                tblintersig.pvalcoef.G1 = array2table(pvalcoefG1, 'VariableNames', labeldim(idsig{f,:}), 'RowNames',coefdataG1{1}.Row);
                tblintersig.pvalcoef.G2 = array2table(pvalcoefG2, 'VariableNames', labeldim(idsig{f,:}), 'RowNames',coefdataG2{1}.Row);

                varargout{1} = tblintersig;
                if savemode
                    save([savepath date '_' labeldata{2} 'ResultsReg' type '.mat'],'tblintersig','-append');
                end
            end
        end
    end
end
end