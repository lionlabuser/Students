%%%%%%%%%%%%%%RANCOVA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
disp('Computing RANCOVA')

data = 'C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Stats\CORR\0,01_0,08\';
load ([data 'workspace.mat'])
%load ([data 'workspacemat.mat'])

savepath='C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Stats\CORR\0,01_0,08\rANCOVA\';
channelmode = 1;
tablegraphmode = 1;
fdrmode = 1;
p = 0.07;

%%CHANNELS%%%%%%%%%%%%%%%%%%%%%
x = 1;
for c = 1:46
    for cc = (c+1):46
        if c <= 23 && cc <= 23
            type{x,1} = 'IntraL';
            x = x + 1;
        elseif c >23 && cc > 23
            type{x,1} = 'IntraR';
            x = x + 1;
        else 
            type{x,1} = 'Inter';
            x = x + 1;
        end
    end
end

withinch = array2table(type,'VariableNames',{'type'});

x = 1;
for c = 1:1035
    labels{x} = ['chp' num2str(c)];
    x = x + 1;
end
betweench = tblch(:,[2 5 6:end]);  
betweench.Properties.VariableNames = {'Group', 'SSE', labels{1,1:end}};

rm = fitrm(betweench,'chp1-chp1035 ~ Group*SSE','WithinDesign',withinch);
rancovatblch = ranova(rm);

clear type labels rm c cc

%%MROI%%%%%%%%%%%%%%%%%%%%%%%
x = 1;
for r = 1:numel(ListMroi)
    for rr = (r+1):numel(ListMroi)
        if r <= numel(ListMroi)/2 && rr <= numel(ListMroi)/2
            type{x,1} = 'IntraL';
            x = x + 1;
        elseif r >numel(ListMroi)/2 && rr > numel(ListMroi)/2
            type{x,1} = 'IntraR';
            x = x + 1;
        else 
            type{x,1} = 'Inter';
            x = x + 1;
        end
    end
end

withinMroi = array2table(type,'VariableNames',{'type'});

x = 1;
for r = 1:91
    labels{x} = ['Mroip' num2str(r)];
    x = x + 1;
end
betweenMroi = tblmeanMroi(:,[2 5 6:end]);  
betweenMroi.Properties.VariableNames = {'Group', 'SSE', labels{1,1:end}};

rm = fitrm(betweenMroi,'Mroip1-Mroip91 ~ Group*SSE','WithinDesign',withinMroi);
rancovatblMroi = ranova(rm);

clear type labels rm r rr

%%RROI%%%%%%%%%%%%%%%%%%%%%%%
x = 1;
for r = 1:numel(ListRroi)
    for rr = (r+1):numel(ListRroi)
        if r <= numel(ListRroi)/2 && rr <= numel(ListRroi)/2
            type{x,1} = 'IntraL';
            x = x + 1;
        elseif r >numel(ListRroi)/2 && rr > numel(ListRroi)/2
            type{x,1} = 'IntraR';
            x = x + 1;
        else 
            type{x,1} = 'Inter';
            x = x + 1;
        end
    end
end

withinRroi = array2table(type,'VariableNames',{'type'});

x = 1;
for r = 1:width(tblmeanRroi)- 5
    labels{x} = ['Rroip' num2str(r)];
    x = x + 1;
end
betweenRroi = tblmeanRroi(:,[2 5 6:end]);  
betweenRroi.Properties.VariableNames = {'Group', 'SSE', labels{1,1:end}};

rm = fitrm(betweenRroi,'Rroip1-Rroip325 ~ Group*SSE','WithinDesign',withinRroi);
rancovatblRroi = ranova(rm);

clear type labels rm r rr

%%FROI%%%%%%%%%%%%%%%%%%%%%%%
x = 1;
for r = 1:numel(ListFroi)
    for rr = (r+1):numel(ListFroi)
        if r <= numel(ListFroi)/2 && rr <= numel(ListFroi)/2
            type{x,1} = 'IntraL';
            x = x + 1;
        elseif r >numel(ListFroi)/2 && rr > numel(ListFroi)/2
            type{x,1} = 'IntraR';
            x = x + 1;
        else 
            type{x,1} = 'Inter';
            x = x + 1;
        end
    end
end

withinFroi = array2table(type,'VariableNames',{'type'});

x = 1;
for r = 1:width(tblmeanFroi)- 5
    labels{x} = ['Froip' num2str(r)];
    x = x + 1;
end
betweenFroi = tblmeanFroi(:,[2 5 6:end]);  
betweenFroi.Properties.VariableNames = {'Group', 'SSE', labels{1,1:end}};

rm = fitrm(betweenFroi,'Froip1-Froip66 ~ Group*SSE','WithinDesign',withinFroi);
rancovatblFroi = ranova(rm);
clear type labels rm r rr

%%AROI%%%%%%%%%%%%%%%%%%%%%%%
x = 1;
for r = 1:numel(ListAroi)
    for rr = (r+1):numel(ListAroi)
        if r <= numel(ListAroi)/2 && rr <= numel(ListAroi)/2
            type{x,1} = 'IntraL';
            x = x + 1;
        elseif r >numel(ListAroi)/2 && rr > numel(ListAroi)/2
            type{x,1} = 'IntraR';
            x = x + 1;
        else 
            type{x,1} = 'Inter';
            x = x + 1;
        end
    end
end

withinAroi = array2table(type,'VariableNames',{'type'});

x = 1;
for r = 1:width(tblmeanAroi)- 5
    labels{x} = ['Aroip' num2str(r)];
    x = x + 1;
end
betweenAroi = tblmeanAroi(:,[2 5 6:end]);  
betweenAroi.Properties.VariableNames = {'Group', 'SSE', labels{1,1:end}};

rm = fitrm(betweenAroi,'Aroip1-Aroip28 ~ Group*SSE','WithinDesign',withinAroi);
rancovatblAroi = ranova(rm);

clear type labels rm r rr
save([savepath date '_resultsRANCOVA.mat'],'rancovatblch','rancovatblMroi','rancovatblRroi','rancovatblFroi','rancovatblAroi');

X = ['Results of RANCOVA saved in ', savepath];
disp(X)
clear X

%%%%%%%%%%RANOVA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Computing RANOVA')
savepath='C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\Stats\CORR\0,01_0,08\rANOVA\';

%%CHANNELS%%%%%%%%%%%%%%%%%%%%%
x = 1;
for c = 1:46
    for cc = (c+1):46
        if c <= 23 && cc <= 23
            type{x,1} = 'IntraL';
            x = x + 1;
        elseif c >23 && cc > 23
            type{x,1} = 'IntraR';
            x = x + 1;
        else 
            type{x,1} = 'Inter';
            x = x + 1;
        end
    end
end

withinch = array2table(type,'VariableNames',{'type'});

x = 1;
for c = 1:1035
    labels{x} = ['chp' num2str(c)];
    x = x + 1;
end
betweench = tblch(:,[2 6:end]);  
betweench.Properties.VariableNames = {'Group', labels{1,1:end}};

rm = fitrm(betweench,'chp1-chp1035 ~ Group','WithinDesign',withinch);
ranovatblch = ranova(rm);

clear type labels rm c cc

%%MROI%%%%%%%%%%%%%%%%%%%%%%%
x = 1;
for r = 1:numel(ListMroi)
    for rr = (r+1):numel(ListMroi)
        if r <= numel(ListMroi)/2 && rr <= numel(ListMroi)/2
            type{x,1} = 'IntraL';
            x = x + 1;
        elseif r >numel(ListMroi)/2 && rr > numel(ListMroi)/2
            type{x,1} = 'IntraR';
            x = x + 1;
        else 
            type{x,1} = 'Inter';
            x = x + 1;
        end
    end
end

withinMroi = array2table(type,'VariableNames',{'type'});

x = 1;
for r = 1:91
    labels{x} = ['Mroip' num2str(r)];
    x = x + 1;
end
betweenMroi = tblmeanMroi(:,[2 6:end]);  
betweenMroi.Properties.VariableNames = {'Group', labels{1,1:end}};

rm = fitrm(betweenMroi,'Mroip1-Mroip91 ~ Group','WithinDesign',withinMroi);
ranovatblMroi = ranova(rm);

clear type labels rm r rr

%%RROI%%%%%%%%%%%%%%%%%%%%%%%
x = 1;
for r = 1:numel(ListRroi)
    for rr = (r+1):numel(ListRroi)
        if r <= numel(ListRroi)/2 && rr <= numel(ListRroi)/2
            type{x,1} = 'IntraL';
            x = x + 1;
        elseif r >numel(ListRroi)/2 && rr > numel(ListRroi)/2
            type{x,1} = 'IntraR';
            x = x + 1;
        else 
            type{x,1} = 'Inter';
            x = x + 1;
        end
    end
end

withinRroi = array2table(type,'VariableNames',{'type'});

x = 1;
for r = 1:width(tblmeanRroi)- 5
    labels{x} = ['Rroip' num2str(r)];
    x = x + 1;
end
betweenRroi = tblmeanRroi(:,[2 6:end]);  
betweenRroi.Properties.VariableNames = {'Group', labels{1,1:end}};

rm = fitrm(betweenRroi,'Rroip1-Rroip325 ~ Group','WithinDesign',withinRroi);
ranovatblRroi = ranova(rm);

clear type labels rm r rr

%%FROI%%%%%%%%%%%%%%%%%%%%%%%
x = 1;
for r = 1:numel(ListFroi)
    for rr = (r+1):numel(ListFroi)
        if r <= numel(ListFroi)/2 && rr <= numel(ListFroi)/2
            type{x,1} = 'IntraL';
            x = x + 1;
        elseif r >numel(ListFroi)/2 && rr > numel(ListFroi)/2
            type{x,1} = 'IntraR';
            x = x + 1;
        else 
            type{x,1} = 'Inter';
            x = x + 1;
        end
    end
end

withinFroi = array2table(type,'VariableNames',{'type'});

x = 1;
for r = 1:width(tblmeanFroi)- 5
    labels{x} = ['Froip' num2str(r)];
    x = x + 1;
end
betweenFroi = tblmeanFroi(:,[2 6:end]);  
betweenFroi.Properties.VariableNames = {'Group', labels{1,1:end}};

rm = fitrm(betweenFroi,'Froip1-Froip66 ~ Group','WithinDesign',withinFroi);
ranovatblFroi = ranova(rm);
clear type labels rm r rr

%%AROI%%%%%%%%%%%%%%%%%%%%%%%
x = 1;
for r = 1:numel(ListAroi)
    for rr = (r+1):numel(ListAroi)
        if r <= numel(ListAroi)/2 && rr <= numel(ListAroi)/2
            type{x,1} = 'IntraL';
            x = x + 1;
        elseif r >numel(ListAroi)/2 && rr > numel(ListAroi)/2
            type{x,1} = 'IntraR';
            x = x + 1;
        else 
            type{x,1} = 'Inter';
            x = x + 1;
        end
    end
end

withinAroi = array2table(type,'VariableNames',{'type'});

x = 1;
for r = 1:width(tblmeanAroi)- 5
    labels{x} = ['Aroip' num2str(r)];
    x = x + 1;
end
betweenAroi = tblmeanAroi(:,[2 6:end]);  
betweenAroi.Properties.VariableNames = {'Group', labels{1,1:end}};

rm = fitrm(betweenAroi,'Aroip1-Aroip28 ~ Group','WithinDesign',withinAroi);
ranovatblAroi = ranova(rm);

clear type labels rm r rr

save([savepath date '_resultsRANOVA.mat'],'ranovatblch','ranovatblMroi','ranovatblRroi','ranovatblFroi','ranovatblAroi');

X = ['Results of RANOVA saved in ', savepath];
disp(X)
clear X