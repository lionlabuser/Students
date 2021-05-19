tmp= PARCOMP(1)
Xconstante = tmp.AUX.data{1};
XHRF = tmp.AUX.data{2};
plot(XHRF )
plot(XHRF )
for ich = 1:size(PARCOMP.data,2)
    y = tmp.data(:,ich);
    [b,bint,r,rint,stats] = regress(y,[XHRF]); %Xconstante,
    
    R2(ich) = stats(1)
end

mean(R2)