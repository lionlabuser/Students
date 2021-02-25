% ROCanalysis/allsubject combine data
GLOBALREMOVALR=1.3;
x.ProcessList=[0]
xch=fNIRSAnalysis_ChannelWisePipeLine(x,0.03);
xch=xch.getdata;
meanxyz=  xch.meanxyz;
R=0;
r=1;c=1;
clear dataall  meandata
for i=1:x.nSubj
l=size( xch.data{r,c,i,side,k,R+1},1);
for k=1:2
dataall(1:min([l 2000]),:,i,k)= xch.data{r,c,i,side,k,R+1}(1:min([l 2000]),:);
end
end
meandata=squeeze(mean(dataall,3));
save testdata  meandata meanxyz
