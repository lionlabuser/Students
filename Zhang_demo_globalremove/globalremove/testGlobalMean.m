nCh=40;
xyz=rand(nCh,3);
polhemus.NFRI_result.OtherC=xyz;
polhemus.result4=zeros(nCh,1);
badCh=zoers(nCh,1);
badCh(1)=1; % for example remove ch 1
bacCh=[]; % n o back channel
bacCh=[1 3 5] ; % 

polhemus.nTR=0; % if you have xyz containing opto location, nTR= all optode number

%data  should be a matrix first dimension is data points, 
% second dim: 1 is  oxy 2 is dexpoy, 
% third dim is channel 
% eg. for our data format
% data0 is our OMM data 
data=data0(:,4:end)
data2(:,1,:)=reshape(data(:,1:3:end),n,1,nch);
data2(:,2,:)=reshape(data(:,2:3:end),n,1,nch);
% at the end the data2 should size  :  1200 ,2,40 if you aquire 1200 data samples
deg=1.3
data2(isnan(data2))=0;
data3=removeglobalY(data2,polhemus,deg,badCh);
nirs_data.oxyData=squeeze(data3(:,1,:));
nirs_data.dxyData=squeeze(data3(:,2,:));
-data3+data2;  %is the global 

