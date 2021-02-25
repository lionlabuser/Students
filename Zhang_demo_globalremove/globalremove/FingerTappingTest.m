% global component removal made by Dr. Xian Zhang from Yale school of medicine
load testdata 
nCh=size(meanxyz,1);
polhemus.NFRI_result.OtherC=meanxyz;
polhemus.nTR=0; % if you have xyz containing opto location, nTR= all optode number
polhemus.result4=zeros(nCh,1);
% just follow the above three line. polhemus is a struction used in our pipe line. you only need these lines here
badCh=[1 3 5] ; % bad channel is 1 3 5
badCh=[]; % no bad channel
n=size(meandata,1); % the number of sample
for i=1:2
data2(:,i,:)=reshape(meandata(:,:,i),n,1,nCh);
%data  should be a matrix first dimension is data points, 
% second dim: 1 is  oxy 2 is dexpoy, 
% third dim is channel 
end
deg=1.3  % the angle of the spatial filter is   sqrt(deg/2)*180/pi
data2(isnan(data2))=0;
data3=removeglobalY(data2,polhemus,deg,badCh);
nirs_data.oxyData=squeeze(data3(:,1,:));
nirs_data.dxyData=squeeze(data3(:,2,:));
%-data3+data2;  %is the global 
figure;plot(mean(meandata(:,:,1),2));
hold on;plot(mean(nirs_data.oxyData(:,:,1),2));
% only a few channel should show finger tap activit, but the  average of raw data show clear finger tapping activity,
%while the clean data does not show finger tapping activity
