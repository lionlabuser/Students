
dirmtg='C:\Data\data_NIRS\ELAN\ANALYSED\Martine_0m\Montage\';
load( [dirmtg 'Global.zone'] ,'-mat'); %zone file
load([dirmtg 'Casque36_ELAN.prj'],'-mat') %prf file

%% help info zone file
%zone.pos(channel#,[x pos, y pos, z pos, channel distance in CM])
%zone.ml (channel#, [Source#, Detector#, ?, HbO(1) or HbR(2)]  

%% create a matrice channel# x channel# to get the distance between each pair of channels

nchan=size(zone.pos,1)/2; %number of channels
   
ChanDistanceBTW=zeros(nchan);
for a=1:nchan
    for b=a+1:nchan
        
        ChanDistanceBTW(a,b)=sqrt( (zone.pos(a,1)-zone.pos(b,1))^2 + (zone.pos(a,2)-zone.pos(b,2))^2+ (zone.pos(a,3)-zone.pos(b,3))^2);
        ChanDistanceBTW(b,a)=ChanDistanceBTW(a,b);
        
        rayon(a,b)=sqrt( zone.pos(a,1)^2 + zone.pos(a,2)^2 + zone.pos(a,3)^2);
        rayon(b,a)=rayon(a,b);
        distanceARC(a,b)=2*rayon(a,b)*asin(ChanDistanceBTW(a,b)/(2*rayon(a,b)));
        distanceARC(b,a)=distanceARC(a,b);
    end
    ChanDistance(a)=zone.pos(a,4); 
    ChanXYZ(a,:)=zone.pos(a,1:3);
end

distanceARC=real(distanceARC);
mrayon=mean(rayon,'all','omitnan');

ChanNum=1:nchan;
ChanHelp=sprintf('ChanNum = Channel ID (#)\nChanDistance = length of each channel in centimeters (source-detector distance)\nChanDistanceBTW = matrice of distance between the midpoint of one channel to the midpoint of another channel (in cm)\nChanXYZ =  x y z coordinates of each channel (in cm)');

%save([dirmtg 'ChannelsDistance.mat'],'Chan*')

