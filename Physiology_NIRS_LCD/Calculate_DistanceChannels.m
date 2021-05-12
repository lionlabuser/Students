
% please, before using the script, copy-paste it somewhere ELSE than the Student Directory (as it refers to the original version of the script

% ChanNum = Channel ID (#)
% ChanSDseparation = length of each channel in centimeters (source-detector separation)
% ChanDistanceLIN = matrice of the linear-direct distance between the midpoint of one channel to the midpoint of another channel (in cm)
% ChanDistanceARC = matrice of the arc distance between the midpoint of one channel to the midpoint of another channel (in cm). The circle is centered at x=0 y=0 z=0
% rayon = In relation to ChanDistanceARC, the radius (rayon) of the circle taken to calculate the distance is stored in <rayon> matrice. It does vary from one pair to the other
% mrayon = average of radius (rayon matrice)
% ChanXYZ =  x y z coordinates of each channel (in cm)

%% to adapt according to your localization
dirmtg='C:\Data\data_NIRS\ELAN\ANALYSED\Martine_0m\Montage\';
load( [dirmtg 'Global.zone'] ,'-mat'); %zone file
load([dirmtg 'Casque36_ELAN.prj'],'-mat') %prf file

%% help info zone file
%zone.pos(channel#,[x pos, y pos, z pos, channel distance in CM])
%zone.ml (channel#, [Source#, Detector#, ?, HbO(1) or HbR(2)]  

% create a matrice channel# x channel# to get the linear distance between each pair of channels
% create a matrice channel# x channel# to get the curvilinear (arc)
% distance between each pair of channels (center at XYZ=0,0,0)


nchan=size(zone.pos,1)/2; %number of channels
   
ChanDistanceLIN=zeros(nchan);
ChanDistanceARC=zeros(nchan);
for a=1:goodHbO%nchan
%     for b=a+1:nchan
%         
%         ChanDistanceLIN(a,b)=sqrt( (zone.pos(a,1)-zone.pos(b,1))^2 + (zone.pos(a,2)-zone.pos(b,2))^2+ (zone.pos(a,3)-zone.pos(b,3))^2);
%         ChanDistanceLIN(b,a)=ChanDistanceLIN(a,b);
%         
%         rayon(a,b)=sqrt( zone.pos(a,1)^2 + zone.pos(a,2)^2 + zone.pos(a,3)^2);
%         rayon(b,a)=rayon(a,b);
%         ChanDistanceARC(a,b)=2*rayon(a,b)*asin(ChanDistanceLIN(a,b)/(2*rayon(a,b)));
%         ChanDistanceARC(b,a)=ChanDistanceARC(a,b);
%     end
    ChanSDseparation(a)=zone.pos(a,4); 
    ChanXYZ(a,:)=zone.pos(a,1:3);
end

ChanDistanceARC=real(ChanDistanceARC);
rayon(rayon==0)=NaN;
mrayon=mean(rayon,'all','omitnan');

ChanNum=1:nchan;


%% calculate proportions of S-D separation
ChanSDseparationINFO={'Source-detector distance' 'Number of channels'};
ChanSDseparationINFO{2,1}={'<1.5 cm'};          ChanSDseparationINFO{2,2}=sum(ChanSDseparation<1.5);
ChanSDseparationINFO{3,1}={'1.5<= x <2.0 cm'};  ChanSDseparationINFO{3,2}=sum([ChanSDseparation<2 & ChanSDseparation>=1.5]);
ChanSDseparationINFO{4,1}={'2.0<= x <2.5 cm'};  ChanSDseparationINFO{4,2}=sum([ChanSDseparation<2.5 & ChanSDseparation>=2]);
ChanSDseparationINFO{5,1}={'2.5<= x <3.0 cm'};  ChanSDseparationINFO{5,2}=sum(ChanSDseparation<3 & ChanSDseparation>=2.5);
ChanSDseparationINFO{6,1}={'3.0<= x <3.5 cm'};  ChanSDseparationINFO{6,2}=sum([ChanSDseparation<3.5 & ChanSDseparation>=3]);
ChanSDseparationINFO{7,1}={'3.5<= x <4.0 cm'};  ChanSDseparationINFO{7,2}=sum(ChanSDseparation<4 & ChanSDseparation>=3.5);
ChanSDseparationINFO{8,1}={'4.0<= x <4.5 cm'};  ChanSDseparationINFO{8,2}=sum([ChanSDseparation<4.5 & ChanSDseparation>=4]);
ChanSDseparationINFO{9,1}={'4.5<= x <5.0 cm'};  ChanSDseparationINFO{9,2}=sum(ChanSDseparation<5 & ChanSDseparation>=4.5);
ChanSDseparationINFO{10,1}={'5.0<= x <5.5 cm'};  ChanSDseparationINFO{10,2}=sum([ChanSDseparation<5.5 & ChanSDseparation>=5]);
ChanSDseparationINFO{11,1}={'>=5.5 cm'};         ChanSDseparationINFO{11,2}=sum(ChanSDseparation>=5.5);

%% save
ChanHelp=sprintf(['ChanNum = Channel ID (#)\n'...
    'ChanSDseparation = length of each channel in centimeters (source-detector separation)\n'...
    'ChanDistanceLIN = matrice of the linear-direct distance between the midpoint of one channel to the midpoint of another channel (in cm)\n'...
    'ChanDistanceARC = matrice of the arc distance between the midpoint of one channel to the midpoint of another channel (in cm). The circle is centered at x=0 y=0 z=0\n'...
    'rayon = In relation to ChanDistanceARC, the radius (rayon) of the circle taken to calculate the distance is stored in <rayon> matrice. It does vary from one pair to the other...\n'...
    'mrayon = average of radius (rayon matrice)\n'...
    'ChanXYZ =  x y z coordinates of each channel (in cm)']);

save([dirmtg 'ChannelsDistance.mat'],'Chan*','rayon','mrayon')

