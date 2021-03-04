function [filteredD,globalComp,individualComp]=spatialsmoothingPCA(D,NC,zonefile,Visualize,kernel)
% SCRIPT based on the article from Noah et al. 2021 (https://doi.org/10.1117/1.NPh.8.1.015004)
% & Zhang et al. 2016 (https://doi.org/10.1117/1.NPh.3.1.015004)
%
%
% INPUTS
%   (mandatory)
%      *D - data matrice (time dimension x channel dimension). If the
%       dimensions are inverted, the current script will permute the
%       dimensions to fit the time x channel dimensions.
%
%      *zonefile - string or character array, 
%
%   (optional)
%       Visualize - [1 or 0] to visualize the relationship between the multiplication factor (kernel) and the distance between channels
%
%       kernel - the size of the kernel (in degree or radian). 
%               smoothing kernel set to 46 deg (or 0.8 rad) - Zhang et al 2017 
%               (https://doi.org/10.1117/1.NPh.4.4.041409)
%               as our distances between channels are all in centimeters, we 
%               need to convert the kernel (in degrees) into cm (arc length)
%               kernel needs to be in radians (between 0 and 2pi)
%
% OUTPUTS
%
%
%
if nargin <3
    Visualize=0;
    if nargin <4
        kernel=0.8;
    end
end

%% Calculate channel distance
load( zonefile ,'-mat'); %zone file

nchan=size(zone.pos,1)/2; %number of channels
Dist=zeros(nchan); 
rayon=zeros(nchan);
for a=1:nchan
    for b=a+1:nchan
        rayon(a,b)=sqrt( zone.pos(a,1)^2 + zone.pos(a,2)^2 + zone.pos(a,3)^2);
        rayon(b,a)=rayon(a,b);
        
        tempdistance=sqrt((zone.pos(a,1)-zone.pos(b,1))^2 + (zone.pos(a,2)-zone.pos(b,2))^2+ (zone.pos(a,3)-zone.pos(b,3))^2);
        Dist(a,b)=2*rayon(a,b)*asin(tempdistance/(2*rayon(a,b))); %distance matrice between channels (around a sphere)
        Dist(b,a)=Dist(a,b);
    end
end
Dist=real(Dist);
rayon(rayon==0)=NaN;
mrayon=mean(rayon,'all','omitnan');  %the mean radius that allowed to measure the arc length distance between channels

%% KERNEL SIZE
% smoothing kernel set to 46 deg (or 0.8 rad) - Zhang et al 2017 (https://doi.org/10.1117/1.NPh.4.4.041409)
% as our distances between channels are all in centimeters, we need to
% convert the kernel (in degrees) into cm (arc length)
% kernel needs to be in radians (between 0 and 2pi)
if kernel < 2*pi %kernel in radiam
    ksigma=kernel*mrayon; %now in cm  
else  %kernel in degree
    ksigma=kernel*pi/180*mrayon; %now in cm  
end

%% VERIFY data configuration
sizeD = size(D);
if sizeD(2) == nchan
    %do nothing
elseif sizeD(1) == nchan
    D=permute(D, [2 1]);
else
    error(['The data matrice <D> doesn''t have the same number of channels as the zone file.\n'...
        'Please check your data matrice <D> (size: %d x %d points) and your zone file %s.\n'],sizeD(1),sizeD(2),zonefile);
end

%% PCA DECOMPOSITION
%SpatialSig = spatial matrice. <spatial pattern of global components>
%             each column = a component // each row = each channel
%TemporalSig = temporal matrice. <temporal pattern of global components>
%             each column = a component  // each row = a time point
%ComponentWeigth = weigth of each component <diagonal matrice>
%             diagonal coordinates (square, eg. (1,1) (2,2)) = each component
squareD = D'*D;
[SpatialSig,ComponentWeigth,~]=svd(squareD);
TemporalSig = D*SpatialSig*inv(ComponentWeigth);
%Xm{c}=TemporalSig(:,c)*ComponentWeigth(c,c)*SpatialSig(:,c)';

%% GAUSSIAN SMOOTHING OF THE SPATIAL MATRICE
SmoothSpatialSig=zeros(size(SpatialSig));

for ci=1:size(SpatialSig,1) %for each channel (row)
        for cj=1:size(SpatialSig,1)
            wij(ci,cj)=exp((-(Dist(ci,cj))^2)/(2*ksigma^2)); %calculate the multiplication factor based on the distance between channels ci and cj
            %lorsque ci==cj >> wij = 1
        end
        %wij(ci,:)=wij(ci,:)./sum(wij(ci,:)); %TEST1 NORMALISATION normalized the sum of multiplication factors to 1
end
wij=(wij./sum(wij,'all')).*size(SpatialSig,1); %NORMALISATION CHOSEN. 
% the sum of each column is approx equal to the sum of each column from the
% initial SpatialSig matrice.

for vv=1:size(SpatialSig,2) %for each component (column)
        for ci=1:size(SpatialSig,1) %for each channel (row)
            for cj=1:54
            SmoothSpatialSig(ci,vv)= SmoothSpatialSig(ci,vv) + wij(ci,cj)*SpatialSig(cj,vv);
            if Visualize
                if vv==1 && ci==16 %help to visualize
                    fprintf('Distance: %0.2f -- kernel: %0.10f\n',Dist(ci,cj),wij(ci,cj));
                end
                
            end
        end
        end
    %(SmoothSpatialSig(:,vv)./sum(SmoothSpatialSig(:,vv))).*sum(SpatialSig(:,vv));%TEST2 NORMALISATION 
end
%(SmoothSpatialSig(:,vv)./sum(SmoothSpatialSig(:,vv))).*sum(SpatialSig(:,vv)); %TEST3 NORMALISATION 

globalComp=TemporalSig*ComponentWeigth*SmoothSpatialSig';
filteredD=D-globalComp;

for vv=1:size(SpatialSig,2) 
individualComp(:,:,vv)=TemporalSig(:,vv)*ComponentWeigth(vv,vv)*SmoothSpatialSig(:,vv)';
end

figure
subplot(3,1,1)
plot(D)
title('Original data')
subplot(3,1,2)
plot(globalComp)
title('Global component')
subplot(3,1,3)
plot(D-globalComp)
title('Derived neuronal component')
end