% SCRIPT based on the article from Noah et al. 2021 (https://doi.org/10.1117/1.NPh.8.1.015004)
% & Zhang et al. 2016 (https://doi.org/10.1117/1.NPh.3.1.015004)

%% PARAMETERS
Visualize=1; %[1 or 0] to visualize the relationship between the multiplication factor (kernel) and the distance between channels
OriginalSpatial=v; % spatial matrice. <spatial pattern of global components>
Dist=ChanDistanceARC; %distance matrice between channels (around a sphere)
rayon=mrayon; %the mean radius that allowed to measure the arc length distance between channels

%% KERNEL SIZE
% smoothing kernel set to 46 deg (or 0.8 rad) - Zhang et al 2017 (https://doi.org/10.1117/1.NPh.4.4.041409)
% as our distances between channels are all in centimeters, we need to
% convert the kernel (in degrees) into cm (arc length)
kernel=0.8*rayon; %now in cm

%% GAUSSIAN SMOOTHING OF THE SPATIAL MATRICE
SmoothSpatial=zeros(size(OriginalSpatial));

for vv=1:size(OriginalSpatial,2) %for each component (column in the Data matrice)
    for ci=1:size(OriginalSpatial,1) %for each channel (row in the Data matrice)
        for cj=1:size(OriginalSpatial,1)
            wij(ci,cj)=exp((-(Dist(ci,cj))^2)/(2*kernel^2)); %calculate the multiplication factor based on the distance between channels ci and cj
            %lorsque ci==cj >> wij = 1
        end
        wij(ci,:)=wij(ci,:)./sum(wij(ci,:)); %normalized the sum of multiplication factors to 1
        
        for cj=1:54
            SmoothSpatial(ci,vv)= SmoothSpatial(ci,vv) + wij(ci,cj)*OriginalSpatial(cj,vv);
            if Visualize
                if vv==1 && ci==16 %help to visualize
                    fprintf('Distance: %0.2f -- kernel: %0.10f\n',Dist(ci,cj),wij(ci,cj));
                end
                disp('===================');
            end
        end
    end
end


GG=u*s*SmoothSpatial';
GG1=u(:,1)*s(1,1)*SmoothSpatial(:,c)';
figure
subplot(3,1,1)
plot(data)
title('Original data')
subplot(3,1,2)
plot(GG)
title('Global component')
subplot(3,1,3)
plot(data-GG)
title('Derived neuronal component')