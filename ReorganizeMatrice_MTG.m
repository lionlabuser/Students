
dirmtg='C:\Data\data_NIRS\ELAN\ANALYSED\Martine_0m\Montage\';
load( [dirmtg 'Global.zone'] ,'-mat'); %zone file
mtg2d = readtable([dirmtg '2D_montage_Martine.xlsx'],'ReadVariableNames',0); % Import the data
mtg2d = table2array(mtg2d); %convert to matrix



%% help info zone file
%zone.pos(channel#,[x pos, y pos, z pos, channel distance in CM])
%zone.ml (channel#, [Source#, Detector#, ?, HbO(1) or HbR(2)]

%% get position on the spatial 2D matrice for each channel
nchan=size(zone.pos,1)/2; %number of channels
for a=1:nchan
    try
        [mtgXY(a,1),mtgXY(a,2)]= find(mtg2d==a);
    catch
        disp(['Problem with Channel # ' num2str(a) ': error'])
    end
end
save([dirmtg 'mtg2d.mat'],'mtg*')

%% ACTIVATION -MATRICE - reorganize according to the spatial matrice
load('C:\Data\data_NIRS\ELAN\ANALYSED\Martine_0m\001-050Group\MatriceBB\All\B16_all.mat','A') %dataset

MAT=A(:,:,1);  %1st dimension = channels
%2nd dimension = time points
%3rd dimension = blocks -- suppress this dimension for the purpose
%of the script. Here for example keeps only the first block


for t=1:size(MAT,2)
    tempMAT=NaN(size(mtg2d));
    for a=1:nchan
        tempMAT(mtgXY(a,1),mtgXY(a,2))=MAT(a,t);
    end
    
    gtempMAT = imgaussfilt(tempMAT,5,'FilterDomain','spatial');
    isc=imagesc(gtempMAT,'AlphaData',~isnan(tempMAT));
    
    cmap=colormap(jet); %blue to red
    cmap=colormap(parula); %blue to yellow
    isc=imagesc(tempMAT,'AlphaData',~isnan(tempMAT)); %color plot matrix
    
    c=colorbar; %legend color bar
        c.FontSize=10;
end


%% Because XYZ position are in centimeters = need to be converted in <pixels> (real integers)
% % %To create a 4d matrice (Xpos, Ypos, Zpos, Time)
% %
% % center =[    0.0095   -0.0002    0.0379]; %in cm; taken from a script in th iomtg... might be useful for a sphere??
% %
% % %will need to convert XYZ pos into real integer to be able to be read
% %
% % Multiplier=3; %choisi un peu arbitrairement, en fonction de la distance actuelle minimale entre les points
% % Pos=round(zone.pos(1:nchan,1:3)*Multiplier);
% % tempPos=sort(Pos);
% % tempDis=tempPos(1:end-1,:)-tempPos(2:end,:);
% % minInt=min(min(abs(tempDis)));
% % minInt=1;
% % % while minInt==0
% % %    Multiplier=Multiplier+50;
% % %    Pos=round(zone.pos(1:nchan,1:3)*Multiplier);
% % % tempPos=sort(Pos);
% % % tempDis=tempPos(1:end-1,:)-tempPos(2:end,:);
% % % minInt=min(min(abs(tempDis)));
% % %
% % % end
% % axeX=min(tempPos(:,1)) : minInt:    max(tempPos(:,1));
% % axeY=min(tempPos(:,2)) : minInt:    max(tempPos(:,2));
% % axeZ=min(tempPos(:,3)) : minInt:    max(tempPos(:,3));
% %
% % for a=1:nchan
% %       [~,realPos(a,1)]=min(abs(axeX-Pos(a,1)));
% %       [~,realPos(a,2)]=min(abs(axeY-Pos(a,2)));
% %       [~,realPos(a,3)]=min(abs(axeZ-Pos(a,3)));
% % end
% %
% %
% % %% manual check for channels that have a similar position (eg. croisement X)
% % pairs=[4 11; 8 15; 9 19; 25 26; 29 30; 36 44 ; 39 46 ; 50 52]; %adjust according to your own montage
% %   %les paires de canaux qui se croisent... et donc que leur position
% %   %pourrait etre la meme en coordonnees 3D
% % for p=1:size(pairs,1)
% %
% %
% % end
% % question=input('Have you check if some channels had the same 3D coordinates?? [1=yes, 0=no]\n');
% %
% % if question=1
% % mat4D=NaN([length(axeX) length(axeY) length(axeZ) size(MAT,2)]);
% %
% % for tt=1:size(MAT,2)
% %
% %     for a=1:nchan
% %
% %         mat4D(realPos(a,1),realPos(a,2),realPos(a,2),1)=MAT(a,t);
% %     end
% % end
% %
% % end
