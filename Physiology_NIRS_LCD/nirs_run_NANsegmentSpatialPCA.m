function out = nirs_run_NANsegmentSpatialPCA(job)
% to conduce AFTER yellow identification of bad segments on the whole
% data + nullify bad intervals (yellow periods = NaN)
% the current function: 
% temporary segment the data according to NaN values
% conduce spatial PCA + save in the selected factors file!

NIRS = [];
load(job.NIRSmat{1});
%use last step of operation
lst = length(NIRS.Dt.fir.pp);
rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
NC = NIRS.Cf.H.C.N;
fs = NIRS.Cf.dev.fs;
fprintf('%s\n','File processed');
load(job.physzone{1}   ,'-mat');
ifile = 1; %Utiliser si on remets chaque stim normaliser dans des fichier différent

    [nirsPATH,~,~] = fileparts(job.NIRSmat{1});
    if exist(fullfile(nirsPATH,'SelectedFactors.mat'),'file')
        load(fullfile(nirsPATH,'SelectedFactors.mat'))
        Prow=length(PARCOMP)+1;
    else
        PARCOMP=struct;
        Prow=1;
    end
    
   if isfield(job,'goodPercent')
       goodPercent=job.goodPercent;
   else
       goodPercent=.5; % if less than 50% of channels are marked 
   end                 % as good, it will segment it as a "block"
                       % to apply spatial PCA

for f=1:numel(rDtp) %Loop over all files of a NIRS.mat
    d = fopen_NIR(rDtp{f},NC);
    
    if size(d,1)==NC
        %do nothing
    elseif size(d,2)==NC
        d=permute(d,[2 1]); %organize the data according to Dim1=channels and Dim2 = time
    else
        error ('The size of your data does not match the number of channels...')
    end
    
    nanid=~isnan(d); %identifies with a 1 all data positions that ARE NOT a NAN value
    nansum=sum(nanid)>(NC*goodPercent); %for each time point, 
                     % check if you have more valid channels than 
                     % the cut-off percentage. If so, return a 1 value (if
                     % not, return a 0 value)
    nansum=[false nansum false];
    nandiff=diff(nansum);
    blockend=find(nandiff==-1)-1;
    blockstart=find(nandiff==1);
    
    if ~(length(blockstart)==length(blockend))
        error('Not the same length of ID position of blocks onsets and ends')
    end
    
    newdata =nan(size(d));
    Xm= nan(size(d));
    for bl=1:length(blockstart)
        tempdata=d(:,blockstart(bl):blockend(bl))';
        temphbo=tempdata(:,1:NC/2);
        CHhbo=1:NC/2;
        temphbr=tempdata(:,1+NC/2:NC);
        
        %check for bad channels
        goodCHhbo=CHhbo(~(sum(isnan(temphbo))));
        if goodCHhbo==CHhbo
            %do nothing
        else
            temphbo=temphbo(:,goodCHhbo);
            temphbr=temphbr(:,goodCHhbo);
        end
        [filteredD,globalComp,~]=spatialPCA(temphbo,goodCHhbo,zone,0,0.8);
        [filteredD2,globalComp2,~]=spatialPCA(temphbr,goodCHhbo,zone,0,0.8); %as hbo and
        %hbr have the same spatial position, easier to take the good hbo list instead of the hbr list
        
        
        newdata(goodCHhbo,blockstart(bl):blockend(bl)) =  filteredD';
        newdata(goodCHhbo+NC/2,blockstart(bl):blockend(bl)) =  filteredD2';
        
        Xm(goodCHhbo,blockstart(bl):blockend(bl))=globalComp';
        Xm(goodCHhbo+NC/2,blockstart(bl):blockend(bl))=globalComp2';
        
    end
    %% WRITE PARCOMP and BETAs
    AUX.data{1}=mean(Xm(1:NC/2,:),'omitnan');
    AUX.data{2}=mean(Xm(1+NC/2:NC,:),'omitnan');
    PARCOMP(Prow).file= f;
    PARCOMP(Prow).filestr =  sprintf('Bl%02.0f', f);
    PARCOMP(Prow).label= ['SpatialPCA_' PARCOMP(Prow).filestr];
    PARCOMP(Prow).type = 'GLM';
    % PARCOMP(Prow).beta = tmpbeta;
    % PARCOMP(Prow).std= 0;
    PARCOMP(Prow).AUX.data= AUX.data;
    % PARCOMP(Prow).AUX.IndividualComp={individualComp, individualComp2};
    PARCOMP(Prow).AUX.label= {'SpatialPCA HbO' 'SpatialPCA HbR' };
    PARCOMP(Prow).data = d ;
    PARCOMP(Prow).Xm = Xm;
    PARCOMP(Prow).dataCORR = newdata;
    % PARCOMP(Prow).indt = [tstart :tstop]; %indice de temps.
    PARCOMP(Prow).listgood =  [1:NC];
    PARCOMP(Prow).module  = lst;
    PARCOMP(Prow).modulestr = NIRS.Dt.fir.pp(lst).pre;
    %PARCOMP(Prow).ComponentToKeep = 1;
    %PARCOMP(Prow).idreg = 1;
    %PARCOMP(Prow).topo =  tmpbeta(PARCOMP(Prow).ComponentToKeep,:);
    Prow=Prow+1;
    clear tmp*
    
end
save(fullfile(nirsPATH,'SelectedFactors.mat'),'PARCOMP');
out.NIRSmat = job.NIRSmat;

end
function [filteredD,globalComp,individualComp]=spatialPCA(D,chlist,z,Visualize,kernel)
% SCRIPT based on the article from Noah et al. 2021 (https://doi.org/10.1117/1.NPh.8.1.015004)
% & Zhang et al. 2016 (https://doi.org/10.1117/1.NPh.3.1.015004)
%
if nargin <3
    Visualize=0;
    if nargin <4
        kernel=0.8;
    end
end

%% Calculate channel distance

nchan=length(chlist); %number of channels
Dist=zeros(nchan);
rayon=zeros(nchan);
for a=1:nchan
    for b=a+1:nchan
        rayon(a,b)=sqrt( z.pos(chlist(a),1)^2 + z.pos(chlist(a),2)^2 + z.pos(chlist(a),3)^2);
        rayon(b,a)=rayon(a,b);
        
        tempdistance=sqrt((z.pos(chlist(a),1)-z.pos(chlist(b),1))^2 + (z.pos(chlist(a),2)-z.pos(chlist(b),2))^2+ (z.pos(chlist(a),3)-z.pos(chlist(b),3))^2);
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
        'Please check your data matrice <D> (size: %d x %d points) and your zone file.\n'],sizeD(1),sizeD(2));
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

%% GAUSSIAN SMOOTHING OF THE SPATIAL MATRICE
SmoothSpatialSig=zeros(size(SpatialSig));

for ci=1:size(SpatialSig,1) %for each channel (row)
    for cj=1:size(SpatialSig,1)
        wij(ci,cj)=exp((-(Dist(ci,cj))^2)/(2*ksigma^2)); %calculate the multiplication factor based on the distance between channels ci and cj
        %lorsque ci==cj >> wij = 1
    end
    %wij(ci,:)=wij(ci,:)./sum(wij(ci,:)); %TEST1
    %   NORMALISATION normalized the sum of multiplication factors to
    %   1. FOR EACH CHANNEL, the sum of its multiplication factors
    %   equals 1. For channels that have a lot of neighbors, it makes
    %   that the weigth of the current channel is smaller compared to
    %   the total weigth. For channels that have fewer neighbors,
    %   their own weigth is larger in proportion to the total.
    
end
wij=(wij./sum(wij,'all')).*size(SpatialSig,1); %NORMALISATION CHOSEN.
% the sum of each column is approx equal to the sum of each column from the
% initial SpatialSig matrice.
% It takes the WHOLE MATRICE to do the normalization. therefore the sum of
% weigths for each channel is not necessarily 1 (the mean sum across channels = 1)

for vv=1:size(SpatialSig,2) %for each component (column)
    for ci=1:size(SpatialSig,1) %for each channel (row)
        for cj=1:size(SpatialSig,1)
            SmoothSpatialSig(ci,vv)= SmoothSpatialSig(ci,vv) + wij(ci,cj)*SpatialSig(cj,vv);
            if Visualize
                if vv==1 && ci==16 %help to visualize
                    fprintf('Distance: %0.2f -- kernel: %0.10f\n',Dist(ci,cj),wij(ci,cj));
                end
                
            end
        end
    end
end
globalComp=TemporalSig*ComponentWeigth*SmoothSpatialSig';
filteredD=D-globalComp;

for vv=1:size(SpatialSig,2)
    individualComp(:,:,vv)=TemporalSig(:,vv)*ComponentWeigth(vv,vv)*SmoothSpatialSig(:,vv)';
end
if Visualize
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
end