function [FSupSup,FSupDeriv,FSupTime,FUniv,T0ij] = Permutations_onesample(ESTAD,INDEP,data,NPERM)
%  ESTAD - Statistic to use
%         1 - t-student with "original" values
%         2 - sum of difference with "original" values
%         4 - t-student with absolute values
%         5 - sum of difference with absolute values
% 
%  INDEP- Indice if groups are independent 
%         0 - DEPENDENT 
%         1 - INDEPENDENT
%
%  data - data matrix : 3D --> (subjects x ROI/channel x time point/interval) 
%
%  NPERM  - Number of permutations to conduct
%%OUTPUT
%   FSupSup : global significant probability that the mean is equal to zero 
%   FSupDeriv: vector (1:numROI) of significant probabilities that the mean is equal to zero 
%       for each ROI/channel with respect of all the time points (multivariate)
%   Fsuptime: vector (1:numtime) of significant probabilities that the mean is equal to zero 
%       for each time point/interval with respect of all ROI/channels (multivariate)
%   Funiv : matriz (numROI x numtime) of significant probabilities that the mean is equal 
%       to zero for each univariate variable (for each combination of ROIxtime)  

[nx1,~,~] = size(data);
%[nx2,ny2,nz2] = size(DatGrup2);

if ESTAD == 4||5
    ValAbs = 1;
else
    ValAbs = 0;
end

if ESTAD == 1||4
         [T0ij,T0_deriv,T0_time,T0_global] = TStudent_oneS(data,ValAbs);
     else
         [T0ij,T0_deriv,T0_time,T0_global] = SumDiff_oneS(data,ValAbs);            
end

[FSupSup,FSupDeriv,FSupTime,FUniv] = Permut_oneS(ESTAD,INDEP,data,T0_global,T0_deriv,T0_time,T0ij,NPERM,ValAbs);
FSupSup = 1 - FSupSup;
FSupDeriv = 1 - FSupDeriv ;
FSupTime  = 1 - FSupTime ;
FUniv= 1- FUniv;
% 
% 

function [FSupSup,FSupDeriv,FSupTime,FUniv] = Permut_oneS(ESTAD,INDEP,Dat,n1,n2,T0_global,T0_deriv,T0_time,T0ij,NPERM,ValAbs)

%global VG1 VG2 
%global TP

[nx1,~,~] = size(data);

v1 = (1:n1)';

FSupSup = 0;
dimderiv = length(T0_deriv);
dimtime  = length(T0_time);
FSupDeriv = zeros(dimderiv,1);
FSupTime = zeros(dimtime,1);
   FUniv = zeros(size(T0ij));
V1 = [];
h = waitbar(0,'Computing Permutations...');
for iperm = 1:NPERM
   %%Calculo del vector de permutaciones
      [vg1,vg2] = GenVecPermDep(nx1);
   V1 = [V1,vg1];
   
   switch(ESTAD)
   case {1,4}
         tmpdat = Dat(vg1,:,:) - Dat(vg2,:,:);
         [TP,TP_deriv,TP_time,TP_global] = TStudent_oneS(tmpdat,ValAbs);
      
      [d,dim_deriv,dim_time] = size(Dat);
      d1 = dim_deriv; d2 = dim_time;
      %Distribucion del SUPSUP
      FSupSup = FSupSup + (TP_global < T0_global);
      %Distribucion del SUPj
      FSupDeriv = FSupDeriv + (TP_deriv < T0_deriv);
      %Distribucion del SUPi
      FSupTime = FSupTime + (TP_time < T0_time);
       FUniv = FUniv + ( TP < T0ij);
      
   case {2,5}
     
      [TP,TP_deriv,TP_time,TP_global] = SumDiff_oneS(Dat,vg1,vg2,ValAbs);
      
      [d,dim_deriv,dim_time] = size(Dat);
      d1 = dim_deriv; d2 = dim_time;
      
      %Distribucion del SUPSUP
      FSupSup = FSupSup + (TP_global < T0_global);
      %Distribucion del SUPj
      FSupDeriv = FSupDeriv + (TP_deriv < T0_deriv);
      %Distribucion del SUPi
      FSupTime = FSupTime + (TP_time < T0_time);
      FUniv = FUniv + ( TP < T0ij);

end

waitbar(iperm/NPERM,h)      
end
close(h)

%Correción propuesta por Pesarin   (puede estar con o sin corrección)
FSupSup  = (FSupSup + 0.5 )./(NPERM + 1);
FSupDeriv  = (FSupDeriv + 0.5 )./(NPERM + 1);
FSupTime  = (FSupTime + 0.5 )./(NPERM + 1);
FUniv  = (FUniv+ 0.5 )./(NPERM+ 1 );

function [vp1,vp2] = GenVecPermDep(N)

perm = randperm(N*2);
vp1 = perm(1:N);

for i = 1:N
    if vp1(i) <= N 
        vp2(i)=vp1(i)+N;
    else
        vp2(i)=vp1(i)-N;
    end
end
vp1=vp1';
vp2=vp2';




function [Tij,T_deriv,T_time,T_global] = SumDiff_oneS(data,ValAbs)

Tij = squeeze(sum(data,1))';
if ValAbs == 1
   Tij = abs(Tij);
end
Tij = Tij';

[~,dim_deriv,dim_time] = size(data);
if dim_deriv == 1
   T_deriv = Tij;
else
   T_deriv = max(Tij')';
end

if dim_time == 1 
   T_time = Tij';
else
   T_time = (max(Tij))';
end

T_global   = max([T_time ; T_deriv]);

function [Tij,T_deriv,T_time,T_global] = TStudent_oneS(data,ValAbs)
   

[nx1,~,~] = size(data);
[M,S] = MeanSigmaMat_corrected(data);
Tij = M./sqrt(S/nx1);
if (ValAbs == 1)
   Tij = abs(Tij);
end
Tij = Tij';

[~,dim_deriv,dim_time] = size(data);
if dim_deriv == 1
   T_deriv = Tij;
else
   T_deriv = max(Tij')';
end

if dim_time == 1   
   T_time = Tij';
else
   T_time = (max(Tij))';
end

T_global   = max([T_time ; T_deriv]);

function [MX,SX] = MeanSigmaMat_corrected(X,v)


MX = mean(X(v,:,:),1,'omitnan');			%%Se suma a travez de la primera dimension
MX = squeeze(MX)';

SX = std(X(v,:,:),0,1,'omitnan').^2;    %%El 0 indica que se normaliza respecto a N-1
SX = squeeze(SX)';


      
