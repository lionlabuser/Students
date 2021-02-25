function [ rebuild1 ax,ay]=removeglobalY(data,polhemus,sigra , badch);
% Global component remove program is made by Dr. Zhang, Xian,  Brain function lab, Yale school of medicine, Psychiatry department
% Please site my paper: "Separation of the global and local components in functional near-infrared spectroscopy signals using principal component spatial filtering"
% and "Signal processing of functional NIRS data acquired during overt speaking" which indicating global component in dexoyHb signal
% xzhang63@gmail.com
nTR=polhemus.nTR;
sz=size(data);
nch0=size(data,3);
nch1=length(polhemus.result4);
nch=min( [ nch0 nch1]);
data(isnan(data))=0;
ind1=1:nch;
badch(badch>nch)=[];
gr=GlobalRemover(polhemus,sigra,badch);
nch=nch-length(badch);
ind1(badch)=[];
for c=1:nch
    temp=data(:,:,ind1(c),:);
    temp=temp(:);
    data1(:,c)=temp;
end
[S,V,D]=svd(data1,'econ');
nPC=size(D,2);
%D1=D;
for i=1:nPC
 D1(:,i) = gr.getGlobal(D(:,i));
end
data2=S*V*D1';
data1=data1-data2;
if nch<nch0
    data1(:,[(nch+1):nch0])=0;
end
for c=1:nch
    temp=data1(:,c);
    switch length(sz)
        case 4
            rebuild(:,:,c,:)=reshape(temp,sz(1),sz(2),1,sz(4));
        case 5
            rebuild(:,:,c,:,:)=reshape(temp,sz(1),sz(2),1,sz(4),sz(5));
        case 3
            rebuild(:,:,c)=reshape(temp,sz(1),sz(2)); 
    end
end
rebuild1=data*0;
switch length(sz)
case 4
    rebuild1(:,:,ind1,:)=rebuild;
case 5 % multiple condition
    rebuild1(:,:,ind1,:,:)=rebuild;
case 3
    rebuild1(:,:,ind1)=rebuild;
end
