



Data=v;
Dist=ChanDistanceBTW;
for vv=1:54
for i=1:54
    for j=1:54
        if i==j
            continue
        end
        wij(i,j)=exp((-(Dist(i,j))^2)/(2*50^2)); %0.8 rad = 46 degre
    end
    wij(i,:)=wij(i,:)./sum(wij(i,:));
    
    for j=1:54
        if i==j
            continue
        end
        beforesum(j)=wij(i,j)*Data(j,vv);
         if vv==1 && i==16
            fprintf('Distance: %0.2f -- kernel: %0.10f\n',Dist(i,j),wij(i,j));
        end
    end
        GlobalComp(i,vv)=sum(beforesum);
        clear beforesum
end
    
end


GG=u*s*GlobalComp';
GG1=u(:,1)*s(1,1)*GlobalComp(:,c)';
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