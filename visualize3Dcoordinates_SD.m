%%visualise 3D coordinates from nirs.mat



%NEED TO LOAD FIRST A NIRS.MAT FILE

%number of sources and detectors
ndetector=NIRS.Cf.H.D.N; 
nsource=NIRS.Cf.H.S.N;

%3D coordinates of sources and detectors
coorD=NIRS.Cf.H.D.r.o.mm.p(:,1:ndetector);
coorS=NIRS.Cf.H.S.r.o.mm.p(:,1:nsource);


%3D plot to help visualize
%blue for detectors / red for sources
% + for odd numbers (impairs)
% o for even numbers (pairs)
figure
grid on
plot3(coorS(1,:),coorS(2,:),coorS(3,:),'color','r');hold on
plot3(coorD(1,:),coorD(2,:),coorD(3,:),'color','b');hold on

plot3(coorS(1,[1:2:nsource]),coorS(2,[1:2:nsource]),coorS(3,[1:2:nsource]),'+','color','r','LineWidth',2);hold on
plot3(coorD(1,[1:2:ndetector]),coorD(2,[1:2:ndetector]),coorD(3,[1:2:ndetector]),'+','color','b','LineWidth',2);hold on

plot3(coorS(1,[2:2:nsource]),coorS(2,[2:2:nsource]),coorS(3,[2:2:nsource]),'o','color','r','LineWidth',2);hold on
plot3(coorD(1,[2:2:ndetector]),coorD(2,[2:2:ndetector]),coorD(3,[2:2:ndetector]),'o','color','b','LineWidth',2);hold on
