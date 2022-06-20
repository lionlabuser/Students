savepath='C:\data\Malnutrition\Resting\NIRS\Analyses\CORRmatrice0,01_0,08\PCAPW\CorrPairC0,25 ExcY\Perm_results500\';

%%Pour matrices%%%
A = gcf;
xlim([0 46]);
ylim([0 46]);
xtickangle(90);

x = xticklabels;
for i = 1:numel(x)
    label = x{i,1};
    modlabels{i,1} = strrep(label(1:end-4),'_',' ');
end
xticklabels(modlabels);
yticklabels(modlabels);
ax = gca;
titre = ax.Title.String;
ax.FontSize = 20;
titremod = titre(34:end-4);
% grid minor
% ax.Layer = 'top';
% ax.MinorGridLineStyle = '-';
pbaspect([1 1 1])
A.WindowState = 'maximized';
savefig(A,[savepath  titremod '.fig']) %'z'
exportgraphics(A,[savepath titremod '.png']) %'z' 
disp(titremod)
clearvars -except savepath
close gcf

%%Pour connectogramme positif%%%
A = gcf;
ax = gca;
ax.XAxis.Color = 'none';
ax.YAxis.Color = 'none';
pbaspect([1 1 1])
A.WindowState = 'maximized';
exportgraphics(A,[savepath 'PosConnectG1-G2p0.05.png']) %ROI
savefig(A,[savepath 'PosConnectG1-G2p0.05.fig']) %ROI
clearvars -except savepath
close gcf

%%Pour connectogramme negatif%%%
A = gcf;
ax = gca;
ax.XAxis.Color = 'none';
ax.YAxis.Color = 'none';
pbaspect([1 1 1])
A.WindowState = 'maximized';
exportgraphics(A,[savepath 'NegConnectG1-G2p0.05.png']) %ROI
savefig(A,[savepath 'NegConnectG1-G2p0.05.fig']) %ROI
clearvars -except savepath
close gcf
