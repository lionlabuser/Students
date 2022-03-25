savepath='C:\data\Malnutrition\Resting\NIRS\Analyses\CORRmatrice0,01_0,08\Channels\PhysioRegressed_SatRespEKG_corr\CorrPair_criteria\Perm500_results\';

%%Pour matrices%%%
A = gcf;
xlim([0 47]);
ylim([0 47]);
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
ax.FontSize = 14;
titremod = titre(34:end-4);
% grid minor
% ax.Layer = 'top';
% ax.MinorGridLineStyle = '-';
pbaspect([1 1 1])
A.WindowState = 'maximized';
savefig(A,[savepath titremod '.fig'])
exportgraphics(A,[savepath titremod '.png'])
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
exportgraphics(A,[savepath 'PosConnectG1-G2p0.05.png'])
savefig(A,[savepath 'PosConnectG1-G2p0.05.fig'])
clearvars -except savepath
close gcf

%%Pour connectogramme negatif%%%
A = gcf;
ax = gca;
ax.XAxis.Color = 'none';
ax.YAxis.Color = 'none';
pbaspect([1 1 1])
A.WindowState = 'maximized';
exportgraphics(A,[savepath 'NegConnectG1-G2p0.05.png'])
savefig(A,[savepath 'NegConnectG1-G2p0.05.fig'])
clearvars -except savepath
close gcf
