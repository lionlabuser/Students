savepath='C:\data\Malnutrition\Resting\NIRS\Analyses préliminaires\CORRmatrice0,01_0,08\Channels\Perm_results N=54\';

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
ax.FontSize = 12;
titremod = titre(34:end-4);
% grid minor
% ax.Layer = 'top';
% ax.MinorGridLineStyle = '-';
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
exportgraphics(A,[savepath 'PosConnectG2-G1p0.05.png'])
savefig(A,[savepath 'PosConnectG2-G1p0.05.fig'])
clearvars -except savepath
close gcf

%%Pour connectogramme negatif%%%
A = gcf;
ax = gca;
ax.XAxis.Color = 'none';
ax.YAxis.Color = 'none';
exportgraphics(A,[savepath 'NegConnectG2-G1p0.05.png'])
savefig(A,[savepath 'NegConnectG2-G1p0.05.fig'])
clearvars -except savepath
close gcf
