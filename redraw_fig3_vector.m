% redraw_fig3_vector.m
% Redraws Figure 3 (determinacy / E-stability region) as pure vector shapes.
% The original (det_estab_mmt8bOccBin3_covid_bT_baa.m, imagesc + colorbar,
% exportgraphics 'vector') stores the heatmap as a JPEG wrapped in a tiling
% pattern fill, which Acrobat can render blank. The 60x60 classification grid
% was recovered by sampling cell centres of that JPEG (fig3_combined_grid.mat;
% 2857 green cells, matching the 2,857 reported in the text).
S = load('fig3_combined_grid.mat');
comb = S.comb; bpi_grid = S.bpi_grid; by_grid = S.by_grid;
N = numel(bpi_grid);
dby = by_grid(2)-by_grid(1); dbp = bpi_grid(2)-bpi_grid(1);
red = [0.85 0.30 0.30]; grn = [0.30 0.70 0.40];
bpi_est = 3.474; by_est = 0.083;

fig = figure('Units','inches','Position',[0.5 0.5 7 5.5],'Color','w');
ax = axes('Position',[0.11 0.22 0.85 0.74]); hold(ax,'on');
patch(ax,[0 3 3 0],[0 0 5 5],red,'EdgeColor','none');

X = []; Y = [];
for i = 1:N
    y0 = max(bpi_grid(i)-dbp/2,0); y1 = min(bpi_grid(i)+dbp/2,5);
    j = 1;
    while j <= N
        if comb(i,j) == 1
            k = j; while k < N && comb(i,k+1) == 1, k = k+1; end
            x0 = max(by_grid(j)-dby/2,0); x1 = min(by_grid(k)+dby/2,3);
            X = [X, [x0;x1;x1;x0]]; Y = [Y, [y0;y0;y1;y1]]; %#ok<AGROW>
            j = k+1;
        else
            j = j+1;
        end
    end
end
patch(ax,X,Y,grn,'EdgeColor',grn,'LineWidth',0.3);

plot(ax,by_est,bpi_est,'k*','MarkerSize',14,'LineWidth',2.5);
yline(ax,1,'k--','LineWidth',1.4);
xlim(ax,[0 3]); ylim(ax,[0 5]); box(ax,'on'); set(ax,'Layer','top','TickDir','out');
xlabel(ax,'$b_y$','Interpreter','latex','FontSize',15);
ylabel(ax,'$b_\pi$','Interpreter','latex','FontSize',15);
text(ax,by_est+0.08,bpi_est,sprintf('Posterior mean $(%.2f,\\, %.2f)$',bpi_est,by_est), ...
    'FontSize',9,'VerticalAlignment','middle','HorizontalAlignment','left','Interpreter','latex');

hg = patch(ax,NaN,NaN,grn,'EdgeColor','none'); hr = patch(ax,NaN,NaN,red,'EdgeColor','none');
lg = legend(ax,[hg hr],{'Determinate & E-stable','Indeterminate/explosive & E-unstable'}, ...
    'Location','southoutside','Orientation','horizontal','Box','off');
lg.FontSize = 10;

outfile = fullfile(pwd,'det_estab_map_mmt8bOccBin3_covid_bT_baa.pdf');
exportgraphics(fig,outfile,'ContentType','vector','BackgroundColor','white');
fprintf('Saved %s\n',outfile);
close(fig);
