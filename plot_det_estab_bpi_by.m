% plot_det_estab_bpi_by.m -- vector figure for det_estab_bpi_by.m (bpi on the
% horizontal axis, by on the vertical axis), same style as plot_det_estab_bpi_phib.m.
S = load('det_estab_bpi_by_results.mat');
cm = S.code_map; bg = S.bpi_grid; yg = S.by_grid;
dx = bg(2)-bg(1); dy = yg(2)-yg(1);
col.det = [0.30 0.70 0.40]; col.ind = [0.85 0.30 0.30];

fig = figure('Units','inches','Position',[0.5 0.5 7 5.5],'Color','w');
ax = axes('Position',[0.11 0.22 0.85 0.74]); hold(ax,'on');
patch(ax,[0 5 5 0],[0 0 3 3],col.ind,'EdgeColor','none');
X = []; Y = [];
for i = 1:numel(yg)
    y0 = max(yg(i)-dy/2,0); y1 = min(yg(i)+dy/2,3);
    j = 1;
    while j <= numel(bg)
        if cm(i,j) == 0
            k = j; while k < numel(bg) && cm(i,k+1) == 0, k = k+1; end
            x0 = max(bg(j)-dx/2,0); x1 = min(bg(k)+dx/2,5);
            X = [X,[x0;x1;x1;x0]]; Y = [Y,[y0;y0;y1;y1]]; %#ok<AGROW>
            j = k+1;
        else
            j = j+1;
        end
    end
end
rects_union_plot(ax,X,Y,col.det);

plot(ax,S.bpi_est,S.by_est,'k*','MarkerSize',14,'LineWidth',2.5);
xline(ax,1,'k--','LineWidth',1.2);
xlim(ax,[0 5]); ylim(ax,[0 3]); box(ax,'on'); set(ax,'Layer','top','TickDir','out');
xlabel(ax,'$b_\pi$','Interpreter','latex','FontSize',15);
ylabel(ax,'$b_y$','Interpreter','latex','FontSize',15);
text(ax,S.bpi_est+0.12,S.by_est+0.12, ...
    sprintf('Posterior mean $(%.2f,\\, %.3f)$',S.bpi_est,S.by_est), ...
    'FontSize',9,'Interpreter','latex','VerticalAlignment','bottom');
text(ax,1.06,2.7,'$b_\pi=1$','FontSize',10,'Interpreter','latex');

h = gobjects(1,2);
h(1) = patch(ax,NaN,NaN,col.det,'EdgeColor','none'); h(2) = patch(ax,NaN,NaN,col.ind,'EdgeColor','none');
lg = legend(ax,h,{'Determinate & E-stable','Indeterminate (E-unstable)'},'Location','southoutside', ...
    'Orientation','horizontal','Box','off');
lg.FontSize = 10;

exportgraphics(fig,fullfile(pwd,'det_estab_bpi_by_rnfixed.pdf'),'ContentType','vector','BackgroundColor','white');
close(fig);
fprintf('Saved det_estab_bpi_by_rnfixed.pdf\n');
