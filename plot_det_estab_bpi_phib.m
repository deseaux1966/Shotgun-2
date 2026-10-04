% plot_det_estab_bpi_phib.m
% Vector figure (no images/patterns, so Acrobat renders it reliably) of the
% determinacy / E-stability map from det_estab_bpi_phib.m: bpi in [0,5],
% phib_val in [0,1], other parameters (incl. by) at the posterior means of
% the sigma_me,rn-fixed production round.
S = load('det_estab_bpi_phib_results.mat');
cm = S.code_map; bg = S.bpi_grid; pg = S.phib_grid;
dx = bg(2)-bg(1); dy = pg(2)-pg(1);
col.det  = [0.30 0.70 0.40];   % determinate & E-stable
col.ind  = [0.85 0.30 0.30];   % indeterminate (BK info 4)
col.expl = [0.93 0.60 0.20];   % explosive / no stable equilibrium (BK info 3)
col.rank = [0.55 0.55 0.55];   % rank condition fails (info 5)

fig = figure('Units','inches','Position',[0.5 0.5 7 5.5],'Color','w');
ax = axes('Position',[0.11 0.22 0.85 0.74]); hold(ax,'on');
patch(ax,[0 5 5 0],[0 0 1 1],col.ind,'EdgeColor','none');   % base: indeterminate

cats = {0,'det'; 3,'expl'; 5,'rank'};
for c = 1:size(cats,1)
    X = []; Y = [];
    for i = 1:numel(pg)
        y0 = max(pg(i)-dy/2,0); y1 = min(pg(i)+dy/2,1);
        j = 1;
        while j <= numel(bg)
            if cm(i,j) == cats{c,1}
                k = j; while k < numel(bg) && cm(i,k+1) == cats{c,1}, k = k+1; end
                x0 = max(bg(j)-dx/2,0); x1 = min(bg(k)+dx/2,5);
                X = [X,[x0;x1;x1;x0]]; Y = [Y,[y0;y0;y1;y1]]; %#ok<AGROW>
                j = k+1;
            else
                j = j+1;
            end
        end
    end
    if ~isempty(X)
        rects_union_plot(ax,X,Y,col.(cats{c,2}));
    end
end

plot(ax,S.bpi_est,S.phib_est,'k*','MarkerSize',14,'LineWidth',2.5);
xline(ax,1.09,'k--','LineWidth',1.2);
xlim(ax,[0 5]); ylim(ax,[0 1]); box(ax,'on'); set(ax,'Layer','top','TickDir','out');
xlabel(ax,'$b_\pi$','Interpreter','latex','FontSize',15);
ylabel(ax,'$\varphi_b$','Interpreter','latex','FontSize',15);
text(ax,S.bpi_est+0.12,S.phib_est+0.035, ...
    sprintf('Posterior mean $(%.2f,\\, %.3f)$',S.bpi_est,S.phib_est), ...
    'FontSize',9,'Interpreter','latex','VerticalAlignment','bottom');
text(ax,1.12,0.55,'$b_\pi\approx1.09$','FontSize',10,'Interpreter','latex');

h = gobjects(1,4);
names = {'Determinate & E-stable','Indeterminate','Explosive (no stable equilibrium)','Rank condition fails'};
fcs = {col.det,col.ind,col.expl,col.rank};
for q = 1:4, h(q) = patch(ax,NaN,NaN,fcs{q},'EdgeColor','none'); end
lg = legend(ax,h,names,'Location','southoutside','Orientation','horizontal','NumColumns',2,'Box','off');
lg.FontSize = 9;

exportgraphics(fig,fullfile(pwd,'det_estab_bpi_phib_rnfixed.pdf'),'ContentType','vector','BackgroundColor','white');
close(fig);
fprintf('Saved det_estab_bpi_phib_rnfixed.pdf\n');
