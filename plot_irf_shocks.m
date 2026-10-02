% plot_irf_shocks.m
% Builds IRF graphs from oo_.irfs produced by mmt8bOccBin3_covid_bT_baa_irf.mod
% (posterior-mean linear IRFs). Run this immediately after that dynare call,
% in the same MATLAB session (needs oo_, M_ in the workspace).

shock_names  = {'epsa','epsz','epsxa','epsg'};
shock_labels = {'Preference wedge (\epsilon_a)', 'Productivity (\epsilon_z)', ...
                 'Deposit technology (\epsilon_{xa})', 'Govt. spending (\epsilon_g)'};

var_names  = {'y','pig','rn','rl','d'};
var_labels = {'Output, y', 'Inflation, \pi', 'Policy rate, r^n', 'Loan rate, r^l', 'Deposits, d'};

ss = struct();
for i = 1:numel(var_names)
    ss.(var_names{i}) = oo_.dr.ys(strcmp(M_.endo_names, var_names{i}));
end

horiz = size(oo_.irfs.([var_names{1} '_' shock_names{1}]), 2);
t = 0:horiz-1;

outdir = 'irf_graphs';
if ~exist(outdir, 'dir'); mkdir(outdir); end

% ---- Combined grid: rows = variables, columns = shocks ----
fig = figure('Position', [50 50 1400 1400], 'Color', 'w');
tl = tiledlayout(numel(var_names), numel(shock_names), 'TileSpacing', 'compact', 'Padding', 'compact');

for r = 1:numel(var_names)
    vn = var_names{r};
    for c = 1:numel(shock_names)
        sn = shock_names{c};
        nexttile;
        irf_level = oo_.irfs.([vn '_' sn]);
        irf_pct = 100 * irf_level / ss.(vn);
        plot(t, irf_pct, 'LineWidth', 1.8, 'Color', [0.05 0.35 0.65]);
        hold on; yline(0, ':', 'Color', [0.5 0.5 0.5]);
        grid on; box on; xlim([0 horiz-1]);
        if r == 1
            title(shock_labels{c}, 'FontSize', 10, 'FontWeight', 'bold');
        end
        if c == 1
            ylabel(var_labels{r}, 'FontSize', 10, 'FontWeight', 'bold');
        end
        if r == numel(var_names)
            xlabel('Quarters');
        end
    end
end
title(tl, 'Impulse Responses at the Posterior Mean (% deviation from steady state)', ...
    'FontSize', 13, 'FontWeight', 'bold');

exportgraphics(fig, fullfile(outdir, 'irf_grid_all_shocks.png'), 'Resolution', 200);
exportgraphics(fig, fullfile(outdir, 'irf_grid_all_shocks.pdf'));

% ---- One figure per shock, 5 stacked panels ----
for c = 1:numel(shock_names)
    sn = shock_names{c};
    figc = figure('Position', [50 50 500 1000], 'Color', 'w');
    tlc = tiledlayout(numel(var_names), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    for r = 1:numel(var_names)
        vn = var_names{r};
        nexttile;
        irf_level = oo_.irfs.([vn '_' sn]);
        irf_pct = 100 * irf_level / ss.(vn);
        plot(t, irf_pct, 'LineWidth', 1.8, 'Color', [0.05 0.35 0.65]);
        hold on; yline(0, ':', 'Color', [0.5 0.5 0.5]);
        grid on; box on; xlim([0 horiz-1]);
        ylabel(var_labels{r}, 'FontSize', 10);
        if r == numel(var_names)
            xlabel('Quarters');
        end
    end
    title(tlc, ['Shock: ' shock_labels{c}], 'Interpreter','tex', 'FontSize', 12, 'FontWeight', 'bold');
    fname = fullfile(outdir, ['irf_' sn '.png']);
    exportgraphics(figc, fname, 'Resolution', 200);
    exportgraphics(figc, fullfile(outdir, ['irf_' sn '.pdf']));
end

fprintf('Saved IRF graphs to %s\n', fullfile(pwd, outdir));
