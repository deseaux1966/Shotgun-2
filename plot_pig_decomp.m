% plot_pig_decomp.m
% Historical shock decomposition of pig (gross quarterly inflation) around
% the post-COVID inflation spike. Run right after
% mmt8bOccBin3_covid_bT_baa_decomp.mod (needs oo_, M_ in the workspace).

idx = strmatch('pig', cellstr(M_.endo_names), 'exact');
shock_names_all = cellstr(M_.exo_names);
nT = size(oo_.shock_decomposition, 3);

% Quarter labels: t=1 is 2009Q1 (per estimation .mod header)
start_year = 2009; start_q = 1;
qnum = (0:nT-1) + (start_q-1);
years = start_year + floor(qnum/4);
quarters = mod(qnum,4) + 1;
qlabels = arrayfun(@(y,q) sprintf('%dQ%d', y, q), years, quarters, 'UniformOutput', false);

% Group small/irrelevant shocks together for readability: keep the 8
% structural shocks separate, lump the 7 measurement-error shocks + the
% initial-condition column into one "Meas. error / init." bar segment.
struct_idx = 1:8; % epsa epsz epsva epsxa epsg epsr epstau epsL
other_idx  = [9:15, 16]; % measurement-error shocks + initial conditions

series = squeeze(oo_.shock_decomposition(idx, struct_idx, :))'; % nT x 8
other_series = sum(squeeze(oo_.shock_decomposition(idx, other_idx, :)), 1)'; % nT x 1
plot_data = [series, other_series];
plot_labels = [shock_names_all(struct_idx); {'meas. error / init.'}];

% Nicer labels
nice_labels = {'\epsilon_a (pref. wedge)','\epsilon_z (productivity)','\epsilon_{va} (shopping-time pref.)', ...
    '\epsilon_{xa} (deposit tech.)','\epsilon_g (govt. spending)','\epsilon_r (MP shock)', ...
    '\epsilon_\tau (tax rule)','\epsilon_L (labor supply)','meas. error / init.'};

% ---- Full-sample stacked bar with spike window highlighted ----
fig = figure('Position',[50 50 1500 650],'Color','w');
window_lo = 50; window_hi = 58; % 2021Q2 - 2023Q2
b = bar(1:nT, plot_data, 'stacked'); hold on;
plot(1:nT, squeeze(oo_.shock_decomposition(idx,17,:)), 'k-', 'LineWidth', 1.8);
xline(window_lo-0.5, '--', 'Color',[0.3 0.3 0.3]);
xline(window_hi+0.5, '--', 'Color',[0.3 0.3 0.3]);
xt = 1:8:nT;
set(gca, 'XTick', xt, 'XTickLabel', qlabels(xt), 'XTickLabelRotation', 45);
legend([nice_labels, {'total (smoothed pig, dev. from SS)'}], 'Location','eastoutside');
ylabel('Deviation of \pi (gross quarterly) from steady state');
title('Historical Shock Decomposition of Inflation, 2009Q1-2025Q4 (posterior mean, linear KF smoother)');
grid on; box on;
xlim([0.5 nT+0.5]);

outdir = 'irf_graphs';
if ~exist(outdir,'dir'); mkdir(outdir); end
exportgraphics(fig, fullfile(outdir,'pig_shock_decomp_full.png'), 'Resolution', 200);
exportgraphics(fig, fullfile(outdir,'pig_shock_decomp_full.pdf'));

% ---- Zoomed window: 2019Q1 - 2024Q4 ----
zoom_lo = find(strcmp(qlabels,'2019Q1'));
zoom_hi = find(strcmp(qlabels,'2024Q4'));
fig2 = figure('Position',[50 50 1300 650],'Color','w');
bar(zoom_lo:zoom_hi, plot_data(zoom_lo:zoom_hi,:), 'stacked'); hold on;
plot(zoom_lo:zoom_hi, squeeze(oo_.shock_decomposition(idx,17,zoom_lo:zoom_hi)), 'k-', 'LineWidth', 2);
xline(window_lo-0.5, '--', 'Color',[0.3 0.3 0.3], 'Label','2021Q2');
xline(window_hi+0.5, '--', 'Color',[0.3 0.3 0.3], 'Label','2023Q2');
xt2 = zoom_lo:2:zoom_hi;
set(gca, 'XTick', xt2, 'XTickLabel', qlabels(xt2), 'XTickLabelRotation', 45);
legend([nice_labels, {'total (smoothed pig, dev. from SS)'}], 'Location','eastoutside');
ylabel('Deviation of \pi (gross quarterly) from steady state');
title('Shock Decomposition of Inflation Around the Post-COVID Spike (2019Q1-2024Q4)');
grid on; box on;
xlim([zoom_lo-0.5 zoom_hi+0.5]);

exportgraphics(fig2, fullfile(outdir,'pig_shock_decomp_zoom.png'), 'Resolution', 200);
exportgraphics(fig2, fullfile(outdir,'pig_shock_decomp_zoom.pdf'));

fprintf('Saved decomposition graphs to %s\n', fullfile(pwd, outdir));
