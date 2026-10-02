% plot_mp_counterfactual.m
% Counterfactual: "what would pig and rn have looked like without the
% monetary-policy shock (epsr)?"
%
% Same additive-decomposition logic as plot_debt_counterfactual.m: since
% the smoother/decomposition is linear (order=1), the counterfactual is
% just the actual smoothed path minus epsr's own contribution. Run right
% after mmt8bOccBin3_covid_bT_baa_decomp.mod (needs oo_, M_ in workspace;
% also needs build_counterfactual, defined at the bottom of
% plot_debt_counterfactual.m -- run that script first, or just re-run this
% one which redefines its own copy below).

exo_names = cellstr(M_.exo_names);
idx_epsr = strmatch('epsr', exo_names, 'exact');

idx_pig = strmatch('pig', cellstr(M_.endo_names), 'exact');
idx_rn  = strmatch('rn',  cellstr(M_.endo_names), 'exact');

nT = size(oo_.shock_decomposition, 3);
start_year = 2009; start_q = 1;
qnum = (0:nT-1) + (start_q-1);
years = start_year + floor(qnum/4);
quarters = mod(qnum,4) + 1;
qlabels = arrayfun(@(y,q) sprintf('%dQ%d', y, q), years, quarters, 'UniformOutput', false);

[pig_actual, pig_cfact, ss_pig] = build_counterfactual_1shock(oo_, idx_pig, idx_epsr);
[rn_actual,  rn_cfact,  ss_rn ] = build_counterfactual_1shock(oo_, idx_rn,  idx_epsr);

% ---- Numbers for the write-up ----
window = 50:58; % 2021Q2 - 2023Q2
covid_window = 45:48; % 2020Q1 - 2020Q4

fprintf('=== Inflation (pig), 2021Q2-2023Q2 ===\n');
fprintf('avg actual dev from SS:        %.5f\n', mean(pig_actual(window)-ss_pig));
fprintf('avg counterfactual dev from SS: %.5f\n', mean(pig_cfact(window)-ss_pig));
fprintf('avg gap (actual - cfact):       %.5f (= avg epsr contribution)\n\n', mean(pig_actual(window)-pig_cfact(window)));

fprintf('=== Inflation (pig), 2020Q1-2020Q4 (COVID) ===\n');
fprintf('avg actual dev from SS:        %.5f\n', mean(pig_actual(covid_window)-ss_pig));
fprintf('avg counterfactual dev from SS: %.5f\n', mean(pig_cfact(covid_window)-ss_pig));
fprintf('avg gap (actual - cfact):       %.5f\n\n', mean(pig_actual(covid_window)-pig_cfact(covid_window)));

fprintf('=== Policy rate (rn), 2021Q2-2023Q2 ===\n');
fprintf('avg actual level:        %.5f\n', mean(rn_actual(window)));
fprintf('avg counterfactual level (rule-consistent, no epsr): %.5f\n', mean(rn_cfact(window)));
fprintf('avg gap (actual - cfact):       %.5f (annualized ~%.2f pp)\n\n', ...
    mean(rn_actual(window)-rn_cfact(window)), 400*mean(rn_actual(window)-rn_cfact(window)));

% ---- Plot: actual vs counterfactual, zoom 2019Q1-2024Q4 ----
outdir = 'irf_graphs';
if ~exist(outdir,'dir'); mkdir(outdir); end

zoom_lo = find(strcmp(qlabels,'2019Q1'));
zoom_hi = find(strcmp(qlabels,'2024Q4'));

fig = figure('Position',[50 50 1300 800],'Color','w');
tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile;
plot(zoom_lo:zoom_hi, pig_actual(zoom_lo:zoom_hi), 'LineWidth', 2.2, 'Color', [0.05 0.35 0.65]); hold on;
plot(zoom_lo:zoom_hi, pig_cfact(zoom_lo:zoom_hi), '--', 'LineWidth', 2.2, 'Color', [0.85 0.33 0.10]);
yline(ss_pig, ':', 'Color',[0.5 0.5 0.5]);
xline(window(1)-0.5,'--','Color',[0.3 0.3 0.3],'Label','2021Q2');
xline(window(end)+0.5,'--','Color',[0.3 0.3 0.3],'Label','2023Q2');
xt2 = zoom_lo:2:zoom_hi;
set(gca,'XTick',xt2,'XTickLabel',qlabels(xt2),'XTickLabelRotation',45);
legend({'actual (smoothed \pi)','counterfactual: \epsilon_r shut off','steady state'}, 'Location','best');
ylabel('\pi (gross quarterly)');
title('Inflation Around the Post-COVID Spike: Actual vs. No-MP-Shock Counterfactual');
grid on; box on; xlim([zoom_lo-0.5 zoom_hi+0.5]);

nexttile;
plot(zoom_lo:zoom_hi, rn_actual(zoom_lo:zoom_hi), 'LineWidth', 2.2, 'Color', [0.05 0.35 0.65]); hold on;
plot(zoom_lo:zoom_hi, rn_cfact(zoom_lo:zoom_hi), '--', 'LineWidth', 2.2, 'Color', [0.85 0.33 0.10]);
yline(ss_rn, ':', 'Color',[0.5 0.5 0.5]);
xline(window(1)-0.5,'--','Color',[0.3 0.3 0.3]);
xline(window(end)+0.5,'--','Color',[0.3 0.3 0.3]);
set(gca,'XTick',xt2,'XTickLabel',qlabels(xt2),'XTickLabelRotation',45);
legend({'actual (smoothed r^n)','counterfactual: rule-consistent (\epsilon_r shut off)','steady state'}, 'Location','best');
ylabel('Policy rate, r^n (gross quarterly)');
xlabel('Quarter');
title('Policy Rate: Actual vs. What the Estimated Rule Alone Would Set');
grid on; box on; xlim([zoom_lo-0.5 zoom_hi+0.5]);

exportgraphics(fig, fullfile(outdir,'mp_counterfactual_zoom.png'), 'Resolution', 200);
exportgraphics(fig, fullfile(outdir,'mp_counterfactual_zoom.pdf'));

% ---- Full sample ----
fig2 = figure('Position',[50 50 1400 800],'Color','w');
tl2 = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile;
plot(1:nT, pig_actual, 'LineWidth', 2, 'Color', [0.05 0.35 0.65]); hold on;
plot(1:nT, pig_cfact, '--', 'LineWidth', 2, 'Color', [0.85 0.33 0.10]);
yline(ss_pig, ':', 'Color',[0.5 0.5 0.5]);
xt = 1:8:nT;
set(gca,'XTick',xt,'XTickLabel',qlabels(xt),'XTickLabelRotation',45);
legend({'actual (smoothed \pi)','counterfactual: \epsilon_r shut off','steady state'}, 'Location','best');
ylabel('\pi (gross quarterly)');
title('Inflation: Actual vs. Counterfactual Without the Monetary-Policy Shock');
grid on; box on; xlim([0.5 nT+0.5]);

nexttile;
plot(1:nT, rn_actual, 'LineWidth', 2, 'Color', [0.05 0.35 0.65]); hold on;
plot(1:nT, rn_cfact, '--', 'LineWidth', 2, 'Color', [0.85 0.33 0.10]);
yline(ss_rn, ':', 'Color',[0.5 0.5 0.5]);
set(gca,'XTick',xt,'XTickLabel',qlabels(xt),'XTickLabelRotation',45);
legend({'actual (smoothed r^n)','counterfactual: rule-consistent (\epsilon_r shut off)','steady state'}, 'Location','best');
ylabel('Policy rate, r^n (gross quarterly)');
xlabel('Quarter');
title('Policy Rate: Actual vs. Rule-Consistent Path (\epsilon_r shut off)');
grid on; box on; xlim([0.5 nT+0.5]);

exportgraphics(fig2, fullfile(outdir,'mp_counterfactual_full.png'), 'Resolution', 200);
exportgraphics(fig2, fullfile(outdir,'mp_counterfactual_full.pdf'));

fprintf('Saved MP counterfactual graphs to %s\n', fullfile(pwd, outdir));

function [actual, cfact, ss] = build_counterfactual_1shock(oo_, idx_var, idx_s1)
    total = squeeze(oo_.shock_decomposition(idx_var,17,:));
    c1 = squeeze(oo_.shock_decomposition(idx_var,idx_s1,:));
    ss = oo_.dr.ys(idx_var);
    actual = ss + total;
    cfact  = ss + (total - c1);
end
