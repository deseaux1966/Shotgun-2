% plot_zlb_regimes.m
% ZLB-binding quarters implied by the estimation: OccBin smoother and
% real-time (filter) regime histories vs. a model-independent check based
% on the linear smoother's notional rate falling below the floor of 1.

S = load('C:\Users\gawater\Documents\GitHub\Shotgun-2\mmt8bOccBin3_covid_bT_baa\mmt8bOccBin3_covid_bT_baa\Output\results_covid_bT_baa_phib5000.mat');
oc = S.oo_.occbin;
T = 68;
bind_smooth = false(T,1); bind_rt = false(T,1);
for t = 1:T
    bind_smooth(t) = oc.smoother.regime_history(t).regime(1) == 1;
    bind_rt(t)     = oc.smoother.realtime_regime_history(t).regime(1) == 1;
end
notional_lin = oc.linear_smoother.SmoothedVariables.rn_notional(1:T);
bind_lin = notional_lin(:) < 1;
notional_occ = oc.smoother.SmoothedVariables.rn_notional(1:T);
rn_occ = oc.smoother.SmoothedVariables.rn(1:T);

qnum = 0:T-1; years = 2009 + floor(qnum/4); qs = mod(qnum,4)+1;
ql = arrayfun(@(y,q) sprintf('%dQ%d',y,q), years, qs, 'UniformOutput', false);

function spells(b, ql)
    d = diff([0; b(:); 0]); s = find(d==1); e = find(d==-1)-1;
    if isempty(s), fprintf('   none\n'); end
    for i = 1:numel(s), fprintf('   %s - %s (%d qtrs)\n', ql{s(i)}, ql{e(i)}, e(i)-s(i)+1); end
end
fprintf('OccBin smoother regime_history:\n');   spells(bind_smooth, ql);
fprintf('OccBin real-time (filter) history:\n'); spells(bind_rt, ql);
fprintf('Linear smoother, notional rn < 1:\n');  spells(bind_lin, ql);
fprintf('Disagreements realtime vs linear: %d qtrs; smoother vs linear: %d qtrs\n', ...
    sum(bind_rt~=bind_lin), sum(bind_smooth~=bind_lin));

fig = figure('Position',[50 50 1300 700],'Color','w');
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile; hold on;
plot(1:T, notional_lin, 'Color',[0.05 0.35 0.65],'LineWidth',1.8);
plot(1:T, notional_occ, '--','Color',[0.85 0.33 0.10],'LineWidth',1.5);
yline(1,'k-'); ylim_ = ylim;
for t = find(bind_lin)', patch([t-.5 t+.5 t+.5 t-.5],[ylim_(1) ylim_(1) ylim_(2) ylim_(2)],[0.5 0.5 0.5],'FaceAlpha',0.25,'EdgeColor','none'); end
set(gca,'XTick',1:8:T,'XTickLabel',ql(1:8:T),'XTickLabelRotation',45);
legend({'notional rate, linear smoother','notional rate, OccBin smoother','ZLB floor = 1'},'Location','best');
ylabel('r^n notional (gross quarterly)'); grid on; box on; xlim([0.5 T+.5]);
title('Notional policy rate vs. ZLB floor (shaded: binding per linear smoother)');
nexttile; hold on;
stairs(1:T, bind_rt+0.04,'LineWidth',1.8,'Color',[0.05 0.35 0.65]);
stairs(1:T, bind_smooth+0.02,'--','LineWidth',1.5,'Color',[0.85 0.33 0.10]);
stairs(1:T, bind_lin,'-','LineWidth',1.2,'Color',[0.2 0.6 0.2]);
set(gca,'XTick',1:8:T,'XTickLabel',ql(1:8:T),'XTickLabelRotation',45,'YTick',[0 1],'YTickLabel',{'relaxed','binding'});
legend({'OccBin real-time (filter)','OccBin smoother','linear smoother notional<1'},'Location','best');
grid on; box on; xlim([0.5 T+.5]); ylim([-0.1 1.15]);
title('ZLB regime by quarter, three measures');
outdir = 'C:\Users\gawater\Documents\GitHub\Shotgun-2\irf_graphs';
exportgraphics(fig, fullfile(outdir,'zlb_regimes.png'),'Resolution',200);
exportgraphics(fig, fullfile(outdir,'zlb_regimes.pdf'));
