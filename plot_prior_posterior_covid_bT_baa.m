% plot_prior_posterior_covid_bT_baa.m
% Prior vs. posterior density for the policy parameters bpi, by, ali,
% phib_val from the mmt8bOccBin3_covid_bT_baa production MCMC run.
% Posterior draws are pooled across both chains after a 50% burn-in
% (matching Dynare's default mh_drop=0.5). Prior densities are the exact
% analytical densities implied by the estimated_params block.

clear; close all;
shotgun_dir = 'C:/Users/gawater/Documents/GitHub/Shotgun 2';
cd(shotgun_dir);

metdir = fullfile(shotgun_dir, 'mmt8bOccBin3_covid_bT_baa', ...
    'mmt8bOccBin3_covid_bT_baa', 'metropolis');
S1 = load(fullfile(metdir, 'mmt8bOccBin3_covid_bT_baa_mh1_blck1.mat'), 'x2');
S2 = load(fullfile(metdir, 'mmt8bOccBin3_covid_bT_baa_mh1_blck2.mat'), 'x2');

% Column order in the raw x2 draws is Dynare's canonical MCMC ordering,
% NOT the estimated_params block's textual order: all 15 shock stderrs
% come first (cols 1-15: epsa,epsz,epsva,epsxa,epsg,epsr,epstau,epsL,
% eps_y_me,eps_d_me,eps_pig_me,eps_rl_me,eps_rn_me,eps_tau_me,eps_bT_me),
% then the 8 structural parameters (cols 16-23: va,eta,bpi,by,ali,br,
% rhoz,phib_val). Verified empirically by matching column means against
% oo_.posterior_mean.parameters (see conversation record) rather than
% trusting the block's textual order, which earlier produced mislabeled
% plots (shock-stderr draws shown as policy-parameter posteriors).
col = struct('bpi',18, 'by',19, 'ali',20, 'phib_val',23);

burn = 0.5;
n1 = size(S1.x2,1); n2 = size(S2.x2,1);
draws1 = S1.x2(round(burn*n1)+1:end, :);
draws2 = S2.x2(round(burn*n2)+1:end, :);
pooled = [draws1; draws2];
fprintf('Pooled post-burn-in draws: %d (chain1) + %d (chain2) = %d\n', ...
    size(draws1,1), size(draws2,1), size(pooled,1));

% Prior specs: {name, dynare_family, mean, std, plot_range}
% Updated for the diffuse-prior re-estimation round: bpi Normal(2.10,0.15)
% -> Gamma(1.00,1.00); by Gamma(0.28,0.10) -> Gamma(1.00,0.70) (reshaped
% off shape=1/Exponential after an earlier attempt put by's mode exactly
% at 0); phib_val mean lowered 0.10 -> 0.05. ali unchanged.
priors = struct( ...
    'bpi',      struct('family','gamma', 'mean',1.00, 'std',1.00, 'range',[0.0 7.0]), ...
    'by',       struct('family','gamma', 'mean',1.00, 'std',0.70, 'range',[0.0 2.0]), ...
    'ali',      struct('family','beta',  'mean',0.49, 'std',0.10, 'range',[0.0 1.0]), ...
    'phib_val', struct('family','beta',  'mean',0.05, 'std',0.05, 'range',[0.0 0.35]) ...
);

titles = struct('bpi','$b_\pi$ (inflation response)', 'by','$b_y$ (output response)', ...
    'ali','$a_i$ (policy-rate level weight)', 'phib_val','$\varphi_b$ (fiscal feedback)');
posterior_means = struct('bpi',3.4742,'by',0.0832,'ali',0.1808,'phib_val',0.0689);

fig = figure('Units','inches','Position',[0.5 0.5 11 8.5],'Color','w');
names = {'bpi','by','ali','phib_val'};
for p = 1:4
    nm = names{p};
    x  = pooled(:, col.(nm));
    pr = priors.(nm);

    subplot(2,2,p); hold on;
    [f_post, xi_post] = ksdensity(x);
    plot(xi_post, f_post, 'b-', 'LineWidth', 1.8);

    xg = linspace(pr.range(1), pr.range(2), 500);
    switch pr.family
        case 'normal'
            fprior = normpdf(xg, pr.mean, pr.std);
        case 'gamma'
            kshape = (pr.mean/pr.std)^2;
            theta  = pr.std^2/pr.mean;
            fprior = gampdf(xg, kshape, theta);
        case 'beta'
            v = pr.std^2;
            a = pr.mean*(pr.mean*(1-pr.mean)/v - 1);
            b = (1-pr.mean)*(pr.mean*(1-pr.mean)/v - 1);
            fprior = betapdf(xg, a, b);
    end
    plot(xg, fprior, 'r--', 'LineWidth', 1.5);

    xline(posterior_means.(nm), 'k:', 'LineWidth', 1.2);
    title(titles.(nm), 'Interpreter','latex', 'FontSize', 12);
    xlabel(nm, 'Interpreter','none'); ylabel('Density');
    if p==1
        legend({'Posterior','Prior','Posterior mean'}, 'Location','best', 'FontSize',8);
    end
    box on; hold off;
end

outfile = fullfile(shotgun_dir, 'prior_posterior_covid_bT_baa.pdf');
if exist('exportgraphics','file') || exist('exportgraphics','builtin')
    try exportgraphics(fig, outfile, 'ContentType','vector','BackgroundColor','white');
    catch, print(fig, outfile, '-dpdf', '-painters'); end
else
    print(fig, outfile, '-dpdf', '-painters');
end
fprintf('Saved: %s\n', outfile);
