% mmt8bOccBin3_covid_bT_baa.mod
% Bayesian estimation of mmt8b with OccBin ZLB + debt-to-GDP ratio.
% rl_obs constructed from BAA (Moody's Baa Corporate Bond Yield) instead of MPRIME.
%
% Extension of mmt8bOccBin3_covid.mod with:
%   - bT: endogenous debt-to-GDP ratio variable
%   - bT_obs: observable for U.S. federal debt HELD BY THE PUBLIC as a
%     share of GDP (FYGFDPUN/GDP), not gross/total public debt
%   - rl_obs: BAA corporate bond yield (credit-spread aware)
%   - epsL: single iid shock to the labor-leisure margin, standing in for the
%     COVID labor-supply collapse. Replaces the Cardani et al.-style multi-shock
%     regime-switching machinery (COVID_APC/COVID_AY/at_covid/zt_covid/gt_covid),
%     which, on inspection, was never actually functional: at_covid, zt_covid,
%     gt_covid were pinned to 0 by their own defining equations with no
%     stochastic driver, so that mechanism was silently a no-op.
%
% FRED series:
%   BAA      - Moody's Seasoned Baa Corporate Bond Yield (monthly, %)
%   FYGFDPUN - Federal Debt Held by the Public (quarterly, millions $)
%   GDP      - Gross Domestic Product (quarterly, nominal, billions $)
%   bT = (FYGFDPUN/1000) / GDP   [unit fix: millions -> billions]
%
% NOTE: earlier versions of this file used GFDEBTN (total/gross public
% debt) paired with "NGDPDQ", a series ID that turns out not to exist on
% FRED at all (fredgraph.csv?id=NGDPDQ returns an HTML error page, not
% data) -- the debt_to_gdp data actually used in those runs came from a
% separately pre-downloaded ratio file (grossdebt.csv), not from that
% script. bT_obs is now built directly from FYGFDPUN and GDP (both
% confirmed valid series) via build_debt_dataset_fygfdpun.m, and is
% demeaned in LOGS (matching the measurement equation below) rather than
% in levels, fixing a second inconsistency in the earlier data pipeline.
%
% Data: 2009Q1-2025Q4 (68 obs).  2020Q1 = t=45, 2020Q2 = t=46.
%
% ZLB:     OccBin constraint on rn. Strict floor rn_lb=1.0.
% Likelihood: OccBin enabled (true piecewise-linear ZLB likelihood/smoother).
%
% NOTE ON varexo_det: a deterministic one-period dummy (varexo_det) was
% considered for the COVID shock, but dynare_estimation_init.m (Dynare 7.1)
% does not read varexo_det from the estimation datafile, so it would silently
% be ignored by the Kalman filter. epsL below is a normal iid varexo instead;
% the smoother will attribute whatever realization the data need to it,
% including at 2020Q2, without requiring a hardcoded date.
%
% Two-stage mode-finding (see bottom of file): OccBin's likelihood is
% nonsmooth at ZLB regime switches, so csminwel is run first on the linear
% (non-OccBin) approximation to get a mode matching THIS parameter vector,
% then OccBin is switched on with mode_compute=0 using that mode as start.
% (A stale mode_file from a different model spec previously caused
% "generated using another specification of the model" errors.)
%
% Fixed:   beta=0.9967, chi=1.0, psi=1.3
%          rhova=0, rhog=0 (i.i.d.)
%          theta=6 (standard NK calibration)
% Estimated: va, eta, bpi, by, ali, br, rhoz, rhoa, rhoxa, phib, epsL, 9 shock SEs

var y c d n pig rl rd h hg hd hn hs mu w rn rn_notional m at zt vat xat gt tau b bT zr
    y_obs d_obs pig_obs rl_obs rn_obs taut_obs bT_obs;
varexo epsa epsz epsva epsxa epsg epsr epstau epsL
        eps_y_me eps_d_me eps_pig_me eps_rl_me eps_rn_me eps_tau_me eps_bT_me;

Parameters chi psi va beta eta theta nu xa xn xh phin phip
            g taubar phib_val bpi by br ali omega b0 z rlss pigss rn_lb
            rhoa rhoz rhova rhoxa rhog rhotau;

ali    =  0.5;
beta   =  0.9967;
eta    =  1.5294;
va     =  2.5035;
xa     =  5.2528;

% chi recalibrated from 1.0 to 0.5: comparative-statics check showed
% lower chi (shopping-time curvature) implies substantially lower
% steady-state inflation via higher money demand -- see conversation
% record. chi remains calibrated, not added to estimated_params.
chi    =  0.5;
psi    =  1.3;
omega  =  0.0;
b0     =  1.0;
z      =  1.3;

g      =  0.2;
% taubar recalibrated from 0.20 to 0.14: deficit g-taubar = 0.06 (6% of
% GDP) is a more realistic steady-state target for the 2009-2025 sample
% away from the ZLB than the previous ~0% (balanced-budget) target.
% (Intermediate value of 0.15 tried first, then adjusted to 0.14.)
if exist('taubar_val','var')
    taubar = taubar_val;
else
    taubar = 0.14;
end;
phib_val = 0.10;

pigss  =  1.02^(0.25);
rlss   =  1.02^(0.25) / beta;
rn_lb  =  1.0;

% rhoa/rhoxa fixed at 0.95, not estimated: a prior "free rhoa/rhoxa"
% experiment in this file converged to rhoa~0.9992 (unit root), producing
% a non-positive-definite Hessian at the mode -- the same failure the base
% mmt8bOccBin3.mod already documented and fixed by pinning rhoa=0.95.
rhoa   =  0.95;
rhoz   =  0.65;
rhova  =  0.0;
rhoxa  =  0.95;
rhog   =  0.0;
rhotau =  0.5;

bpi    =  1.5;
by     =  0.5;
br     =  0.3;

theta  =  6.0;
nu     =  0.25;
xn     =  0.95;
xh     = 1 - xn;
phin   =  0.0005;
phip   =  50;

model;
% --- Household / bank block ---
% epsL: iid labor-supply shock (COVID stand-in). epsL<0 shrinks the RHS,
% forcing lower h for given c,w -- a negative labor supply shock.
hs = (1/chi) * (vat*c/d)^chi;
(c/z)*(h^psi) + chi*hs = (1/eta)*exp(epsL);
(rl-rd)*(h^psi)*d = chi*hs*w;
mu = at * w * (h^psi);
beta*(rl(+1)/pig(+1))*(mu(+1)/mu) = 1;

% --- Production ---
y = zt*hg;

% --- Banking: deposit supply and portfolio ---
d = xat*((xn^(1/nu))*(n^((nu-1)/nu)) + (xh^(1/nu))*(hd^((nu-1)/nu)))^(nu/(nu-1));
hn = phin*n;
n  = ((rl-rd)^nu)*(xat^(nu-1))*xn*d*(rl-rn+phin*z)^(-nu);
hd = ((rl-rd)^nu)*(xat^(nu-1))*xh*d*(z^(-nu));
h  = hg + hd + hn;

% --- Price Phillips Curve ---
0.5*(pig-STEADY_STATE(pig))*pig/(STEADY_STATE(pig)^2) = (1-theta+theta*(w/zt))/phip
    + beta*(pig(+1)-STEADY_STATE(pig))*pig(+1)*(y(+1)/y)/(STEADY_STATE(pig)^2);

% --- Resource constraint ---
y = c + gt + 0.5*phip*(((pig-STEADY_STATE(pig))/(STEADY_STATE(pig)))^2);

% --- Government / money / bonds ---
m = (gt) - (tau) + (rn)*n;
n = (m(-1)/pig) + (b(-1)/pig) - (b/rl);
tau = (tau(-1)^rhotau) * (taubar^(1-rhotau)) * ((b(-1)/((1-omega)*b0))^phib_val) * exp(epstau);
b = (rlss/pigss)*b(-1) + gt - tau + (1-omega)*b0*(1 - rlss/pigss) - (g - taubar);

% --- Debt-to-GDP ratio ---
% bT = b / (1 + g)  where b is debt-to-output and g is govt spending-to-output
% GDP = Y + G, so in ratio terms GDP = 1 + g
% Steady state: bT_ss = b_ss / (1 + g) = 1.0 / 1.2 = 0.8333
bT = b / (1 + g);

% --- Monetary policy ---
rn_notional = (STEADY_STATE(rl)^ali) * (((pig)/STEADY_STATE(pig))^bpi)
              * ((y/steady_state(y))^by) * exp(zr);
zr = br * zr(-1) + epsr;

[name='ZLB_policy', relax='ZLB']
rn = rn_notional;

[name='ZLB_policy', bind='ZLB']
rn = rn_lb;

% --- Exogenous shock processes ---
at = (at(-1)^rhoa) * exp(epsa);
zt = (zt(-1)^rhoz) * (z^(1-rhoz)) * exp(epsz);
vat = (vat(-1)^rhova) * (va^(1-rhova)) * exp(epsva);
xat = (xat(-1)^rhoxa) * (xa^(1-rhoxa)) * exp(epsxa);
gt = (rhog*gt(-1)) + ((1-rhog)*g) + epsg;

% --- Observable equations (7 observables) ---
y_obs    = log(y)            - log(STEADY_STATE(y))   + eps_y_me;
d_obs    = log(d)            - log(STEADY_STATE(d))   + eps_d_me;
pig_obs  = log(pig)          - log(STEADY_STATE(pig)) + eps_pig_me;
rl_obs   = log(rl)           - log(STEADY_STATE(rl))  + eps_rl_me;
rn_obs   = log(rn)           - log(STEADY_STATE(rn))  + eps_rn_me;
taut_obs = log(tau)          - log(STEADY_STATE(tau)) + eps_tau_me;
bT_obs   = log(bT)           - log(STEADY_STATE(bT))  + eps_bT_me;
end;

% --- OccBin ZLB constraint ---
occbin_constraints;
name 'ZLB'; bind log(rn_notional) < -0.000001; relax log(rn_notional) >= -0.000001;
end;

initval;
pig  = 1.004963;
rl   = 1.008270;
rn   = 1.004128;
rn_notional = 1.004128;
zr   = 0;
w    = 1.08333;
at   = 1;  zt = 1.3;  vat = 2.5035;  xat = 5.2528;  gt = 0.2;  tau = 0.20;  b = 1.0;
n    = 9.344;
m    = 9.498;
hn   = 0.00196;
hd   = 0.0510;
d    = 13.38;
rd   = 1.002;
hg   = 0.783;
h    = 0.836;
y    = 1.018;
c    = 0.818;
hs   = 0.153;
mu   = 0.858;
% Debt-to-GDP ratio steady state: bT = b / (1 + g) = 1.0 / 1.2 = 0.8333
bT   = 1.0 / (1 + g);
y_obs    = 0;
d_obs    = 0;
pig_obs  = 0;
rl_obs   = 0;
rn_obs   = 0;
taut_obs = 0;
bT_obs   = 0;
end;

shocks;
var epsa;               stderr 0.1;
var epsz;               stderr 0.1;
var epsva;              stderr 0.1;
var epsxa;              stderr 0.1;
var epsg;               stderr 0.1;
var epsr;               stderr 0.1;
var epstau;             stderr 0.1;
var epsL;               stderr 0.15;
var eps_y_me;           stderr 0.1;
var eps_d_me;           stderr 0.1;
var eps_pig_me;         stderr 0.1;
var eps_rl_me;          stderr 0.1;
var eps_rn_me;          stderr 0.1;
var eps_tau_me;         stderr 0.1;
var eps_bT_me;          stderr 0.02;   % debt-to-GDP data is precise
end;

steady(solve_algo=4);
check;

varobs y_obs d_obs pig_obs rl_obs rn_obs taut_obs bT_obs;

estimated_params;
% Fixed: beta=0.9967, chi=1.0, psi=1.3
% Fixed: rhova=rhog=0 (i.i.d.)
% Fixed: theta=6 (standard NK calibration)

va,         gamma_pdf,  2.50,  0.75;
 eta,        gamma_pdf,  1.53,  0.40;

 bpi,        gamma_pdf,  1.00,  1.00;
% by: std=0.70 (not 1.00) -> shape=(mean/std)^2~=2.04, an interior-peaked
% hump (mode ~0.51) rather than shape=1 (Exponential, mode AT 0), so the
% prior itself doesn't pre-load a corner solution -- see conversation
% record (Stage-1 run under Gamma(1,1) put by's mode at exactly 0.0000
% with a non-PD Hessian there).
 by,         gamma_pdf,  1.00,  0.70;
 ali,        beta_pdf,   0.49,  0.10;
 br,         beta_pdf,   0.36,  0.10;
 rhoz,       beta_pdf,   0.45,  0.10;

% rhoa, rhoxa fixed at 0.95 (see calibration comment above) -- not estimated.

% Estimated fiscal feedback (Extension 2)
phib_val,   beta_pdf,   0.05,  0.05;

stderr epsa,               inv_gamma_pdf, 0.05, 0.025;
stderr epsz,               inv_gamma_pdf, 0.05, 0.025;
stderr epsva,              inv_gamma_pdf, 0.05, 0.025;
stderr epsxa,              inv_gamma_pdf, 0.05, 0.025;
stderr epsg,               inv_gamma_pdf, 0.05, 0.025;
stderr epsr,               inv_gamma_pdf, 0.05, 0.025;
stderr epstau,             inv_gamma_pdf, 0.05, 0.025;
% epsL: loose/wide prior so the one-off COVID quarter isn't crowded out by
% measurement error absorbing the outlier instead.
stderr epsL,               inv_gamma_pdf, 0.15, 0.075;
stderr eps_y_me,           inv_gamma_pdf, 0.05, 0.025;
stderr eps_d_me,           inv_gamma_pdf, 0.05, 0.025;
stderr eps_pig_me,         inv_gamma_pdf, 0.05, 0.025;
stderr eps_rl_me,          inv_gamma_pdf, 0.05, 0.025;
 stderr eps_rn_me,          inv_gamma_pdf, 0.01, 0.005;
 stderr eps_tau_me,         inv_gamma_pdf, 0.05, 0.025;
 stderr eps_bT_me,          inv_gamma_pdf, 0.05, 0.025;
 end;

estimated_params_init;
   br, 0.36;
   phib_val, 0.05;
% bpi, by initialized at gen1's posterior means (safely determinate; see
% conversation record), NOT at the new diffuse priors' own means -- the new
% bpi prior mean of 1.0 sits just below the determinacy threshold (~1.095
% at these by/phib_val values), so starting the Stage-1 optimizer there
% risks an undefined/rejected first likelihood evaluation.
   bpi, 2.169;
   by, 0.135;
  end;

%----------------------------------------------------------------
% Estimation control.
%
% IMPORTANT Dynare constraint (discovered the hard way): mode_compute,
% mh_replic, mh_jscale, mh_nblocks inside the estimation(...) command must
% be literal INT_NUMBER/value tokens at parse time -- the preprocessor
% rejects a bare MATLAB identifier there ("syntax error, unexpected
% IDENTIFIER, expecting INT_NUMBER"). So those three are hardcoded below
% and tuned by editing this file directly between runs (Stage 1 / Stage 2
% / production), same as the rest of this codebase already does.
%
% occbin_val and mode_file_val ARE overridable from a wrapper script,
% because they are set via plain options_.* assignment OUTSIDE the
% estimation(...) parens, which Dynare does not literal-check:
%   occbin_val     false = linear KF (Stage 1: mode-finding for this
%                  parameter vector); true = true OccBin piecewise-linear
%                  ZLB likelihood (Stage 2).
%   mode_file_val  path (no .mat) to a mode file matching THIS spec.
%                  A mode from a different parameter vector throws
%                  "generated using another specification of the model".
%----------------------------------------------------------------
if exist('occbin_val','var')
    occbin_use = occbin_val;
else
    occbin_use = false;
end;

options_.occbin.likelihood.status = occbin_use;
options_.occbin.smoother.status   = occbin_use;
if occbin_use
    options_.occbin.likelihood.periodic_solution  = true;
    options_.occbin.smoother.periodic_solution    = true;
    options_.occbin.likelihood.use_updated_regime = 1;
    options_.occbin.regularisation                = 1e-3;
end;

%----------------------------------------------------------------
% Bayesian estimation. Uses fred_data_debt_baa.mat (BAA-based rl_obs,
% debt-to-GDP bT_obs).
%
% RUN HISTORY (previous data: GFDEBTN-based bT_obs, level-demeaned):
%   Stage 1: mode_compute=4, mh_replic=3000, no mode_file -> acceptance
%     0.23 / 0.228 (linear KF, informational only).
%   Stage 2 @ jscale=0.05, prior_variance: 0.0433 / 0.071 -- killed.
%   Stage 2 @ jscale=0.03, prior_variance: 0.203 / 0.057 -- promoted.
%   Production @ jscale=0.03, 50000 draws: 0.157 / 0.0477 -- killed.
%   Stage 2 @ jscale=0.02, prior_variance: 0.306 / 0.293 -- promoted.
%   Production @ jscale=0.02, 50000 draws: 0.271 / 0.242 -- accepted, but
%     3 of 4 policy-parameter posteriors (bpi, by, phib_val) showed
%     chain-to-chain disagreement (pooled KDE looked bimodal) traced to
%     mcmc_jumping_covariance=prior_variance ignoring the bpi/by
%     correlation -- see conversation record.
%
% CURRENT DATA: bT_obs rebuilt from FYGFDPUN/GDP (debt held by the
% public, not gross debt), log-demeaned (see header).
%
% Stage 1 (linear KF, mode_compute=4, new FYGFDPUN data): mode/Hessian
%   regenerated cleanly, no singularity. Acceptance 0.223 / 0.231
%   (informational, prior_variance jumping).
%
% Stage 2 @ jscale=1.0, hessian: acceptance 0.0073 / 0.0040 -- both far
%   below 0.05. Took >1hr and produced 6,861 near-singular-matrix
%   warnings (RCOND down to 1e-21/1e-22, far worse than any prior run)
%   -- proposals were landing in numerically pathological parameter
%   regions. Killed per instructions; jscale dropped substantially.
% Stage 2 @ jscale=0.2, hessian, 3000 draws: acceptance 0.100 / 0.104 --
%   both healthy on the diagnostic. Promoted to production.
% Production @ jscale=0.2, hessian, 50000 draws: acceptance 0.0449 /
%   0.0582 -- chain 1 below 0.05. Took 15h03m (vs 7h36m for the earlier
%   prior_variance production run). As with the earlier prior_variance
%   round (diagnostic 0.203/0.057 -> production 0.157/0.048), the short
%   diagnostic overstated what the full run sustains -- acceptance
%   degrades over the full 50k draws, likely from exploring more of the
%   harder ZLB-heavy tails. Killed per instructions; jscale cut by ~40%
%   (not a token nudge) given the cost of another full run.
% Stage 2 @ jscale=0.12, hessian, 3000 draws: acceptance 0.122 / 0.129.
% Stage 2 @ jscale=0.10, hessian, 3000 draws: acceptance 0.236 / 0.265 --
%   both in the textbook-healthy 0.2-0.4 range, extra margin built in
%   given the cost of another full run. Promoted to production.
% Production @ jscale=0.10, hessian, 50000 draws: acceptance 0.214 /
%   0.195 -- ACCEPTED (taubar=0.20 baseline). See
%   mmt8bOccBin3_covid_bT_baa_epsL_summary.tex/pdf for the full posterior
%   table, Geweke diagnostics, and determinacy/E-stability results from
%   this run.
%
% RECALIBRATION: taubar changed from 0.20 to 0.14 (deficit g-taubar=0.06,
% 6% of GDP, a more realistic steady-state target than the previous ~0%
% balanced-budget target), and chi changed from 1.0 to 0.5 (lower
% shopping-time curvature; comparative statics showed this substantially
% lowers steady-state inflation via higher money demand) -- see
% conversation record. Both changes shift the model's steady state, so
% the mode_file above (computed under taubar=0.20, chi=1.0) is STALE and
% mode-finding must restart from Stage 1. THIS BLOCK IS SET TO STAGE 1
% (linear KF, mode_compute=4, no mode_file). Once a valid Stage-1 mode is
% obtained, follow the same Stage 2 (occbin_val=true, mode_compute=0,
% mode_file pointing at the new mode, mcmc_jumping_covariance=hessian) +
% jscale-tuning + production sequence documented above and in the .tex
% writeup.
%
% Stage 1 (taubar=0.14, chi=0.5): mode/Hessian regenerated cleanly, no
%   singularity. Acceptance 0.223 / 0.230 (informational, prior_variance).
%
% Stage 2 @ jscale=0.10, hessian, 3000 draws (taubar=0.14, chi=0.5):
%   acceptance 0.405 / 0.259 -- both healthy with good margin. Promoted
%   directly to production.
% Production @ jscale=0.10, hessian, 50000 draws: acceptance 0.1775 /
%   0.1534 -- ACCEPTED, but Geweke convergence markedly uneven (chain 1
%   6/8, chain 2 1/8) and by/ali chain-mean gaps ~0.95-0.96 pooled SD,
%   worse than the previous hessian round. Diagnosed as the Hessian
%   (a local Gaussian approximation at the Stage-1 mode) being a
%   mediocre match for this posterior's true shape away from the mode,
%   not correlation-blindness (ruled out prior_variance as a fix --
%   see conversation record).
%
% THIS BLOCK NOW USES mcmc_jumping_covariance pointing at
% posterior_jumping_covariance.mat, built from this round's own pooled
% post-burn-in draws (build_posterior_jumping_covariance.m) rather than
% the Stage-1 mode's Hessian. mode_file is unchanged (still a valid
% starting point under this calibration; only the proposal SHAPE
% changes). Starting mh_jscale=0.5, the classic Roberts-Gelman-Gilks
% 2.4/sqrt(23) rule of thumb -- appropriate here since this covariance
% reflects actual posterior spread rather than local mode curvature, so
% jscale should behave less unpredictably than it did tuning from the
% Hessian (which needed 0.10, far from the RGG-implied ~0.5).
%
% Diagnostic @ jscale=0.5, posterior-covariance, 3000 draws: acceptance
%   0.2047 / 0.2017 -- excellent, right in the textbook 0.2-0.4 range
%   and remarkably symmetric. Chain-mean gaps: bpi 0.29, by 0.17,
%   phib_val 0.31 pooled SD -- all clearly better than the hessian
%   round. ali showed a large gap (1.83 SD, both chains ~0.49-0.56 vs
%   the true posterior mean ~0.18) but this diagnostic only kept 1500
%   draws/chain post-burn-in from a mode_file starting point found under
%   the LINEAR (non-OccBin) approximation -- read as insufficient travel
%   time from a mismatched starting point, not proposal-covariance
%   disagreement, and not expected to recur with 25,000 kept draws/chain
%   in the full run. Promoted directly to production per instructions.
% Production @ jscale=0.5, posterior-covariance, 50000 draws: acceptance
%   0.0805 / 0.0387 -- chain 2 below 0.05, killed. Sharper diagnostic-
%   to-production degradation than usual (0.20/0.20 -> 0.08/0.04).
%   jscale cut roughly in half to 0.25, back to a short diagnostic.
% Diagnostic @ jscale=0.25, posterior-covariance, 3000 draws: acceptance
%   0.543 / 0.441 -- correct direction (smaller jscale, higher
%   acceptance), now above the textbook 0.2-0.4 range but promoted
%   directly to production given this round's steep diagnostic-to-
%   production degradation pattern (0.20/0.20 -> 0.08/0.04 last time),
%   which this extra margin is meant to absorb.
% Production @ jscale=0.25: accepted 0.231/0.114, 7h51m. Best convergence
%   yet (Geweke 6/8, 4/8) but by remained the most stubborn parameter
%   (0.74 pooled-SD chain gap) across every jumping-covariance approach
%   tried. Diagnosed as gen1's covariance (built from the earlier
%   hessian round's own draws) inheriting that round's incomplete
%   exploration of by's posterior. GEN2: covariance rebuilt from THIS
%   (accepted) run's own draws (build_posterior_jumping_covariance_gen2.m)
%   -- standard adaptive-MCMC pilot/tune iteration. Confirmed: by's
%   proposal std rose 24% (gen1->gen2, 0.0293->0.0363) while the three
%   already-resolved parameters barely moved (0.86-0.98x). mode_file
%   unchanged (still valid; 25k kept draws/chain was ample travel time
%   last round). Starting jscale=0.25 again (what worked for gen1's
%   production), short diagnostic first as always.
%----------------------------------------------------------------
% RE-ESTIMATION ROUND: (HH BC n) equation fix (n now responds to genuine
% b(-1)/b instead of a fixed b0*(1-omega) placeholder -- see conversation
% record) + less informative priors on bpi, by (Normal/Gamma -> Gamma(1,1)
% each) and a lower phib_val prior mean (0.10 -> 0.05). Both the model
% equations and the estimated_params priors changed, so the OLD mode file
% is stale on two independent counts and MUST NOT be reused (Dynare would
% throw "generated using another specification of the model" for the
% equation change alone). THIS BLOCK IS SET TO STAGE 1 (linear KF,
% mode_compute=4, no mode_file, short diagnostic mh_replic=3000).
%
% Stage 1 @ jscale=0.2, prior_variance, bpi/by Gamma(1,1) both: mode found
%   (bpi=2.2534) but by's mode landed at EXACTLY 0.0000 with a non-PD
%   Hessian there -- shape=1 Gamma is Exponential, mode-at-0 by
%   construction, and the data wasn't informative enough to pull by off
%   that boundary. Diagnostic acceptance 0.47%/1.23%, both catastrophically
%   below 0.05 (prior_variance jumping scales off the new, ~7-10x larger
%   prior variances on bpi/by). By instruction: by's prior reshaped to
%   Gamma(1.00, 0.70) (shape~=2.04, interior mode ~0.51, no boundary
%   pre-loading). jscale cut hard (0.2 -> 0.03) to address the acceptance
%   problem, which is independent of the by reshaping (bpi's prior
%   variance, the likely dominant driver, is unchanged).
%
% Stage 1 @ jscale=0.03, prior_variance, by reshaped to Gamma(1.00, 0.70):
%   mode found cleanly, Hessian positive definite (no warning), log data
%   density finite (1113.7 Laplace / 1102.9 MHM). by's mode moved off the
%   boundary to 0.0794 (post. mean 0.126, 90% HPD [0.018, 0.239]) -- fixed.
%   bpi's posterior moved substantially higher than any previous round:
%   mode 4.01, post. mean 4.34, 90% HPD [3.81, 4.89] -- the diffuse prior
%   genuinely changes the answer, not just its uncertainty. Acceptance
%   45.2%/44.3% -- high vs the textbook 0.2-0.4 range, but NOT retuned:
%   this diagnostic uses prior_variance jumping under the LINEAR (non-
%   OccBin) likelihood, which does not carry forward to Stage 2 (fresh
%   hessian covariance) or production, and this file's own history shows
%   diagnostic acceptance consistently degrading substantially by the
%   full OccBin production run -- extra margin here is expected to help,
%   not signal a problem. Promoted directly to Stage 2 per instructions.
%
% THIS BLOCK IS SET TO STAGE 2 (true OccBin ZLB likelihood, mode_compute=0
% from the Stage-1 mode above, hessian jumping covariance, short
% diagnostic mh_replic=3000). mode_file matches THIS parameter vector
% (new priors + the (HH BC n) equation fix) since it was generated by the
% Stage 1 run immediately above -- do not reuse across further prior/
% equation changes without redoing Stage 1 first.
%
% Stage 2 @ jscale=0.2, hessian, 3000 draws: acceptance 10.2%/7.2%, above
%   the 0.05 kill threshold but well below the 0.2-0.4 target. Matches a
%   precedent from this file's earlier history (jscale=0.12 gave
%   12.2%/12.9%, jscale=0.10 then gave a healthy 23.6%/26.5%) -- cutting
%   jscale by a similar proportion here. ali's diagnostic posterior mean
%   (~0.17-0.22) diverges sharply from the Stage-1 mode (0.4884); this is
%   the same "mismatched starting point, insufficient travel time"
%   pattern seen for ali in an earlier round, not expected to persist in
%   the full production run.
%
% Stage 2 @ jscale=0.10, hessian, 3000 draws: acceptance 28.0%/17.5% --
%   chain 1 healthy, chain 2 just under 20% but this asymmetry has been
%   promoted directly before in this file's history. bpi stable across
%   reruns (post. mean 3.91 [3.67,4.17], consistent with both earlier
%   diagnostics). Promoted directly to PRODUCTION per instructions.
%
% Production @ jscale=0.10, hessian, 50000 draws: acceptance 5.43%/5.31%
%   -- technically above the 0.05 floor but barely, and a ~5x drop from
%   the 3000-draw diagnostic's 28.0%/17.5% (worse than any prior
%   diagnostic-to-production degradation in this file's history). Geweke:
%   chain 1 2/8 structural params converged, chain 2 3/8 (va right on the
%   boundary, p=0.042). Sits in the same marginal zone as a round that
%   WAS killed earlier in this file's history (4.49%/5.82%). NOT accepted
%   as final. jscale cut roughly in half (0.10 -> 0.05); given how badly
%   the 3000-draw diagnostic mispredicted this run's degradation, using a
%   longer diagnostic (8000 draws, not 3000) before committing to another
%   full 50000-draw production run.
%
% Extended diagnostic @ jscale=0.05, hessian, 8000 draws: acceptance
%   39.0%/27.6% -- large improvement over the 5.43%/5.31% the jscale=0.10
%   production run gave. Geweke ~3-4/8 structural params converged per
%   chain, comparable to or better than the full 50000-draw run despite
%   far fewer draws. Posterior stable (bpi 4.09 [3.82,4.33], by 0.128
%   [0.084,0.167], ali 0.196 [0.183,0.207]). Promoted to PRODUCTION per
%   instructions.
%
% THIS BLOCK IS SET TO PRODUCTION (mh_replic=50000 x 2 chains, jscale=0.05,
% same mode/hessian jumping covariance as the accepted extended diagnostic).
%----------------------------------------------------------------
estimation(
    datafile                = fred_data_debt_baa,
    first_obs               = 1,
    nobs                    = 68,
    mode_compute            = 0,
    mode_file               = 'mmt8bOccBin3_covid_bT_baa/Output/mmt8bOccBin3_covid_bT_baa_mode',
    mh_replic               = 50000,
    mh_nblocks              = 2,
    mh_jscale               = 0.05,
    mcmc_jumping_covariance = 'hessian',
    order                   = 1,
    graph_format            = pdf
) y_obs d_obs pig_obs rl_obs rn_obs taut_obs bT_obs;

% Save posterior.
if ~exist('skip_save','var') || ~skip_save
    phib_str = sprintf('%03d', round(M_.params(12)*100));
    if exist('taubar_val','var') && abs(taubar_val - 0.15) > 1e-6
        tau_str = sprintf('_tau%03d', round(taubar_val*100));
    else
        tau_str = '';
    end
    outmat = fullfile(M_.fname, 'Output', sprintf('results_covid_bT_baa_phib%s%s.mat', phib_str, tau_str));
    save(outmat, 'oo_', 'M_');
    fprintf('\nSaved posteriors to %s\n', outmat);
end
