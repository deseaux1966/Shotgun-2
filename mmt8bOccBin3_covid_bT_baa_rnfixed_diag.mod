% mmt8bOccBin3_covid_bT_baa_rnfixed.mod -- eps_rn_me FIXED (calibrated) at 0.0015 rather than
% estimated, following profile_rnsd.m's finding that this value (with the other policy/shock
% parameters re-optimized) recovers the true ZLB dates. Copy of mmt8bOccBin3_covid_bT_baa.mod;
% original files untouched.
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
b = (rlss/pigss)*b(-1) + gt - tau + b0*(1 - rlss/pigss) - (g - taubar);

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
var eps_rn_me;          stderr 0.0015;   % FIXED (calibrated), not estimated -- see rn_error_jscale_investigation_summary.md
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
 stderr eps_tau_me,         inv_gamma_pdf, 0.05, 0.025;
 stderr eps_bT_me,          inv_gamma_pdf, 0.05, 0.025;
 end;

estimated_params_init;
stderr epsa, 0.0172095;
stderr epsz, 0.0226426;
stderr epsva, 0.0391502;
stderr epsxa, 0.0362417;
stderr epsg, 0.0300289;
stderr epsr, 0.0256;
stderr epstau, 0.0306352;
stderr epsL, 0.0688822;
stderr eps_y_me, 0.0205987;
stderr eps_d_me, 0.0183915;
stderr eps_pig_me, 0.0135105;
stderr eps_rl_me, 0.0106643;
stderr eps_tau_me, 0.0196986;
stderr eps_bT_me, 0.0220863;
va, 2.889863;
eta, 1.1362798;
bpi, 3.646;
by, 0.0901;
ali, 0.1831;
br, 0.271;
rhoz, 0.4430581;
phib_val, 0.068925701;
end;

options_.occbin.likelihood.status = true;
options_.occbin.smoother.status   = false;
options_.occbin.likelihood.periodic_solution  = true;
options_.occbin.likelihood.use_updated_regime = 1;
options_.occbin.regularisation                = 1e-3;

% Short diagnostic MCMC (800 draws x 2 chains). The mode_compute=4 Hessian was singular/NaN
% (plausible given OccBin's nonsmooth likelihood), so this uses prior_variance jumping
% covariance as a fallback -- the same starting choice this project used before ever
% adopting Hessian-based jumping. mh_jscale is a first guess; retune based on acceptance.
estimation(
    datafile                = fred_data_debt_baa,
    first_obs               = 1,
    nobs                    = 68,
    mode_compute            = 0,
    mode_file               = 'mmt8bOccBin3_covid_bT_baa_rnfixed/Output/mmt8bOccBin3_covid_bT_baa_rnfixed_mode',
    mh_replic               = 800,
    mh_nblocks              = 2,
    mh_jscale               = 0.3,
    mcmc_jumping_covariance = 'prior_variance',
    order                   = 1,
    nograph
) y_obs d_obs pig_obs rl_obs rn_obs taut_obs bT_obs;
