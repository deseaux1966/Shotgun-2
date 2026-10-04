% mmt8bOccBin3_covid_bT_baa_softzlb.mod -- DIAGNOSTIC ONLY.
% Smooth (logistic) approximation of the ZLB floor, replacing OccBin's
% piecewise-linear occbin_constraints entirely, to test whether a
% nonlinear-perturbation + particle-filter estimation is a viable
% alternative to OccBin for this model. Copy of
% mmt8bOccBin3_covid_bT_baa_rnfixed_prod.mod; original files untouched.
%
% MOTIVATION: OccBin's guess-and-verify regime-consistency search has been
% the source of repeated likelihood-evaluation failures in this project
% (smoother non-convergence at the posterior mean, false ZLB-binding flags
% in 2017-2019, mh_jscale/eps_rn_me pathologies). A smooth policy rule has
% no kink, hence no regime search and no discrete jump in the likelihood
% surface -- see rn_error_jscale_investigation_summary.md and Section 2.5
% of mmt8bOccBin3_covid_bT_baa_epsL_summary.tex for the OccBin-side history
% this is meant to be an alternative to.
%
% SMOOTHING FUNCTIONAL FORM: rn = rn_lb + (rn_notional-rn_lb) / (1+exp(-kappa_zlb*(rn_notional-rn_lb)))
% This is a logistic soft-max of (rn_notional, rn_lb): as kappa_zlb -> infinity
% it converges to max(rn_notional, rn_lb), i.e. the exact OccBin floor.
% kappa_zlb=2000 below is NOT taken from a specific paper -- the two
% candidate references discussed (Andreasen & Kronborg 2022 QE 13,
% 1171-1202 "extended perturbation" w/ a smooth Taylor rule; Richter &
% Throckmorton 2015 BEJM 15(1), 157-182 on smoothing for global solution
% convergence) were only available to me as scanned/protected PDFs, so I
% could not extract their exact functional form or calibrated smoothing
% parameter. kappa_zlb is therefore my own construction, chosen so the
% logistic transition width (~10/kappa_zlb in gross-quarterly-rate units)
% is about 0.005 (~2 pp annualized) -- a free choice to sensitivity-check,
% not a verified literature value.
%
% PLAN: (1) steady+check the smooth model at current posterior-mean
% parameters (no ZLB kink so this is now a single, non-regime-switching
% determinacy question); (2) confirm a 2nd-order pruned perturbation
% solution exists; (3) evaluate (not estimate -- mh_replic=0) the
% particle-filter likelihood once at the posterior mean to see whether it
% returns a finite value without error, before committing to a full
% particle-filter MCMC run.
%
% Posterior-mean parameter values below are read directly from
% mmt8bOccBin3_covid_bT_baa_rnfixed_prod/Output/..._results.mat
% (oo_.posterior_mean.parameters / .shocks_std), i.e. the actual production
% OccBin estimation this whole investigation is about, NOT the Stage-1
% starting values in mmt8bOccBin3_covid_bT_baa_rnfixed_prod.mod's own
% estimated_params_init block (those are the chunk-0 starting point, not
% the final mean).

var y c d n pig rl rd h hg hd hn hs mu w rn rn_notional m at zt vat xat gt tau b bT zr
    y_obs d_obs pig_obs rl_obs rn_obs taut_obs bT_obs;
varexo epsa epsz epsva epsxa epsg epsr epstau epsL;

% NOTE: measurement error is declared below via the standard Dynare
% "var <obsname>; stderr ...;" mechanism (which populates M_.H) rather than
% as explicit eps_*_me exogenous shocks inside the model equations (as the
% OccBin production file does). initial_estimation_checks.m:116-120
% requires particle-filter estimation to have measurement error registered
% in M_.H specifically -- an eps_*_me shock loading through Sigma_e is
% invisible to that check even though it is econometrically equivalent
% additive iid noise on the observable. This is a mechanical reformulation
% for particle-filter compatibility, not a substantive change to the model.

Parameters chi psi va beta eta theta nu xa xn xh phin phip
            g taubar phib_val bpi by br ali omega b0 z rlss pigss rn_lb kappa_zlb
            rhoa rhoz rhova rhoxa rhog rhotau;

ali    =  0.184101;
beta   =  0.9967;
eta    =  1.187837;
va     =  3.068990;
xa     =  5.2528;

chi    =  0.5;
psi    =  1.3;
omega  =  0.0;
b0     =  1.0;
z      =  1.3;

g      =  0.2;
taubar =  0.14;
phib_val = 0.091204;

pigss  =  1.02^(0.25);
rlss   =  1.02^(0.25) / beta;
rn_lb  =  1.0;
kappa_zlb = 2000;

rhoa   =  0.95;
rhoz   =  0.455941;
rhova  =  0.0;
rhoxa  =  0.95;
rhog   =  0.0;
rhotau =  0.5;

bpi    =  3.355493;
by     =  0.082551;
br     =  0.252438;

theta  =  6.0;
nu     =  0.25;
xn     =  0.95;
xh     = 1 - xn;
phin   =  0.0005;
phip   =  50;

model;
% --- Household / bank block ---
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
bT = b / (1 + g);

% --- Monetary policy: notional (unconstrained) Taylor rule, unchanged ---
rn_notional = (STEADY_STATE(rl)^ali) * (((pig)/STEADY_STATE(pig))^bpi)
              * ((y/steady_state(y))^by) * exp(zr);
zr = br * zr(-1) + epsr;

% --- Smooth (logistic) ZLB floor: replaces OccBin's occbin_constraints.
% As kappa_zlb -> infinity this equation -> rn = max(rn_notional, rn_lb).
rn = rn_lb + (rn_notional - rn_lb) / (1 + exp(-kappa_zlb*(rn_notional - rn_lb)));

% --- Exogenous shock processes ---
at = (at(-1)^rhoa) * exp(epsa);
zt = (zt(-1)^rhoz) * (z^(1-rhoz)) * exp(epsz);
vat = (vat(-1)^rhova) * (va^(1-rhova)) * exp(epsva);
xat = (xat(-1)^rhoxa) * (xa^(1-rhoxa)) * exp(epsxa);
gt = (rhog*gt(-1)) + ((1-rhog)*g) + epsg;

% --- Observable equations (7 observables). Measurement error enters via
% M_.H (declared in the shocks block below), not as an explicit shock here.
y_obs    = log(y)            - log(STEADY_STATE(y));
d_obs    = log(d)            - log(STEADY_STATE(d));
pig_obs  = log(pig)          - log(STEADY_STATE(pig));
rl_obs   = log(rl)           - log(STEADY_STATE(rl));
rn_obs   = log(rn)           - log(STEADY_STATE(rn));
taut_obs = log(tau)          - log(STEADY_STATE(tau));
bT_obs   = log(bT)           - log(STEADY_STATE(bT));
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
var epsa;               stderr 0.018665;
var epsz;               stderr 0.024515;
var epsva;              stderr 0.039137;
var epsxa;              stderr 0.037542;
var epsg;               stderr 0.029707;
var epsr;               stderr 0.021563;
var epstau;             stderr 0.029190;
var epsL;               stderr 0.058791;
% Measurement error (M_.H), matching the posterior-mean eps_*_me stderrs
% from the production OccBin run (eps_rn_me fixed at 0.0015 there too).
var y_obs;               stderr 0.022137;
var d_obs;               stderr 0.019422;
var pig_obs;             stderr 0.012380;
var rl_obs;              stderr 0.010414;
var rn_obs;              stderr 0.0015;
var taut_obs;            stderr 0.020180;
var bT_obs;              stderr 0.025222;
end;

% --- Step 1: steady state + BK check at posterior-mean parameters ---
steady(solve_algo=4);
check;

% --- Step 2: confirm a pruned 2nd-order perturbation solution exists ---
stoch_simul(order=2, pruning, noprint, nograph, nocorr, nofunctions);

varobs y_obs d_obs pig_obs rl_obs rn_obs taut_obs bT_obs;

% --- Step 3: single particle-filter likelihood evaluation at the
% posterior mean (mode_compute=0, no mode_file => xparam1 = the
% estimated_params_init values below; mh_replic=0 => no MCMC; cova_compute=0
% => skip the (expensive, noisy-for-a-particle-filter) numerical Hessian).
estimated_params;
va,         gamma_pdf,  2.50,  0.75;
eta,        gamma_pdf,  1.53,  0.40;
bpi,        gamma_pdf,  1.00,  1.00;
by,         gamma_pdf,  1.00,  0.70;
ali,        beta_pdf,   0.49,  0.10;
br,         beta_pdf,   0.36,  0.10;
rhoz,       beta_pdf,   0.45,  0.10;
phib_val,   beta_pdf,   0.05,  0.05;

stderr epsa,               inv_gamma_pdf, 0.05, 0.025;
stderr epsz,               inv_gamma_pdf, 0.05, 0.025;
stderr epsva,              inv_gamma_pdf, 0.05, 0.025;
stderr epsxa,              inv_gamma_pdf, 0.05, 0.025;
stderr epsg,               inv_gamma_pdf, 0.05, 0.025;
stderr epsr,               inv_gamma_pdf, 0.05, 0.025;
stderr epstau,             inv_gamma_pdf, 0.05, 0.025;
stderr epsL,               inv_gamma_pdf, 0.15, 0.075;
end;

estimated_params_init;
stderr epsa, 0.018665;
stderr epsz, 0.024515;
stderr epsva, 0.039137;
stderr epsxa, 0.037542;
stderr epsg, 0.029707;
stderr epsr, 0.021563;
stderr epstau, 0.029190;
stderr epsL, 0.058791;
va, 3.068990;
eta, 1.187837;
bpi, 3.355493;
by, 0.082551;
ali, 0.184101;
br, 0.252438;
rhoz, 0.455941;
phib_val, 0.091204;
end;

options_.particle.status = true;
options_.particle.filter_algorithm = 'sis';
options_.particle.number_of_particles = 5000;

estimation(
    datafile        = fred_data_debt_baa,
    first_obs       = 1,
    nobs            = 68,
    mode_compute    = 0,
    mh_replic       = 0,
    cova_compute    = 0,
    order           = 2,
    pruning,
    nograph
) y_obs d_obs pig_obs rl_obs rn_obs taut_obs bT_obs;
