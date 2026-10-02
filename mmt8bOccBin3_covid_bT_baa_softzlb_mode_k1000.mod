% mmt8bOccBin3_covid_bT_baa_softzlb_mode.mod -- DIAGNOSTIC ONLY (Option 2).
% Mode search under the smooth-ZLB + particle-filter likelihood, for direct
% comparison against OccBin's own posterior mode. Copy of
% mmt8bOccBin3_covid_bT_baa_softzlb.mod (see that file for the smoothing
% rationale and the H-matrix measurement-error reformulation); original
% files untouched.
%
% WHY A MODE SEARCH SPECIFICALLY: non_linear_dsge_likelihood.m's particle
% filters (sequential_importance_particle_filter.m:72 etc.) hard-reset
% MATLAB's RNG to a fixed seed=0 on every call
% (set_dynare_seed_local_options([],false,'default')), confirmed directly
% by re-evaluating the likelihood at the same parameters under three
% different set_dynare_seed() calls and getting bit-identical results to
% 12 decimals. This "common random numbers" behavior makes the likelihood
% a smooth, reproducible function of parameters -- good for an optimizer,
% but NOT valid for MH/MCMC in the pseudo-marginal sense (Andrieu & Roberts
% 2009 requires fresh, independent likelihood-estimator noise at every
% iteration). So: use this property for what it's actually good for (mode
% finding), not for an MCMC chain.
%
% STARTING VALUES: OccBin's own posterior MODE (not mean), read directly
% from mmt8bOccBin3_covid_bT_baa_rnfixed_prod/Output/..._rnfixed_prod_mh_mode.mat
% (fval=1291.297, i.e. log posterior ~1291.30 under OccBin/Kalman). Starting
% the smooth-model optimizer from OccBin's own answer makes "how far does it
% move" the direct object of interest, rather than conflating a different
% starting point's basin of attraction with a genuine difference in what
% the smoothed model prefers. eps_*_me measurement-error entries from that
% same mode vector are used as the (calibrated, not re-estimated) M_.H
% values here, for consistency with the diagnostic file's H-matrix approach.
%
% PARTICLE COUNT REDUCED TO 1000 (from 5000 in the diagnostic file) purely
% for mode-search speed -- a single evaluation is ~5x cheaper, and csminwel
% needs on the order of hundreds of evaluations (numerical gradient at each
% iteration). The CRN property means this is still a well-defined,
% reproducible fixed target to optimize, not added noise -- but it is a
% cruder (1000-particle) approximation of the likelihood than the 5000-
% particle diagnostic evaluation, so the resulting mode/fval is a fast
% first pass, not a final production number.

var y c d n pig rl rd h hg hd hn hs mu w rn rn_notional m at zt vat xat gt tau b bT zr
    y_obs d_obs pig_obs rl_obs rn_obs taut_obs bT_obs;
varexo epsa epsz epsva epsxa epsg epsr epstau epsL;

Parameters chi psi va beta eta theta nu xa xn xh phin phip
            g taubar phib_val bpi by br ali omega b0 z rlss pigss rn_lb kappa_zlb
            rhoa rhoz rhova rhoxa rhog rhotau;

% Parameters initialized at OccBin's posterior MODE (not mean).
ali    =  0.180891;
beta   =  0.9967;
eta    =  1.259586;
va     =  2.746560;
xa     =  5.2528;

chi    =  0.5;
psi    =  1.3;
omega  =  0.0;
b0     =  1.0;
z      =  1.3;

g      =  0.2;
taubar =  0.14;
phib_val = 0.093842;

pigss  =  1.02^(0.25);
rlss   =  1.02^(0.25) / beta;
rn_lb  =  1.0;
kappa_zlb = 1000;

rhoa   =  0.95;
rhoz   =  0.502569;
rhova  =  0.0;
rhoxa  =  0.95;
rhog   =  0.0;
rhotau =  0.5;

bpi    =  3.279499;
by     =  0.072119;
br     =  0.255415;

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
b = (rlss/pigss)*b(-1) + gt - tau + b0*(1 - rlss/pigss) - (g - taubar);

% --- Debt-to-GDP ratio ---
bT = b / (1 + g);

% --- Monetary policy: notional (unconstrained) Taylor rule, unchanged ---
rn_notional = (STEADY_STATE(rl)^ali) * (((pig)/STEADY_STATE(pig))^bpi)
              * ((y/steady_state(y))^by) * exp(zr);
zr = br * zr(-1) + epsr;

% --- Smooth (logistic) ZLB floor: replaces OccBin's occbin_constraints.
rn = rn_lb + (rn_notional - rn_lb) / (1 + exp(-kappa_zlb*(rn_notional - rn_lb)));

% --- Exogenous shock processes ---
at = (at(-1)^rhoa) * exp(epsa);
zt = (zt(-1)^rhoz) * (z^(1-rhoz)) * exp(epsz);
vat = (vat(-1)^rhova) * (va^(1-rhova)) * exp(epsva);
xat = (xat(-1)^rhoxa) * (xa^(1-rhoxa)) * exp(epsxa);
gt = (rhog*gt(-1)) + ((1-rhog)*g) + epsg;

% --- Observable equations. Measurement error enters via M_.H. ---
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
var epsa;               stderr 0.019501;
var epsz;               stderr 0.024089;
var epsva;              stderr 0.046312;
var epsxa;              stderr 0.042058;
var epsg;               stderr 0.029470;
var epsr;               stderr 0.017955;
var epstau;             stderr 0.031459;
var epsL;               stderr 0.063620;
% Measurement error (M_.H), set at OccBin's posterior-MODE eps_*_me values
% (calibrated here, not re-estimated -- see softzlb.mod for why M_.H is
% used instead of explicit shocks). eps_rn_me fixed at 0.0015 at OccBin's
% mode too (it was calibrated, not estimated, there).
var y_obs;               stderr 0.022753;
var d_obs;               stderr 0.016978;
var pig_obs;             stderr 0.011867;
var rl_obs;              stderr 0.010884;
var rn_obs;              stderr 0.0015;
var taut_obs;            stderr 0.019902;
var bT_obs;              stderr 0.022693;
end;

steady(solve_algo=4);
check;

varobs y_obs d_obs pig_obs rl_obs rn_obs taut_obs bT_obs;

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

% Starting values below are OccBin's own posterior MODE (same starting
% point used for the original kappa_zlb=2000, 1000-particle cold-start run,
% for a like-for-like comparison of how kappa alone shifts the result).
estimated_params_init;
stderr epsa, 0.019501;
stderr epsz, 0.024089;
stderr epsva, 0.046312;
stderr epsxa, 0.042058;
stderr epsg, 0.029470;
stderr epsr, 0.017955;
stderr epstau, 0.031459;
stderr epsL, 0.063620;
va, 2.746560;
eta, 1.259586;
bpi, 3.279499;
by, 0.072119;
ali, 0.180891;
br, 0.255415;
rhoz, 0.502569;
phib_val, 0.093842;
end;

options_.particle.status = true;
options_.particle.filter_algorithm = 'sis';
options_.particle.number_of_particles = 1000;

% mode_compute=4 (csminwel, gradient-based) triggers an INTERACTIVE prompt
% here -- Dynare itself warns "Particle filtering introduces discontinuities
% in the objective function w.r.t the parameters. Thus, should use a
% non-gradient based optimizer" and blocks for a keyboard choice, which
% hangs a non-interactive/scripted run. mode_compute=8 (Dynare's own
% Nelder-Mead simplex) is chosen directly here to avoid that prompt, per
% Dynare's own recommendation printed at that prompt.
estimation(
    datafile        = fred_data_debt_baa,
    first_obs       = 1,
    nobs            = 68,
    mode_compute    = 8,
    mh_replic       = 0,
    cova_compute    = 0,
    order           = 2,
    pruning,
    nograph
) y_obs d_obs pig_obs rl_obs rn_obs taut_obs bT_obs;
