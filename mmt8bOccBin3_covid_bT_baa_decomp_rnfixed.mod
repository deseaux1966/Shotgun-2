% mmt8bOccBin3_covid_bT_baa_decomp_rnfixed.mod -- same as _decomp.mod but at the POSTERIOR MEAN of the
% rnfixed production round (eps_rn_me fixed at 0.0015; values from
% mmt8bOccBin3_covid_bT_baa_rnfixed_prod/Output/*_results.mat, oo_.posterior_mean). Original _decomp.mod untouched.
% mmt8bOccBin3_covid_bT_baa_decomp.mod
% Historical shock decomposition at the posterior MEAN, to identify which
% structural shocks drove the post-COVID inflation spike (pig_obs, roughly
% 2021Q2-2023Q2 per ShotgunPCestCovid.pdf discussion of "elevated inflation").
%
% Uses calib_smoother (Kalman smoother at fixed, posterior-mean parameters,
% no MCMC) + shock_decomposition on the actual estimation dataset
% (fred_data_debt_baa.mat), same posterior means as
% mmt8bOccBin3_covid_bT_baa_irf.mod.
%
% NOTE: this uses the LINEAR (non-OccBin) Kalman smoother -- same
% approximation as this codebase's own "Stage 1" mode-finding step, not the
% full OccBin piecewise-linear ZLB smoother used for the accepted
% estimation's likelihood. CAUTION: the target window (2021Q2-2023Q2) is NOT
% away from the binding ZLB -- the verified regime path (tex Section 2.5) has
% the floor binding through 2022Q2, so 2021Q2-2022Q2 are binding quarters in
% which this linear smoother is infeasible (rn<1). See check_zlb_sensitivity_v2.m.
%
% Data indexing: t=1 is 2009Q1 (per the estimation .mod header: 2020Q1=t=45,
% 2020Q2=t=46). So 2021Q2 = t=50, 2023Q2 = t=58.

var y c d n pig rl rd h hg hd hn hs mu w rn rn_notional m at zt vat xat gt tau b bT zr
    y_obs d_obs pig_obs rl_obs rn_obs taut_obs bT_obs;
varexo epsa epsz epsva epsxa epsg epsr epstau epsL
        eps_y_me eps_d_me eps_pig_me eps_rl_me eps_rn_me eps_tau_me eps_bT_me;

Parameters chi psi va beta eta theta nu xa xn xh phin phip
            g taubar phib_val bpi by br ali omega b0 z rlss pigss rn_lb
            rhoa rhoz rhova rhoxa rhog rhotau;

% --- Calibrated (fixed) parameters, unchanged from the estimation .mod ---
beta   =  0.9967;
chi    =  0.5;
psi    =  1.3;
omega  =  0.0;
b0     =  1.0;
z      =  1.3;

g      =  0.2;
taubar =  0.14;

pigss  =  1.02^(0.25);
rlss   =  1.02^(0.25) / beta;
rn_lb  =  1.0;

rhoa   =  0.95;
rhova  =  0.0;
rhoxa  =  0.95;
rhog   =  0.0;
rhotau =  0.5;

theta  =  6.0;
nu     =  0.25;
xn     =  0.95;
xh     = 1 - xn;
phin   =  0.0005;
phip   =  50;

% --- Estimated parameters, set to POSTERIOR MEAN (Table 5 / *_mean.mat) ---
va       =  3.068990;
eta      =  1.187837;
bpi      =  3.355493;
by       =  0.082551;
ali      =  0.184101;
br       =  0.252438;
rhoz     =  0.455941;
phib_val =  0.091204;
xa       =  5.2528;

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

% --- Monetary policy (RELAX regime only -- linear smoother, no ZLB) ---
rn_notional = (STEADY_STATE(rl)^ali) * (((pig)/STEADY_STATE(pig))^bpi)
              * ((y/steady_state(y))^by) * exp(zr);
zr = br * zr(-1) + epsr;

rn = rn_notional;

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

% --- Shock standard deviations, set to POSTERIOR MEAN ---
shocks;
var epsa;               stderr 0.018665;
var epsz;               stderr 0.024515;
var epsva;              stderr 0.039137;
var epsxa;              stderr 0.037542;
var epsg;               stderr 0.029707;
var epsr;               stderr 0.021563;
var epstau;             stderr 0.029190;
var epsL;               stderr 0.058791;
var eps_y_me;           stderr 0.022137;
var eps_d_me;           stderr 0.019422;
var eps_pig_me;         stderr 0.012380;
var eps_rl_me;          stderr 0.010414;
var eps_rn_me;          stderr 0.0015;
var eps_tau_me;         stderr 0.020180;
var eps_bT_me;          stderr 0.025222;
end;

steady(solve_algo=4);
check;

varobs y_obs d_obs pig_obs rl_obs rn_obs taut_obs bT_obs;

calib_smoother(datafile=fred_data_debt_baa, first_obs=1, nobs=68) y_obs d_obs pig_obs rl_obs rn_obs taut_obs bT_obs pig;

shock_decomposition(datafile=fred_data_debt_baa, first_obs=1, nobs=68, nograph) pig_obs pig;
