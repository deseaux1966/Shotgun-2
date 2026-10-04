% mmt8bOccBin3_covid_bT_baa_bkonly.mod
% Stripped-down version of mmt8bOccBin3_covid_bT_baa.mod for BK/E-stability
% analysis, following the same convention as mmt8bOccBin3_bkonly.mod:
% no estimation, no OccBin constraints -- just steady + check on the
% RELAX-regime linear model (rn = rn_notional unconditionally, no ZLB
% floor). det_estab_mmt8bOccBin3_covid_bT_baa.m calls this to get M_,
% oo_, options_ without triggering MCMC estimation.
%
% All shocks (including epsL, the COVID labor-supply shock) are kept,
% matching mmt8bOccBin3_bkonly.mod's convention: purely contemporaneous
% iid varexo shocks with no own law of motion do not enter the
% endogenous-variable Jacobian, so they cannot affect the BK eigenvalue
% count or the E-stability spectral radius. Stripping OccBin's
% occbin_constraints/BIND-regime equation is what actually matters for
% isolating "the non-ZLB equilibrium" -- not which shocks are declared.

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

% chi recalibrated from 1.0 to 0.5, matching the production .mod.
chi    =  0.5;
psi    =  1.3;
omega  =  0.0;
b0     =  1.0;
z      =  1.3;

g      =  0.2;
% Recalibrated from 0.20 to 0.14 to match the production .mod's new
% deficit target (g-taubar=0.06, 6% of GDP). This default is overwritten
% by det_estab_mmt8bOccBin3_covid_bT_baa.m with the posterior mean value
% loaded from the results file anyway.
taubar =  0.14;
phib_val = 0.10;

pigss  =  1.02^(0.25);
rlss   =  1.02^(0.25) / beta;
rn_lb  =  1.0;

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
hs = (1/chi) * (vat*c/d)^chi;
(c/z)*(h^psi) + chi*hs = (1/eta)*exp(epsL);
(rl-rd)*(h^psi)*d = chi*hs*w;
mu = at * w * (h^psi);
beta*(rl(+1)/pig(+1))*(mu(+1)/mu) = 1;

y = zt*hg;

d = xat*((xn^(1/nu))*(n^((nu-1)/nu)) + (xh^(1/nu))*(hd^((nu-1)/nu)))^(nu/(nu-1));
hn = phin*n;
n  = ((rl-rd)^nu)*(xat^(nu-1))*xn*d*(rl-rn+phin*z)^(-nu);
hd = ((rl-rd)^nu)*(xat^(nu-1))*xh*d*(z^(-nu));
h  = hg + hd + hn;

0.5*(pig-STEADY_STATE(pig))*pig/(STEADY_STATE(pig)^2) = (1-theta+theta*(w/zt))/phip
    + beta*(pig(+1)-STEADY_STATE(pig))*pig(+1)*(y(+1)/y)/(STEADY_STATE(pig)^2);

y = c + gt + 0.5*phip*(((pig-STEADY_STATE(pig))/(STEADY_STATE(pig)))^2);

m = (gt) - (tau) + (rn)*n;
n = (m(-1)/pig) + (b(-1)/pig) - (b/rl);
tau = (tau(-1)^rhotau) * (taubar^(1-rhotau)) * ((b(-1)/((1-omega)*b0))^phib_val) * exp(epstau);
b = (rlss/pigss)*b(-1) + gt - tau + (1-omega)*b0*(1 - rlss/pigss) - (g - taubar);

bT = b / (1 + g);

rn_notional = (STEADY_STATE(rl)^ali) * (((pig)/STEADY_STATE(pig))^bpi)
              * ((y/steady_state(y))^by) * exp(zr);
zr = br * zr(-1) + epsr;

% RELAX regime only: rn = rn_notional (no ZLB binding, no occbin_constraints)
rn = rn_notional;

at   = (at(-1)^rhoa)    * exp(epsa);
zt   = (zt(-1)^rhoz)    * (z^(1-rhoz))   * exp(epsz);
vat  = (vat(-1)^rhova)  * (va^(1-rhova)) * exp(epsva);
xat  = (xat(-1)^rhoxa)  * (xa^(1-rhoxa)) * exp(epsxa);
gt   = (rhog*gt(-1))    + ((1-rhog)*g)   + epsg;

y_obs    = log(y)   - log(STEADY_STATE(y))   + eps_y_me;
d_obs    = log(d)   - log(STEADY_STATE(d))   + eps_d_me;
pig_obs  = log(pig) - log(STEADY_STATE(pig)) + eps_pig_me;
rl_obs   = log(rl)  - log(STEADY_STATE(rl))  + eps_rl_me;
rn_obs   = log(rn)  - log(STEADY_STATE(rn))  + eps_rn_me;
taut_obs = log(tau) - log(STEADY_STATE(tau)) + eps_tau_me;
bT_obs   = log(bT)  - log(STEADY_STATE(bT))  + eps_bT_me;
end;

initval;
pig  = 1.004963;
rl   = 1.008270;
rn   = 1.004128;
rn_notional = 1.004128;
zr   = 0;
w    = 1.08333;
at   = 1;  zt = 1.3;  vat = 2.5035;  xat = 5.2528;  gt = 0.2;  tau = 0.20;  b = 1.0;
n    = 3.927;
m    = 3.943;
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
var epsa;   stderr 0.1;
var epsz;   stderr 0.1;
var epsva;  stderr 0.1;
var epsxa;  stderr 0.1;
var epsg;   stderr 0.1;
var epsr;   stderr 0.1;
var epstau; stderr 0.1;
var epsL;   stderr 0.1;
var eps_y_me;   stderr 0.1;
var eps_d_me;   stderr 0.1;
var eps_pig_me; stderr 0.1;
var eps_rl_me;  stderr 0.1;
var eps_rn_me;  stderr 0.1;
var eps_tau_me; stderr 0.1;
var eps_bT_me;  stderr 0.1;
end;

steady(solve_algo=4);
check;

% Silent first-order solve so oo_.dr.order_var (and ghx) are fully
% populated before external MATLAB scripts call resol().
stoch_simul(order=1, noprint, nograph, nocorr, nofunctions);
