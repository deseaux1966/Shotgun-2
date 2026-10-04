% mmt8bOccBin3_covid_bT_baa_softzlb_kappasweep.mod -- DIAGNOSTIC ONLY.
% Steady-state-only sweep over kappa_zlb (the logistic ZLB-smoothing
% sharpness parameter introduced in mmt8bOccBin3_covid_bT_baa_softzlb.mod).
% No particle filter, no estimation -- just steady(); check(); at OccBin's
% posterior-mode parameters, for each kappa_zlb value, to see (a) how much
% the steady-state "shaving" artifact (rn_ss vs rn_notional_ss) shrinks as
% kappa sharpens, and (b) whether the Blanchard-Kahn order/rank conditions
% remain verified throughout. kappa_zlb is set from the MATLAB workspace
% variable kappa_zlb_val if present (driven by a wrapper script that loops
% over values), following this project's existing taubar_val convention;
% otherwise defaults to 2000. Copy of mmt8bOccBin3_covid_bT_baa_softzlb.mod;
% original files untouched.

var y c d n pig rl rd h hg hd hn hs mu w rn rn_notional m at zt vat xat gt tau b bT zr
    y_obs d_obs pig_obs rl_obs rn_obs taut_obs bT_obs;
varexo epsa epsz epsva epsxa epsg epsr epstau epsL;

Parameters chi psi va beta eta theta nu xa xn xh phin phip
            g taubar phib_val bpi by br ali omega b0 z rlss pigss rn_lb kappa_zlb
            rhoa rhoz rhova rhoxa rhog rhotau;

% Parameters at OccBin's posterior MODE (same point used in _softzlb_mode.mod).
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
if exist('kappa_zlb_val','var')
    kappa_zlb = kappa_zlb_val;
else
    kappa_zlb = 2000;
end;

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

% --- Smooth (logistic) ZLB floor: replaces OccBin's occbin_constraints.
rn = rn_lb + (rn_notional - rn_lb) / (1 + exp(-kappa_zlb*(rn_notional - rn_lb)));

at = (at(-1)^rhoa) * exp(epsa);
zt = (zt(-1)^rhoz) * (z^(1-rhoz)) * exp(epsz);
vat = (vat(-1)^rhova) * (va^(1-rhova)) * exp(epsva);
xat = (xat(-1)^rhoxa) * (xa^(1-rhoxa)) * exp(epsxa);
gt = (rhog*gt(-1)) + ((1-rhog)*g) + epsg;

y_obs    = log(y)            - log(STEADY_STATE(y));
d_obs    = log(d)            - log(STEADY_STATE(d));
pig_obs  = log(pig)          - log(STEADY_STATE(pig));
rl_obs   = log(rl)           - log(STEADY_STATE(rl));
rn_obs   = log(rn)           - log(STEADY_STATE(rn));
taut_obs = log(tau)          - log(STEADY_STATE(tau));
bT_obs   = log(bT)           - log(STEADY_STATE(bT));
end;

varobs y_obs d_obs pig_obs rl_obs rn_obs taut_obs bT_obs;

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
