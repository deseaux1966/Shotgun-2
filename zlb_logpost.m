function [lp, ll, info] = zlb_logpost(x, rh, mode_flag)
% Log posterior of the OccBin model with a GIVEN regime history rh (fixed,
% not inferred), using a time-varying-matrix Kalman filter built from
% occbin.check_regimes. Needs Dynare globals from a prior dynare run.
global M_ options_ oo_ bayestopt_ estim_params_ dataset_
info = 0; ll = -Inf; lp = -Inf;
x = x(:);
if any(x(1:15) <= 0), return, end
M = set_all_parameters(x, estim_params_, M_);
Q = M.Sigma_e; H = M.H; if isempty(H), H = 0; end
opt = options_;
try
    [T,R,SS,info,dr,M.params,TTx,RRx,CCx,~,~] = occbin.dynare_resolve(M,opt,oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state,[],'restrict');
catch
    info = 1; return
end
if info(1), return, end
opts_regime.regime_history = rh; opts_regime.binding_indicator = [];
[TT,RR,CC] = occbin.check_regimes(TTx,RRx,CCx,opts_regime,M,opt,dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state);
Y = dataset_.data'; % nvar x nobs
mf = bayestopt_.mf; nobs = size(Y,2); m = size(T,1);
Z = zeros(numel(mf), m); for i=1:numel(mf), Z(i,mf(i)) = 1; end
a = zeros(m,1);
P = lyapunov_symm(T, R*Q*R', options_.lyapunov_fixed_point_tol, options_.lyapunov_complex_threshold, 1);
ll = 0;
for t = 1:nobs
    v = Y(:,t) - Z*a; F = Z*P*Z' + H;
    ll = ll - 0.5*(log(det(F)) + v'*(F\v) + numel(v)*log(2*pi));
    K = P*Z'/F; ah = a + K*v; Ph = P - K*Z*P;
    Tt = TT(:,:,t+1); Rt = RR(:,:,t+1); Ct = CC(:,t+1);
    a = Tt*ah + Ct; P = Tt*Ph*Tt' + Rt*Q*Rt';
end
lpr = priordens(x, bayestopt_.pshape, bayestopt_.p6, bayestopt_.p7, bayestopt_.p3, bayestopt_.p4);
lp = ll + lpr;
end
