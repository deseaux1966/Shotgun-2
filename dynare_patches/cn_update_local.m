function U = cn_update_local(U, rho)
Z.init = randn(size(U.init));
Z.eps  = randn(size(U.eps));
Z.resample = randn(size(U.resample));
s = sqrt(1-rho^2);
U.init = rho*U.init + s*Z.init;
U.eps  = rho*U.eps  + s*Z.eps;
U.resample = rho*U.resample + s*Z.resample;
end
