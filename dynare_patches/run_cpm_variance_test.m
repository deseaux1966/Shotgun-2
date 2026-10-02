% run_cpm_variance_test.m -- correlated-pseudo-marginal variance diagnostic.
% Requires ReducedForm, Y, start, options_, M_, state_variance_rank, n_innov, T
% already in the base workspace (built by the companion setup script run
% interactively before this). Not a function file (needs the caller's
% workspace), but keeps its own local functions below (R2016b+ script-with-
% local-functions support) to avoid MATLAB's "no function definitions in
% evaluated-code context" restriction.

Nvals = [1000 5000];
rhovals = [0 0.9 0.99 0.999];
nrep = 10;

opts = options_.particle;
threads = options_.threads;

fprintf('\n%8s %8s %10s %10s %10s\n','N','rho','mean(diff)','var(diff)','var(single)');
for Ni = 1:numel(Nvals)
    N = Nvals(Ni);
    opts.number_of_particles = N;
    for ri = 1:numel(rhovals)
        rho = rhovals(ri);
        diffs = zeros(nrep,1);
        singles = zeros(nrep,1);
        for r=1:nrep
            U1 = draw_U_local(state_variance_rank, n_innov, N, T);
            L1 = sequential_importance_particle_filter_cpm(ReducedForm, Y, start, opts, threads, options_, M_, U1);
            U2 = cn_update_local(U1, rho);
            L2 = sequential_importance_particle_filter_cpm(ReducedForm, Y, start, opts, threads, options_, M_, U2);
            diffs(r) = L2 - L1;
            singles(r) = L1;
        end
        fprintf('%8d %8.3f %10.4f %10.4f %10.4f\n', N, rho, mean(diffs), var(diffs), var(singles));
    end
end

function U = draw_U_local(state_variance_rank, n_innov, N, T)
    U.init = randn(state_variance_rank, N);
    U.eps  = randn(n_innov, N, T);
    U.resample = randn(1, T);
end

function U = cn_update_local(U, rho)
    Z.init = randn(size(U.init));
    Z.eps  = randn(size(U.eps));
    Z.resample = randn(size(U.resample));
    s = sqrt(1-rho^2);
    U.init = rho*U.init + s*Z.init;
    U.eps  = rho*U.eps  + s*Z.eps;
    U.resample = rho*U.resample + s*Z.resample;
end
