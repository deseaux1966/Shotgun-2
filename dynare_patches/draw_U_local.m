function U = draw_U_local(state_variance_rank, n_innov, N, T)
U.init = randn(state_variance_rank, N);
U.eps  = randn(n_innov, N, T);
U.resample = randn(1, T);
end
