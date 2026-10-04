% det_estab_bpi_by.m
% Determinacy / E-stability map over (bpi, by): bpi in [0,5], by in [0,3],
% all other parameters (incl. phib_val) at the posterior means of the
% production round with sigma_me,rn fixed at 0.0015. Same machinery as
% det_estab_bpi_phib.m (RELAX-regime linear model, BK via resol, E-stability =
% spectral radius of the predetermined-state block of ghx < 1). Also reports
% the spectral radius at the posterior mean and at the phib_val 90% HPD bounds.

addpath('C:/dynare/7.2/matlab');
dynare mmt8bOccBin3_covid_bT_baa_bkonly noclearall nolog
if isempty(options_.qz_criterium), options_.qz_criterium = 1 + 1e-6; end

R  = load(fullfile('mmt8bOccBin3_covid_bT_baa_rnfixed_prod','Output', ...
    'mmt8bOccBin3_covid_bT_baa_rnfixed_prod_results.mat'));
pm = R.oo_.posterior_mean.parameters;
pnames = strtrim(cellstr(M_.param_names));
for nm = fieldnames(pm)'
    idx = find(strcmp(pnames, nm{1}), 1);
    if ~isempty(idx), M_.params(idx) = pm.(nm{1}); end
end
i_bpi = find(strcmp(pnames,'bpi'),1); i_by = find(strcmp(pnames,'by'),1);
i_phib = find(strcmp(pnames,'phib_val'),1);
bpi_est = pm.bpi; by_est = pm.by; phib_est = pm.phib_val;

bpi_grid = linspace(0,5,101);
by_grid  = linspace(0,3,61);
Nb = numel(bpi_grid); Ny = numel(by_grid);
code_map = NaN(Ny,Nb); rho_map = NaN(Ny,Nb); estab_map = zeros(Ny,Nb);
for i = 1:Ny
    for j = 1:Nb
        M_.params(i_bpi) = bpi_grid(j); M_.params(i_by) = by_grid(i);
        try
            [dr_tmp, info_tmp] = resol(0, M_, options_, oo_.dr, oo_.steady_state, ...
                oo_.exo_steady_state, oo_.exo_det_steady_state);
            code_map(i,j) = info_tmp(1);
            if info_tmp(1) == 0
                Phi = dr_tmp.ghx(M_.nstatic+1 : M_.nstatic+M_.nspred, :);
                rho_map(i,j) = max(abs(eig(Phi))); estab_map(i,j) = double(rho_map(i,j) < 1);
            end
        catch
            code_map(i,j) = -99;
        end
    end
    if mod(i,10)==0, fprintf('  by row %d / %d\n', i, Ny); end
end
M_.params(i_bpi) = bpi_est; M_.params(i_by) = by_est;

u = unique(code_map(~isnan(code_map)))';
for c = u, fprintf('  info code %3d : %d points\n', c, sum(code_map(:)==c)); end
fprintf('  determinate: %d, det & E-stable: %d, of %d; mismatches det vs E-stable: %d\n', ...
    sum(code_map(:)==0), sum(estab_map(:)==1), Ny*Nb, sum((code_map(:)==0) ~= (estab_map(:)==1)));

% threshold in bpi along by = by_est (bisection) and spectral radii at phib points
f = @(b) resolcode(M_,options_,oo_,i_bpi,b);
lo = 0.5; hi = 2;
for it = 1:30, mid = (lo+hi)/2; if f(mid)==0, hi = mid; else lo = mid; end, end
bpi_thresh = hi; fprintf('  bpi threshold at by=%.4f: %.5f\n', by_est, bpi_thresh);
phib_pts = [R.oo_.posterior_hpdinf.parameters.phib_val, phib_est, R.oo_.posterior_hpdsup.parameters.phib_val];
rho_phib = zeros(1,3);
for k = 1:3
    M_.params(i_phib) = phib_pts(k);
    [dr_tmp, info_tmp] = resol(0, M_, options_, oo_.dr, oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);
    Phi = dr_tmp.ghx(M_.nstatic+1 : M_.nstatic+M_.nspred, :); rho_phib(k) = max(abs(eig(Phi)));
    fprintf('  phib=%.4f: info=%d spectral radius=%.4f\n', phib_pts(k), info_tmp(1), rho_phib(k));
end
M_.params(i_phib) = phib_est;
save('det_estab_bpi_by_results.mat','code_map','rho_map','estab_map','bpi_grid','by_grid','bpi_est','by_est','phib_est','bpi_thresh','phib_pts','rho_phib');

function c = resolcode(M_,options_,oo_,ib,b)
M_.params(ib) = b;
[~,info] = resol(0,M_,options_,oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state);
c = info(1);
end
