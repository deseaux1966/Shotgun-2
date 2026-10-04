% det_estab_bpi_phib.m
% Determinacy / E-stability map over (bpi, phib_val): bpi in [0,5],
% phib_val in [0,1], all other parameters (including by) at the posterior
% means of the production round with sigma_me,rn fixed at 0.0015
% (mmt8bOccBin3_covid_bT_baa_rnfixed_prod). Same machinery as
% det_estab_mmt8bOccBin3_covid_bT_baa.m: RELAX-regime linear model
% (mmt8bOccBin3_covid_bT_baa_bkonly.mod), BK determinacy from resol(), and
% E-stability as spectral radius of the predetermined-state block of ghx < 1.
% resol() re-solves the steady state at each grid point (phib_val enters the
% tax rule, so b_ss and tau_ss move with it).

addpath('C:/dynare/7.2/matlab');
dynare mmt8bOccBin3_covid_bT_baa_bkonly noclearall nolog

% options_.qz_criterium is left empty by the .mod (no stoch_simul/estimation
% option sets it); dyn_first_order_solver then crashes (line 216) whenever the
% BK count is wrong instead of returning info=3/4. Use Dynare's default.
if isempty(options_.qz_criterium), options_.qz_criterium = 1 + 1e-6; end

R  = load(fullfile('mmt8bOccBin3_covid_bT_baa_rnfixed_prod','Output', ...
    'mmt8bOccBin3_covid_bT_baa_rnfixed_prod_results.mat'));
pm = R.oo_.posterior_mean.parameters;
pnames = strtrim(cellstr(M_.param_names));
for nm = fieldnames(pm)'
    idx = find(strcmp(pnames, nm{1}), 1);
    if ~isempty(idx), M_.params(idx) = pm.(nm{1}); fprintf('  %-9s = %.6f\n', nm{1}, pm.(nm{1})); end
end
i_bpi  = find(strcmp(pnames,'bpi'),1);
i_phib = find(strcmp(pnames,'phib_val'),1);
bpi_est = pm.bpi; phib_est = pm.phib_val; by_est = pm.by;

bpi_grid  = linspace(0,5,101);
phib_grid = linspace(0,1,51);
Nb = numel(bpi_grid); Np = numel(phib_grid);
code_map  = NaN(Np,Nb);   % rows: phib, cols: bpi
rho_map   = NaN(Np,Nb);   % spectral radius of predetermined block (determinate points)
estab_map = zeros(Np,Nb);

for i = 1:Np
    for j = 1:Nb
        M_.params(i_bpi)  = bpi_grid(j);
        M_.params(i_phib) = phib_grid(i);
        try
            [dr_tmp, info_tmp] = resol(0, M_, options_, oo_.dr, oo_.steady_state, ...
                oo_.exo_steady_state, oo_.exo_det_steady_state);
            code_map(i,j) = info_tmp(1);
            if info_tmp(1) == 0
                Phi = dr_tmp.ghx(M_.nstatic+1 : M_.nstatic+M_.nspred, :);
                rho_map(i,j)   = max(abs(eig(Phi)));
                estab_map(i,j) = double(rho_map(i,j) < 1);
            end
        catch
            code_map(i,j) = -99;
        end
    end
    if mod(i,10)==0, fprintf('  phib row %d / %d\n', i, Np); end
end
M_.params(i_bpi) = bpi_est; M_.params(i_phib) = phib_est;

u = unique(code_map(~isnan(code_map)))';
for c = u, fprintf('  info code %3d : %d points\n', c, sum(code_map(:)==c)); end
fprintf('  determinate: %d, determinate & E-stable: %d, of %d\n', ...
    sum(code_map(:)==0), sum(estab_map(:)==1), Np*Nb);
save('det_estab_bpi_phib_results.mat','code_map','rho_map','estab_map','bpi_grid','phib_grid','bpi_est','phib_est','by_est','pm');
