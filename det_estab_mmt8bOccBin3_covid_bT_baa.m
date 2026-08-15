% det_estab_mmt8bOccBin3_covid_bT_baa.m
% Grid search over (bpi, by) for determinacy and E-stability in the
% mmt8bOccBin3_covid_bT_baa model (bT debt/GDP, BAA-based rl, epsL COVID
% labor-supply shock; true OccBin ZLB likelihood estimation), holding all
% other parameters at their posterior means from the production MCMC run.
%
% Unlike det_estab_mmt8bOccBin3.m (which swept 3 fixed/calibrated phib
% values), phib_val is now an estimated parameter with a single posterior
% mean -- so this script runs one grid at that mean, plus a sensitivity
% check at phib_val's 90% HPD bounds [0.0695, 0.1639] (flagged since the
% upper bound sits close to the ~0.153 destabilizing threshold found in
% the calibrated Section 4 analysis of the original model).
%
% bpi in [0, 5],  by in [0, 3]  (bpi range extended from [0,3] to include
% the current round's posterior mean, ~3.47, which the diffuse prior
% moved outside the old [0,3] grid)
%
% Determinacy : BK condition from resol on the RELAX regime (standard
%               linear model without the ZLB constraint).  info(1)==0
%               means exactly as many explosive eigenvalues as
%               forward-looking variables -> unique stable REE.
%
% E-stability : Spectral radius of the predetermined state-to-state
%               block of ghx < 1.  For linear RE models, determinacy and
%               E-stability of the MSV solution coincide (Evans &
%               Honkapohja 2001, Prop. 10.3; McCallum 2007).

clear; close all;
shotgun_dir = 'C:/Users/gawater/Documents/GitHub/Shotgun 2';
cd(shotgun_dir);
addpath('C:/dynare/7.1/matlab');

%% 1. Compile model (RELAX-regime linear solution used for BK check)
dynare mmt8bOccBin3_covid_bT_baa_bkonly noclearall nolog
cd(shotgun_dir);

%% 2. Load posterior means from the production MCMC run
rfile = fullfile(shotgun_dir, 'mmt8bOccBin3_covid_bT_baa', ...
    'mmt8bOccBin3_covid_bT_baa', 'Output', 'results_covid_bT_baa_phib5000.mat');
fprintf('Loading: %s\n', rfile);
R  = load(rfile);
pm = R.oo_.posterior_mean;

pnames = strtrim(cellstr(M_.param_names));
idx_bpi      = find(strcmp(pnames, 'bpi'),      1);
idx_by       = find(strcmp(pnames, 'by'),       1);
idx_phib_val = find(strcmp(pnames, 'phib_val'), 1);

param_fields = fieldnames(pm.parameters);
for k = 1:numel(param_fields)
    nm  = param_fields{k};
    idx = find(strcmp(pnames, nm), 1);
    if ~isempty(idx)
        M_.params(idx) = pm.parameters.(nm);
        fprintf('  %-10s posterior mean = %.4f\n', nm, pm.parameters.(nm));
    end
end

bpi_est      = pm.parameters.bpi;
by_est       = pm.parameters.by;
phib_val_est = pm.parameters.phib_val;

% phib_val 90% HPD bounds (from ESTIMATION RESULTS table; see the
% production log). Used only for the sensitivity check in Section 4 below.
phib_val_lo = 0.0301;
phib_val_hi = 0.1005;

%% 3. Grid search over (bpi, by) at the posterior-mean phib_val
N        = 60;
bpi_grid = linspace(0.01, 5.0, N);
by_grid  = linspace(0.00, 3.0, N);

graphdir = fullfile(shotgun_dir, 'mmt8bOccBin3_covid_bT_baa_bkonly', 'graphs');
if ~exist(graphdir, 'dir'); mkdir(graphdir); end

fprintf('\n============================================================\n');
fprintf(' Posterior-mean scenario: phib_val = %.4f\n', phib_val_est);
fprintf('============================================================\n');

det_map   = NaN(N, N);
estab_map = NaN(N, N);

for i = 1:N
    for j = 1:N
        M_.params(idx_bpi)      = bpi_grid(i);
        M_.params(idx_by)       = by_grid(j);
        M_.params(idx_phib_val) = phib_val_est;
        try
            [dr_tmp, info_tmp] = resol(0, M_, options_, oo_.dr, ...
                oo_.steady_state, oo_.exo_steady_state, ...
                oo_.exo_det_steady_state);
            code = info_tmp(1);
            det_map(i,j) = code;
            if code == 0
                Phi = dr_tmp.ghx(M_.nstatic+1 : M_.nstatic+M_.nspred, :);
                rho = max(abs(eig(Phi)));
                estab_map(i,j) = double(rho < 1);
            else
                estab_map(i,j) = 0;
            end
        catch
            det_map(i,j)   = 4;
            estab_map(i,j) = 0;
        end
    end
    if mod(i,10)==0, fprintf('  row %d / %d\n', i, N); end
end

M_.params(idx_bpi)      = bpi_est;
M_.params(idx_by)       = by_est;
M_.params(idx_phib_val) = phib_val_est;

n_det   = sum(det_map(:)==0);
n_estab = sum(estab_map(:)==1);
fprintf('  Determinate: %d / %d\n', n_det, N*N);
fprintf('  E-stable:    %d / %d\n', n_estab, N*N);

cdet = NaN(N,N);
cdet(det_map==0)                                 =  1;
cdet(det_map==3)                                 = -1;
cdet(det_map~=0 & det_map~=3 & ~isnan(det_map))  =  0;
cdet(isnan(det_map))                             =  0;

%% 4. Check whether the determinate and E-stable regions coincide exactly,
% then combine into a single two-region figure (as they are expected to,
% by Evans & Honkapohja 2001 Prop. 10.3, for this class of model).
region_det   = (det_map == 0);
region_estab = (estab_map == 1);
n_mismatch = sum(region_det(:) ~= region_estab(:));
fprintf('  Cells where determinate ~= E-stable: %d / %d\n', n_mismatch, N*N);

% combined = 1 where determinate AND E-stable, 0 otherwise (covers
% "indeterminate" and "explosive" together with "E-unstable" as one region).
combined = double(region_det & region_estab);

fig = figure('Units','inches','Position',[0.5 0.5 7 5.5],'Color','w');
imagesc(by_grid, bpi_grid, combined, [-0.5 1.5]);
set(gca, 'YDir','normal');
colormap(gca, [0.85 0.30 0.30; 0.30 0.70 0.40]);
cb = colorbar; cb.Ticks = [0 1];
cb.TickLabels = {'Indeterminate/explosive & E-unstable','Determinate & E-stable'};
hold on;
plot(by_est, bpi_est, 'k*', 'MarkerSize', 14, 'LineWidth', 2.5);
yline(1, 'k--', 'LineWidth', 1.4);
xlabel('$b_y$','Interpreter','latex','FontSize',15);
ylabel('$b_\pi$','Interpreter','latex','FontSize',15);
xlim([0 3]); ylim([0 5]);
text(by_est + 0.08, bpi_est, ...
    sprintf('Posterior mean $(%.2f,\\, %.2f)$', bpi_est, by_est), ...
    'FontSize', 9, 'VerticalAlignment', 'middle', ...
    'HorizontalAlignment', 'left', 'BackgroundColor', 'none', ...
    'Interpreter', 'latex');

outfile = fullfile(shotgun_dir, 'det_estab_map_mmt8bOccBin3_covid_bT_baa.pdf');
if exist('exportgraphics','file') || exist('exportgraphics','builtin')
    try exportgraphics(fig, outfile, 'ContentType','vector','BackgroundColor','white');
    catch, print(fig, outfile, '-dpdf', '-painters'); end
else
    print(fig, outfile, '-dpdf', '-painters');
end
fprintf('  Saved: %s\n', outfile);

% Report at estimated means
[~, i_e] = min(abs(bpi_grid - bpi_est));
[~, j_e] = min(abs(by_grid  - by_est));
fprintf('  At posterior mean: BK=%d (0=det), E-stable=%d\n', ...
    det_map(i_e,j_e), estab_map(i_e,j_e));

close all;

%% 7. Sensitivity check: phib_val at its 90% HPD bounds, (bpi,by) at posterior means
fprintf('\n============================================================\n');
fprintf(' Sensitivity: phib_val at 90%% HPD bounds [%.4f, %.4f]\n', phib_val_lo, phib_val_hi);
fprintf(' (bpi, by) held at posterior means (%.4f, %.4f)\n', bpi_est, by_est);
fprintf('============================================================\n');
M_.params(idx_bpi) = bpi_est;
M_.params(idx_by)  = by_est;
for phib_test = [phib_val_lo, phib_val_est, phib_val_hi]
    M_.params(idx_phib_val) = phib_test;
    try
        [dr_tmp, info_tmp] = resol(0, M_, options_, oo_.dr, ...
            oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);
        code = info_tmp(1);
        if code == 0
            Phi = dr_tmp.ghx(M_.nstatic+1 : M_.nstatic+M_.nspred, :);
            rho = max(abs(eig(Phi)));
            estab = rho < 1;
        else
            rho = NaN; estab = false;
        end
        fprintf('  phib_val = %.4f : BK=%d (0=det), spectral radius=%.4f, E-stable=%d\n', ...
            phib_test, code, rho, estab);
    catch ME
        fprintf('  phib_val = %.4f : FAILED (%s)\n', phib_test, ME.message);
    end
end
M_.params(idx_phib_val) = phib_val_est;

fprintf('\nDone.\n');
