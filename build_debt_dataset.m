% build_debt_dataset.m
% Read grossdebt.csv, convert pct to decimal, extract 2009Q1-2025Q4,
% demean, and combine with fred_data_occbin to produce fred_data_debt.mat.

clear; close all; clc;

shotgun_dir = 'C:/Users/gawater/Documents/GitHub/Shotgun 2';
cd(shotgun_dir);

% --- Load debt-to-GDP data ---
fprintf('Reading grossdebt.csv...\n');
opts = detectImportOptions('grossdebt.csv');
opts = setvartype(opts, 'date', 'datetime');
opts = setvaropts(opts, 'date', 'InputFormat', 'MM/dd/yyyy');
debt_data = readtable('grossdebt.csv', opts);
fprintf('  Debt-to-GDP: %d observations (%s to %s)\n', ...
    height(debt_data), string(debt_data.date(1)), string(debt_data.date(end)));

% --- Convert percentage to decimal ---
debt_data.debt_to_gdp = debt_data.debt_to_gdp / 100;
fprintf('  Converted to decimal: mean = %.4f, std = %.4f\n', ...
    mean(debt_data.debt_to_gdp), std(debt_data.debt_to_gdp));

% --- Extract sample 2009Q1-2025Q4 ---
sample_start = datetime(2009, 1, 1);
sample_end   = datetime(2025, 10, 1);  % 2025Q4 = Oct 2025

idx = debt_data.date >= sample_start & debt_data.date <= sample_end;
bT_raw = debt_data.debt_to_gdp(idx);
n_debt = length(bT_raw);

fprintf('  Sample: %d observations (%s to %s)\n', n_debt, ...
    string(debt_data.date(idx(1))), string(debt_data.date(idx(end))));
fprintf('  Raw bT: mean = %.4f, std = %.4f\n', mean(bT_raw), std(bT_raw));

% Demean (same treatment as other observables)
bT_obs = bT_raw - mean(bT_raw);

% --- Load existing observables ---
fprintf('\nLoading fred_data_occbin.mat...\n');
existing = load('fred_data_occbin.mat');
y_obs    = existing.y_obs;
d_obs    = existing.d_obs;
pig_obs  = existing.pig_obs;
rl_obs   = existing.rl_obs;
rn_obs   = existing.rn_obs;
taut_obs = existing.taut_obs;
n_orig = length(y_obs);
fprintf('  Original: %d observations\n', n_orig);

% --- Match sample lengths ---
if n_orig ~= n_debt
    fprintf('WARNING: Length mismatch (%d vs %d). Truncating to %d.\n', ...
        n_orig, n_debt, min(n_orig, n_debt));
    n_use = min(n_orig, n_debt);
    y_obs    = y_obs(1:n_use);
    d_obs    = d_obs(1:n_use);
    pig_obs  = pig_obs(1:n_use);
    rl_obs   = rl_obs(1:n_use);
    rn_obs   = rn_obs(1:n_use);
    taut_obs = taut_obs(1:n_use);
    bT_obs   = bT_obs(1:n_use);
else
    n_use = n_debt;
end

% --- Assemble ---
new_data = table(...
    y_obs(1:n_use), d_obs(1:n_use), pig_obs(1:n_use), ...
    rl_obs(1:n_use), rn_obs(1:n_use), taut_obs(1:n_use), ...
    bT_obs(1:n_use), ...
    'VariableNames', {'y_obs','d_obs','pig_obs','rl_obs','rn_obs','taut_obs','bT_obs'});

save('fred_data_debt.mat', 'y_obs', 'd_obs', 'pig_obs', 'rl_obs', 'rn_obs', 'taut_obs', 'bT_obs');

fprintf('\nSaved fred_data_debt.mat (%d obs, 7 vars)\n', n_use);
fprintf('\nbT stats (demeaned):\n');
fprintf('  Mean:  %.10f\n', mean(bT_obs));
fprintf('  Std:   %.6f\n', std(bT_obs));
fprintf('  Min:   %.4f\n', min(bT_obs));
fprintf('  Max:   %.4f\n', max(bT_obs));
fprintf('\nDone.\n');
