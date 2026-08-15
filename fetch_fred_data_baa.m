% fetch_fred_data_baa.m
% Builds fred_data_debt_baa.mat: identical to fred_data_debt.mat except
% rl_obs is constructed from the BAA corporate bond yield instead of MPRIME.
%
% BAA (Moody's Seasoned Baa Corporate Bond Yield, FRED, monthly, %).
%   Advantages over MPRIME:
%     - Has genuine credit-spread variation (widens in recessions, narrows
%       in recoveries) rather than being locked at FF + 300bp by convention.
%     - Average over 2009-2025 (~5.1%) aligns closer to the model SS rl
%       than MPRIME average (~4.3%).
%   Caveats:
%     - Corporate bond rate: includes term premium and liquidity effects
%       absent from a short-term bank lending rate.
%     - More volatile than a bank loan rate.
%
% rl_obs = log(1 + BAA_quarterly/400) - mean_sample[log(1 + BAA_quarterly/400)]
%        = demeaned log gross quarterly BAA rate.
%
% All other 7 observables (y, d, pig, rn, taut, bT, g) are taken unchanged
% from fred_data_debt.mat.
%
% REQUIRED FILE (one browser click):
%   BAA.csv : https://fred.stlouisfed.org/graph/fredgraph.csv?id=BAA
%   Save as BAA.csv in:  C:\Users\gawater\Documents\GitHub\Shotgun\
%
% Sample: 2009Q1 - 2025Q4 (nobs = 68, matching fred_data_debt.mat).

cd('C:\Users\gawater\Documents\GitHub\Shotgun 2');

%% ---- Step 1: Load base debt data (7 unchanged observables) ---------------
fprintf('=== Step 1: Loading fred_data_debt.mat ===\n');
base = load('fred_data_debt.mat');
required = {'y_obs','d_obs','pig_obs','rl_obs','rn_obs','taut_obs','bT_obs'};
for k = 1:numel(required)
    if ~isfield(base, required{k})
        error('Missing field "%s" in fred_data_debt.mat.', required{k});
    end
end
T = length(base.y_obs);
sample_dates = datetime(2009,1,1) + calquarters(0:T-1)';
fprintf('T = %d  (%s - %s)\n', T, ...
    char(sample_dates(1),'yyyy-QQQ'), char(sample_dates(end),'yyyy-QQQ'));

%% ---- Step 2: Load BAA.csv ----------------------------------------------
fprintf('\n=== Step 2: Loading BAA.csv ===\n');
if ~isfile('BAA.csv')
    error(['BAA.csv not found.\n' ...
           'Download from:\n' ...
           '  https://fred.stlouisfed.org/graph/fredgraph.csv?id=BAA\n' ...
           'Save as BAA.csv in:  %s'], pwd);
end
baa_tbl = read_fred_csv('BAA.csv');
fprintf('BAA: %d rows  (%s to %s)\n', height(baa_tbl), ...
    char(baa_tbl.Date(1),'yyyy-MM-dd'), char(baa_tbl.Date(end),'yyyy-MM-dd'));

%% ---- Step 3: Quarterly average -----------------------------------------
fprintf('\n=== Step 3: Quarterly averaging ===\n');
baa_q = quarterly_avg(baa_tbl, sample_dates);   % Tx1, annualised %

n_nan = sum(isnan(baa_q));
if n_nan > 0
    warning('%d quarter(s) have no BAA data - check CSV coverage.', n_nan);
end
fprintf('BAA quarterly: mean = %.3f%%,  range = [%.3f, %.3f]\n', ...
    mean(baa_q,'omitnan'), min(baa_q,[],'omitnan'), max(baa_q,[],'omitnan'));

%% ---- Step 4: rl_obs from BAA -------------------------------------------
fprintf('\n=== Step 4: Computing rl_obs (BAA) ===\n');

% Gross quarterly rate: r_q = 1 + r_ann/400
baa_q_gross = 1 + baa_q / 400;
log_baa     = log(baa_q_gross);
rl_obs_baa  = log_baa - mean(log_baa, 'omitnan');   % demeaned log-deviation

fprintf('rl_obs (BAA):   range = [%.6f, %.6f],  std = %.6f\n', ...
    min(rl_obs_baa), max(rl_obs_baa), std(rl_obs_baa,'omitnan'));
fprintf('rl_obs (MPRIME, original): std = %.6f\n', std(base.rl_obs(:),'omitnan'));

% Correlation with original MPRIME-based rl_obs
fprintf('Corr(BAA rl_obs, MPRIME rl_obs): %.4f\n', ...
    corr(rl_obs_baa(:), base.rl_obs(:), 'rows','pairwise'));

%% ---- Step 5: Save ------------------------------------------------------
fprintf('\n=== Step 5: Saving fred_data_debt_baa.mat ===\n');
y_obs    = base.y_obs(:);
d_obs    = base.d_obs(:);
pig_obs  = base.pig_obs(:);
rl_obs   = rl_obs_baa(:);          % <-- BAA replaces MPRIME here
rn_obs   = base.rn_obs(:);
taut_obs = base.taut_obs(:);
bT_obs   = base.bT_obs(:);
save('fred_data_debt_baa.mat', ...
    'y_obs','d_obs','pig_obs','rl_obs','rn_obs','taut_obs','bT_obs');
fprintf('Saved fred_data_debt_baa.mat  (%d obs, 7 series)\n', T);

%% =========================================================================
function tbl = read_fred_csv(filepath)
    txt   = fileread(filepath);
    txt   = strrep(txt, char(13), '');
    lines = strsplit(txt, newline);
    lines = strtrim(lines);
    lines = lines(~cellfun(@isempty, lines));
    if contains(lower(lines{1}), 'date')
        lines = lines(2:end);
    end
    n = numel(lines);
    dates  = NaT(n,1);
    values = NaN(n,1);
    for k = 1:n
        p = strsplit(lines{k}, ',');
        if numel(p) >= 2
            try
                dates(k) = datetime(strtrim(p{1}), 'InputFormat','yyyy-MM-dd');
            catch; continue; end
            v = strtrim(p{2});
            if ~isempty(v) && ~strcmp(v,'.')
                values(k) = str2double(v);
            end
        end
    end
    ok  = ~isnat(dates);
    tbl = table(dates(ok), values(ok), 'VariableNames', {'Date','Value'});
end

function q_avg = quarterly_avg(tbl, sample_dates)
    T     = length(sample_dates);
    q_avg = NaN(T,1);
    for t = 1:T
        q_start = sample_dates(t);
        q_end   = q_start + calquarters(1) - caldays(1);
        mask    = tbl.Date >= q_start & tbl.Date <= q_end & ~isnan(tbl.Value);
        if any(mask)
            q_avg(t) = mean(tbl.Value(mask));
        end
    end
end
