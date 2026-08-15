% build_debt_dataset_fygfdpun.m
% Rebuild bT_obs using FYGFDPUN (Federal Debt Held by the Public) instead
% of GFDEBTN (Gross/Total Public Debt), and fix a data/model inconsistency
% found while writing up the estimation summary: the model's measurement
% equation is a LOG deviation, bT_obs = log(bT) - log(bTbar), but the
% previous data-construction scripts (build_debt_dataset.m,
% fetch_fred_debt_gdp.m) demeaned the raw ratio in LEVELS instead. This
% script demeans in logs, matching the model.
%
% FRED series:
%   FYGFDPUN - Federal Debt Held by the Public (quarterly, millions $)
%   GDP      - Gross Domestic Product (quarterly, nominal, billions $)
%   bT = (FYGFDPUN/1000) / GDP    [unit conversion: millions -> billions]
%
% Note: NGDPDQ, referenced in fetch_fred_debt_gdp.m's header comments and
% in mmt8bOccBin3_covid_bT_baa.mod's header, is not a valid FRED series
% (confirmed: fredgraph.csv?id=NGDPDQ returns an HTML error page, not
% data). The debt_to_gdp data actually used in all prior estimation runs
% came from grossdebt.csv (a pre-downloaded ratio, matching
% build_debt_dataset.m's input format), not from fetch_fred_debt_gdp.m.
%
% Sample: 2009Q1-2025Q4 (68 obs), matching the other observables.

clear; close all; clc;
shotgun_dir = 'C:/Users/gawater/Documents/GitHub/Shotgun 2';
cd(shotgun_dir);

%% --- Fetch FYGFDPUN and GDP directly (no API key needed) ---
fredR = @(id) ['https://fred.stlouisfed.org/graph/fredgraph.csv?id=' id];

fygfdpun = getFREDcsv(fredR('FYGFDPUN'));
gdp      = getFREDcsv(fredR('GDP'));

fprintf('FYGFDPUN: %d obs (%s to %s)\n', height(fygfdpun), ...
    string(fygfdpun.date(1)), string(fygfdpun.date(end)));
fprintf('GDP:      %d obs (%s to %s)\n', height(gdp), ...
    string(gdp.date(1)), string(gdp.date(end)));

%% --- Align to the estimation sample: 2009Q1 - 2025Q4 (68 quarters) ---
sample_start = datetime(2009,1,1);
sample_end   = datetime(2025,10,1);   % 2025Q4
qtrs = (sample_start:calquarters(1):sample_end)';
T = numel(qtrs);
fprintf('Target sample: %d quarters (%s to %s)\n', T, string(qtrs(1)), string(qtrs(end)));

fyg_q = alignQ(fygfdpun, qtrs);
gdp_q = alignQ(gdp,      qtrs);

n_nan = sum(isnan(fyg_q)) + sum(isnan(gdp_q));
if n_nan > 0
    error('Missing observations in sample window (%d NaNs) - check series coverage.', n_nan);
end

%% --- Compute ratio (unit fix: FYGFDPUN is millions $, GDP is billions $) ---
bT_raw = (fyg_q / 1000) ./ gdp_q;
fprintf('\nbT_raw (debt held by public / GDP): mean=%.4f min=%.4f max=%.4f\n', ...
    mean(bT_raw), min(bT_raw), max(bT_raw));

%% --- Demean in LOGS, matching the model's measurement equation ---
log_bT = log(bT_raw);
bT_obs = log_bT - mean(log_bT);
fprintf('bT_obs (log-demeaned): mean=%.10f std=%.6f min=%.4f max=%.4f\n', ...
    mean(bT_obs), std(bT_obs), min(bT_obs), max(bT_obs));

%% --- Combine with the other 6 observables (unchanged) ---
fprintf('\nLoading fred_data_debt_baa.mat for y,d,pig,rl(BAA),rn,taut...\n');
base = load('fred_data_debt_baa.mat');
required = {'y_obs','d_obs','pig_obs','rl_obs','rn_obs','taut_obs'};
for k = 1:numel(required)
    if ~isfield(base, required{k})
        error('Missing field "%s" in fred_data_debt_baa.mat.', required{k});
    end
end
n_orig = length(base.y_obs);
if n_orig ~= T
    error('Length mismatch: existing observables have %d obs, new bT has %d.', n_orig, T);
end

y_obs    = base.y_obs(:);
d_obs    = base.d_obs(:);
pig_obs  = base.pig_obs(:);
rl_obs   = base.rl_obs(:);
rn_obs   = base.rn_obs(:);
taut_obs = base.taut_obs(:);
bT_obs   = bT_obs(:);

%% --- Save (overwrite fred_data_debt_baa.mat in place) ---
save('fred_data_debt_baa.mat', 'y_obs','d_obs','pig_obs','rl_obs','rn_obs','taut_obs','bT_obs');
fprintf('\nSaved fred_data_debt_baa.mat (%d obs, bT_obs now from FYGFDPUN/GDP, log-demeaned)\n', T);

%% ========== Local functions ==========
function tbl = getFREDcsv(url)
    tmp = [tempdir 'fred_' num2str(randi(99999)) '.csv'];
    cmd = sprintf('curl -s -L --max-time 60 -o "%s" "%s"', tmp, url);
    [st, out] = system(cmd);
    if st ~= 0
        error('curl failed (exit %d): %s\nURL: %s', st, out, url);
    end
    T = readtable(tmp, 'VariableNamingRule','preserve');
    delete(tmp);
    vn = T.Properties.VariableNames;
    dates  = datetime(T.(vn{1}));
    values = T.(vn{2});
    if iscell(values) || ischar(values), values = str2double(values); end
    valid = ~isnan(values);
    tbl = table(dates(valid), values(valid), 'VariableNames', {'date','value'});
end

function v = alignQ(tbl, qtrs)
    v = NaN(numel(qtrs),1);
    for i = 1:numel(qtrs)
        idx = find(tbl.date == qtrs(i), 1);
        if ~isempty(idx)
            v(i) = tbl.value(idx);
        end
    end
end
