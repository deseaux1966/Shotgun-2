% check_zlb_sensitivity.m
% Sensitivity check: how much of each shock's contribution to pig (and bT)
% over 2021Q2-2023Q2 comes from innovations realized in the ZLB quarters
% (2020Q1-2021Q1, t=45:49) vs. earlier vs. inside the window itself.
% Contributions are linear in the smoothed innovations, so each shock's
% path is simulated three times, with only one innovation-date block active.

if ~exist('oo_','var') || ~isfield(oo_,'shock_decomposition')
    addpath('C:\dynare\7.1\matlab');
    dynare mmt8bOccBin3_covid_bT_baa_decomp.mod
end

exo = cellstr(M_.exo_names);
nT  = size(oo_.shock_decomposition,3);
nx  = numel(exo);
E = zeros(nT, nx);
for k = 1:nx
    E(:,k) = oo_.SmoothedShocks.(exo{k})(:);
end

idx_pig = strmatch('pig', cellstr(M_.endo_names), 'exact');
idx_bT  = strmatch('bT',  cellstr(M_.endo_names), 'exact');
window  = 50:58;
blocks  = struct('name',{'pre-2020 (t<45)','ZLB qtrs 2020Q1-2021Q1 (t=45:49)','in-window 2021Q2+ (t>=50)'}, ...
                 'idx',{1:44, 45:49, 50:nT});

ys = oo_.dr.ys;
shocks_show = {'epsa','epsz','epsva','epsxa','epsg','epsr','epstau','epsL'};

for vv = {'pig','bT'}
    iv = strmatch(vv{1}, cellstr(M_.endo_names), 'exact');
    fprintf('\n=== %s: avg contribution over 2021Q2-2023Q2, by origin of innovation ===\n', vv{1});
    fprintf('%-8s %12s %12s %12s %12s %12s\n','shock','pre-2020','ZLB qtrs','in-window','sum','shockdecomp');
    for s = 1:numel(shocks_show)
        k = strmatch(shocks_show{s}, exo, 'exact');
        vals = zeros(1,3);
        for b = 1:3
            ex_ = zeros(nT, nx);
            ex_(blocks(b).idx, k) = E(blocks(b).idx, k);
            y = simult_(M_, options_, ys, oo_.dr, ex_, 1);
            dev = y(iv, 2:end)' - ys(iv);
            vals(b) = mean(dev(window));
        end
        ref = mean(squeeze(oo_.shock_decomposition(iv,k,window)));
        fprintf('%-8s %12.5f %12.5f %12.5f %12.5f %12.5f\n', shocks_show{s}, vals, sum(vals), ref);
    end
end
