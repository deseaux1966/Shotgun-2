% check_zlb_sensitivity_v2.m
% Re-split of each shock's contribution to pig / bT / rn over sub-windows of
% 2021Q2-2023Q2, by ORIGIN of the innovation, with ZLB windows consistent
% with Section 2.5 (verified regime path binds 2009Q1-2016Q4 and 2020Q2-2022Q2).
% t=1 is 2009Q1: 2020Q1=45, 2020Q2=46, 2021Q1=49, 2021Q2=50, 2022Q2=54,
% 2022Q3=55, 2023Q2=58.
% Origin blocks (nest the old three-block table: old in-window = blocks 3+4):
%   1: t<=44            pre-2020
%   2: t=45:49          2020Q1-2021Q1 (old "ZLB quarters")
%   3: t=50:54          2021Q2-2022Q2 (in window, floor binding per regime path)
%   4: t=55:nT          2022Q3+       (post-liftoff)
% Requires oo_, M_, options_ in workspace from a *_decomp*.mod run.

exo = cellstr(M_.exo_names);
nT  = size(oo_.shock_decomposition,3);
nx  = numel(exo);
E = zeros(nT, nx);
for k = 1:nx
    E(:,k) = oo_.SmoothedShocks.(exo{k})(:);
end
ys = oo_.dr.ys;
blk = {1:44, 45:49, 50:54, 55:nT};
blkname = {'pre-2020','2020Q1-21Q1','21Q2-22Q2(bind)','22Q3+(post)'};
windows = {50:58, 50:54, 55:58};
winname = {'2021Q2-2023Q2','2021Q2-2022Q2 (binding)','2022Q3-2023Q2 (post-liftoff)'};
shocks_show = {'epsa','epsz','epsva','epsxa','epsg','epsr','epstau','epsL'};

results = struct();
for vv = {'pig','bT','rn'}
    iv = strmatch(vv{1}, cellstr(M_.endo_names), 'exact');
    C = zeros(nT, numel(shocks_show), 4);
    for s = 1:numel(shocks_show)
        k = strmatch(shocks_show{s}, exo, 'exact');
        for b = 1:4
            ex_ = zeros(nT, nx);
            ex_(blk{b}, k) = E(blk{b}, k);
            y = simult_(M_, options_, ys, oo_.dr, ex_, 1);
            C(:,s,b) = y(iv, 2:end)' - ys(iv);
        end
    end
    results.(vv{1}).C = C;
    for w = 1:numel(windows)
        fprintf('\n=== %s: avg contribution over %s, by origin of innovation ===\n', vv{1}, winname{w});
        fprintf('%-8s %12s %12s %16s %12s %12s %12s\n','shock',blkname{:},'sum','shockdecomp');
        for s = 1:numel(shocks_show)
            k = strmatch(shocks_show{s}, exo, 'exact');
            vals = squeeze(mean(C(windows{w},s,:),1))';
            ref = mean(squeeze(oo_.shock_decomposition(iv,k,windows{w})));
            fprintf('%-8s %12.5f %12.5f %16.5f %12.5f %12.5f %12.5f\n', shocks_show{s}, vals, sum(vals), ref);
        end
        tot = mean(squeeze(oo_.shock_decomposition(iv,17,windows{w})));
        fprintf('actual smoothed dev from SS (avg): %.5f\n', tot);
    end
end

fprintf('\n=== smoothed epsr by sub-window ===\n');
kr = strmatch('epsr', exo, 'exact');
for w = 1:numel(windows)
    fprintf('%-30s mean epsr = %.5f\n', winname{w}, mean(E(windows{w},kr)));
end
iv = strmatch('rn', cellstr(M_.endo_names), 'exact');
fprintf('\n=== smoothed rn level (gross quarterly) by sub-window ===\n');
for w = 1:numel(windows)
    fprintf('%-30s mean rn = %.5f\n', winname{w}, ys(iv) + mean(squeeze(oo_.shock_decomposition(iv,17,windows{w}))));
end
save('zlb_sens_v2_last.mat','results','blkname','winname','shocks_show');
