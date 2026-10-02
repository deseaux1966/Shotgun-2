% profile_rnsd05.m -- re-optimise (bpi, by, stderr epsr, ali, br) with stderr eps_rn_me FIXED at 0.0005.
% Other parameters at the posterior means. Then re-run the OccBin filter/smoother at the optimum.
restoredefaultpath; addpath('C:\dynare\7.2\matlab'); cd('C:\Users\gawater\Documents\GitHub\Shotgun-2');
evalc('dynare rn05_prof.mod noclearall');
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info
if isempty(options_.qz_criterium), options_.qz_criterium = 1 + 1e-6; end
bounds = prior_bounds(bayestopt_, options_.prior_trunc);
nm = cellstr(bayestopt_.name); ix = @(n) find(strcmp(nm, n));
x0 = get_all_parameters(estim_params_, M_); x0 = x0(:);
f = @(x) dsge_likelihood(x, dataset_, dataset_info, options_, M_, estim_params_, bayestopt_, bounds, oo_.dr, oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);
lp = @(x) -f(x);                       % dsge_likelihood returns minus log posterior
fprintf('start (sd_rn=0.0005, others at posterior mean): logpost %.2f\n', lp(x0));
free = [ix('bpi') ix('by') ix('stderr epsr') ix('ali') ix('br')];
obj = @(z) penal(-lp(setfree(x0, free, z)));
z = x0(free); opts = optimset('MaxFunEvals', 500, 'MaxIter', 500, 'TolX', 1e-5, 'TolFun', 1e-3, 'Display', 'off');
for pass = 1:3
    [z, fv] = fminsearch(obj, z, opts);
    fprintf('pass %d: logpost %.2f  bpi=%.3f by=%.4f sd_epsr=%.4f ali=%.4f br=%.3f\n', pass, -fv, z(1), z(2), z(3), z(4), z(5));
end
xb = setfree(x0, free, z);
fprintf('\nOptimised values -> baseline comparison: (sd_rn 0.0029, posterior mean) logpost 1266.56\n');
save('profile_rnsd05_results.mat','xb','z');

% OccBin filter/smoother regime path at the optimum
tmpl = fileread('rn05_template.txt');
rep = {'stderr epsr, 0.023622;', sprintf('stderr epsr, %.7f;', xb(ix('stderr epsr')));
       'bpi, 3.47416;', sprintf('bpi, %.6f;', xb(ix('bpi')));
       'by, 0.0832422;', sprintf('by, %.7f;', xb(ix('by')));
       'ali, 0.180836;', sprintf('ali, %.7f;', xb(ix('ali')));
       'br, 0.262394;', sprintf('br, %.6f;', xb(ix('br')))};
for k = 1:size(rep,1), assert(contains(tmpl, rep{k,1})); tmpl = strrep(tmpl, rep{k,1}, rep{k,2}); end
fid = fopen('rn05_run.mod','w'); fwrite(fid,tmpl); fclose(fid);
evalc('dynare rn05_run.mod noclearall');
T = 68; q = 0:T-1; yr = 2009+floor(q/4); qq = mod(q,4)+1; ql = arrayfun(@(y,k) sprintf('%dQ%d',y,k),yr,qq,'UniformOutput',false);
sm = oo_.occbin.smoother; h = sm.realtime_regime_history; b = false(T,1); for t=1:T, b(t) = h(t).regime(1)==1; end
d = diff([0;b;0]); st = find(d==1); en = find(d==-1)-1; sp = '';
for i=1:numel(st), sp = [sp sprintf(' %s-%s', ql{st(i)}, ql{en(i)})]; end
fprintf('regime path at optimum (warning_flag %d): %s\n', sm.warning_flag(1), sp);
fprintf('flagged share 2017-19: %.2f ; data-implied ZLB quarters flagged: %d of 41\n', mean(b(yr'>=2017&yr'<=2019)), sum(b(1:32))+sum(b(46:54)));

function v = penal(v), if ~isfinite(v) || v > 1e5, v = 1e5; end, end
function x = setfree(x, idx, z), x(idx) = z(:); end
