% profile_rnsd.m -- for fixed stderr eps_rn_me in SDS, re-optimise (bpi, by, stderr epsr, ali, br),
% then evaluate the OccBin regime path at the optimum (backing off toward the start if the
% Dynare run at the exact optimum hits the edge of the feasible region).
restoredefaultpath; addpath('C:\dynare\7.2\matlab'); cd('C:\Users\gawater\Documents\GitHub\Shotgun-2');
SDS = [0.001 0.0015];
pm = struct('epsr',0.023622,'bpi',3.47416,'by',0.0832422,'ali',0.180836,'br',0.262394);
T = 68; q = 0:T-1; yr = 2009+floor(q/4); qq = mod(q,4)+1; ql = arrayfun(@(y,k) sprintf('%dQ%d',y,k),yr,qq,'UniformOutput',false);
tp = fileread('rnsd_prof_tmpl.txt'); tr = fileread('rnsd_run_tmpl.txt');
fill = @(t, sd, v) strrep(strrep(strrep(strrep(strrep(strrep(t,'RNSD',sprintf('%.14g',sd)),'EPSR',sprintf('%.14g',v(3))),'BPI',sprintf('%.14g',v(1))),'BY',sprintf('%.14g',v(2))),'ALI',sprintf('%.14g',v(4))),'BR',sprintf('%.14g',v(5)));
v0 = [pm.bpi pm.by pm.epsr pm.ali pm.br];
for sd = SDS
    fid = fopen('rnsd_prof.mod','w'); fwrite(fid, fill(tp, sd, v0)); fclose(fid);
    evalc('dynare rnsd_prof.mod noclearall');
    global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info
    if isempty(options_.qz_criterium), options_.qz_criterium = 1 + 1e-6; end
    bounds = prior_bounds(bayestopt_, options_.prior_trunc);
    nm = cellstr(bayestopt_.name); ix = @(n) find(strcmp(nm, n));
    x0 = get_all_parameters(estim_params_, M_); x0 = x0(:);
    lp = @(x) -dsge_likelihood(x, dataset_, dataset_info, options_, M_, estim_params_, bayestopt_, bounds, oo_.dr, oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);
    free = [ix('bpi') ix('by') ix('stderr epsr') ix('ali') ix('br')];
    obj = @(z) penal(-lp(setfree(x0, free, z)));
    lp0 = lp(x0);
    z = x0(free); opts = optimset('MaxFunEvals', 500, 'MaxIter', 500, 'TolX', 1e-5, 'TolFun', 1e-3, 'Display', 'off');
    for pass = 1:3, [z, fv] = fminsearch(obj, z, opts); end
    fprintf('\nSD=%.4f: logpost at posterior-mean params %.2f -> re-optimised %.2f | bpi=%.3f by=%.4f sd_epsr=%.4f ali=%.4f br=%.3f\n', sd, lp0, -fv, z(1), z(2), z(3), z(4), z(5));
    ok = false;
    for shrink = [1 0.9 0.75 0.5]
        zz = x0(free) + shrink*(z - x0(free));
        fid = fopen('rnsd_run.mod','w'); fwrite(fid, fill(tr, sd, zz(:)')); fclose(fid);
        try
            s = evalc('dynare rnsd_run.mod noclearall'); ok = true; break
        catch
        end
    end
    if ~ok, fprintf('   Dynare run failed at all shrink levels\n'); continue, end
    k = regexp(s,'log posterior \(or likelihood\): ([-\d\.]+)','tokens');
    sm = oo_.occbin.smoother; h = sm.realtime_regime_history; b = false(T,1); for t=1:T, b(t) = h(t).regime(1)==1; end
    d = diff([0;b;0]); st = find(d==1); en = find(d==-1)-1; sp = '';
    for i=1:numel(st), sp = [sp sprintf(' %s-%s', ql{st(i)}, ql{en(i)})]; end
    fprintf('   regime path (params shrunk to %.2f of the way to the optimum; Dynare logpost %s): %s\n', shrink, k{1}{1}, sp);
    fprintf('   flagged share 2017-19: %.2f ; data-implied ZLB quarters flagged: %d of 41\n', mean(b(yr'>=2017&yr'<=2019)), sum(b(1:32))+sum(b(46:54)));
end

function v = penal(v), if ~isfinite(v) || v > 1e5, v = 1e5; end, end
function x = setfree(x, idx, z), x(idx) = z(:); end
