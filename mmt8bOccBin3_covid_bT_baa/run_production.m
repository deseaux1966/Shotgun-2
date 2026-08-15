% run_production.m
% Full production Bayesian estimation: true OccBin ZLB likelihood,
% mode_compute=0 from the validated Stage 1 mode, mh_jscale=0.03
% (tuned in Stage 2 to clear the 0.05 acceptance-ratio floor on both
% chains), mh_replic=50000 x 2 chains.
addpath('C:\dynare\7.1\matlab');
cd('C:\Users\gawater\Documents\GitHub\Shotgun 2\mmt8bOccBin3_covid_bT_baa');
diary('production_epsL_occbin.log');
fprintf('=== Production: OccBin ZLB, epsL COVID shock, 50000x2 draws, jscale=0.03 ===\n');
occbin_val = true;
try
    dynare mmt8bOccBin3_covid_bT_baa.mod nograph;
    fprintf('=== Production complete ===\n');
catch ME
    fprintf('=== Production FAILED: %s ===\n', ME.message);
end
diary off;
exit;
