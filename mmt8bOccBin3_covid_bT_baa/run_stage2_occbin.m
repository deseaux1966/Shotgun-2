% run_stage2_occbin.m
% Stage 2: true OccBin ZLB piecewise-linear likelihood, mode_compute=0
% starting from the Stage 1 linear-KF mode (which matches this exact
% parameter vector -- a mismatched mode_file throws "generated using
% another specification of the model").
addpath('C:\dynare\7.1\matlab');
cd('C:\Users\gawater\Documents\GitHub\Shotgun 2\mmt8bOccBin3_covid_bT_baa');
diary('stage2_occbin.log');
fprintf('=== Stage 2: OccBin ZLB likelihood, mode_compute=0 from Stage 1 mode ===\n');
occbin_val = true;
skip_save  = true;
try
    dynare mmt8bOccBin3_covid_bT_baa.mod nograph;
    fprintf('=== Stage 2 complete ===\n');
catch ME
    fprintf('=== Stage 2 FAILED: %s ===\n', ME.message);
end
diary off;
exit;
