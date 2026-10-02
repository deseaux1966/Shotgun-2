% run_rnfixed_diag.m -- short diagnostic (800 draws x 2 chains) for
% mmt8bOccBin3_covid_bT_baa_rnfixed_diag.mod (eps_rn_me fixed at 0.0015).
restoredefaultpath; addpath('C:\dynare\7.2\matlab');
cd('C:\Users\gawater\Documents\GitHub\Shotgun-2');
tic;
dynare mmt8bOccBin3_covid_bT_baa_rnfixed_diag.mod
fprintf('\nDIAGNOSTIC RUN FINISHED in %.1f minutes\n', toc/60);
