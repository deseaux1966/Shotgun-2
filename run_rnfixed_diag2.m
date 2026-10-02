% run_rnfixed_diag2.m -- short diagnostic (800 draws x 2 chains), mh_jscale=0.01, for
% mmt8bOccBin3_covid_bT_baa_rnfixed_diag2.mod (eps_rn_me fixed at 0.0015).
restoredefaultpath; addpath('C:\dynare\7.2\matlab');
cd('C:\Users\gawater\Documents\GitHub\Shotgun-2');
tic;
dynare mmt8bOccBin3_covid_bT_baa_rnfixed_diag2.mod
fprintf('\nDIAGNOSTIC RUN FINISHED in %.1f minutes\n', toc/60);
