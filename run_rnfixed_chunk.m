% run_rnfixed_chunk.m -- one production chunk (5000 draws x 2 chains) of
% mmt8bOccBin3_covid_bT_baa_rnfixed_prod.mod (eps_rn_me fixed at 0.0015, mh_jscale=0.01).
% Expects the workspace variable loadflag (0 = first chunk, starts new chains; 1 = append).
restoredefaultpath; addpath('C:\dynare\7.2\matlab');
cd('C:\Users\gawater\Documents\GitHub\Shotgun-2');
tic;
dynare('mmt8bOccBin3_covid_bT_baa_rnfixed_prod.mod', sprintf('-DLOAD=%d', loadflag));
fprintf('\nCHUNK FINISHED in %.1f minutes\n', toc/60);
