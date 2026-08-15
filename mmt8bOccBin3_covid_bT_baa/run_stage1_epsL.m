% run_stage1_epsL.m
% Stage 1: linear-KF posterior mode for the epsL (simple labor-supply
% shock) version of mmt8bOccBin3_covid_bT_baa.mod. OccBin is off here so
% csminwel has a smooth likelihood; this mode is then reused (mode_file)
% as the Stage 2 starting point once OccBin is switched on.
%
% mode_compute/mh_replic/mh_jscale are hardcoded literals inside the .mod
% file's estimation(...) command (Dynare requires literal tokens there);
% edit them directly in the .mod for different stages.
addpath('C:\dynare\7.1\matlab');
cd('C:\Users\gawater\Documents\GitHub\Shotgun 2\mmt8bOccBin3_covid_bT_baa');
diary('stage1_epsL.log');
fprintf('=== Stage 1: linear KF mode-finding + short MCMC diagnostic (epsL model) ===\n');
occbin_val = false;
skip_save  = true;
try
    dynare mmt8bOccBin3_covid_bT_baa.mod nograph;
    fprintf('=== Stage 1 complete ===\n');
catch ME
    fprintf('=== Stage 1 FAILED: %s ===\n', ME.message);
end
diary off;
exit;
