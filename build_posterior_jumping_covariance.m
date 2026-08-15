% build_posterior_jumping_covariance.m
% Builds a custom MCMC jumping covariance file from this round's own
% pooled posterior draws (taubar=0.14, chi=0.5, hessian-covariance run),
% rather than continuing to reuse the Stage-1 mode's local Hessian.
%
% IMPORTANT (verified against set_mcmc_jumping_covariance.m in the
% Dynare 7.1 source): despite the required variable name
% "jumping_covariance", Dynare treats the loaded matrix as a
% Hessian-like precision object and INVERTS it to get the actual
% proposal covariance used by the sampler:
%   hess = jumping_covariance;
%   invhess = inv(hess./(hsd*hsd'))./(hsd*hsd');   % this is what's used
% So to make the sampler actually propose from the empirical posterior
% covariance Sigma_post, we must save inv(Sigma_post) (the empirical
% precision matrix) under the name "jumping_covariance", not Sigma_post
% itself.

clear; close all;
shotgun_dir = 'C:/Users/gawater/Documents/GitHub/Shotgun 2';
cd(shotgun_dir);

metdir = fullfile(shotgun_dir, 'mmt8bOccBin3_covid_bT_baa', ...
    'mmt8bOccBin3_covid_bT_baa', 'metropolis');
S1 = load(fullfile(metdir, 'mmt8bOccBin3_covid_bT_baa_mh1_blck1.mat'), 'x2');
S2 = load(fullfile(metdir, 'mmt8bOccBin3_covid_bT_baa_mh1_blck2.mat'), 'x2');

burn = 0.5;
n1 = size(S1.x2,1); n2 = size(S2.x2,1);
draws1 = S1.x2(round(burn*n1)+1:end, :);
draws2 = S2.x2(round(burn*n2)+1:end, :);
pooled = [draws1; draws2];
fprintf('Pooled post-burn-in draws: %d x %d parameters\n', size(pooled,1), size(pooled,2));

Sigma_post = cov(pooled);
fprintf('Empirical posterior covariance: %d x %d, rank=%d\n', ...
    size(Sigma_post,1), size(Sigma_post,2), rank(Sigma_post));

% Check positive definiteness of Sigma_post itself (sanity check).
[~, p] = chol(Sigma_post);
if p ~= 0
    error('Empirical posterior covariance is not positive definite (chol failed at row %d).', p);
end
fprintf('Sigma_post is positive definite.\n');

jumping_covariance = inv(Sigma_post);   %#ok<MINV>  % precision matrix; Dynare re-inverts this internally

% Verify Dynare's own positive-definiteness check will pass on what we save.
[~, p2] = chol(jumping_covariance);
if p2 ~= 0
    error('inv(Sigma_post) is not positive definite (chol failed at row %d) -- Dynare would reject this file.', p2);
end
fprintf('inv(Sigma_post) (the saved "jumping_covariance") is positive definite -- will pass Dynare''s check.\n');

outfile = fullfile(shotgun_dir, 'mmt8bOccBin3_covid_bT_baa', 'posterior_jumping_covariance.mat');
save(outfile, 'jumping_covariance');
fprintf('Saved: %s\n', outfile);
