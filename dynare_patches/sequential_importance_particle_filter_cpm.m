function [LIK,lik] = sequential_importance_particle_filter_cpm(ReducedForm,Y,start,ParticleOptions,ThreadsOptions, options_, M_, U)
% [LIK,lik] = sequential_importance_particle_filter_cpm(ReducedForm,Y,start,ParticleOptions,ThreadsOptions, options_, M_, U)
%
% CORRELATED-PSEUDO-MARGINAL TEST VARIANT of
% sequential_importance_particle_filter.m. Every source of randomness the
% stock filter draws internally (initial particle cloud, per-period
% structural innovations, and the single kitagawa resampling uniform used
% every period under this project's default options -- systematic
% resampling status, kitagawa method) is instead taken from the externally
% supplied struct U:
%   - U.init     [state_variance_rank x N]          standard normal
%   - U.eps      [n_struct_innov x N x T]            standard normal
%   - U.resample [1 x T]                             standard normal
%                (converted to a uniform via normcdf before use, so it can
%                 be updated by an ordinary Crank-Nicolson step in normal
%                 space along with everything else)
%
% Given (theta, U), this function is a PURE, deterministic map -- no
% internal randn()/rand() calls at all. This is what correlated
% pseudo-marginal MCMC (Deligiannidis, Doucet & Pitt, 2018, "The
% Correlated Pseudo-Marginal Method," JRSS-B 80(5), 839-870) needs: an
% outer sampler proposes theta' via the usual random walk AND proposes
% U' = rho*U + sqrt(1-rho^2)*Z (Z fresh iid N(0,1), same shape as U) via a
% Crank-Nicolson update, so the two log-likelihood estimates entering the
% MH ratio are evaluated under highly correlated randomness -- the
% difference has much lower variance than two fully independent estimates,
% even though each individual estimate's own variance is unchanged.
%
% This is a TEST/DIAGNOSTIC harness only (see
% run_cpm_variance_test.m in this directory): it checks whether
% correlating U actually buys the variance reduction the method promises
% on THIS model, at feasible particle counts, before committing to
% building a full custom MH driver around it. Pruning (options_.particle
% order>1 state augmentation) is not supported here -- the live project
% options have ParticleOptions.pruning=0, so this is not a live limitation
% for this test, only for any future extension of it.
%
% Copyright structure and non-random-generation logic follow
% sequential_importance_particle_filter.m (Copyright (C) 2011-2026 Dynare
% Team, GPLv3+); only the source of randomness is changed.

if isempty(start)
    start = 1;
end

sample_size = size(Y,2);
number_of_state_variables = length(ReducedForm.mf0);
number_of_observed_variables = length(ReducedForm.mf1);
number_of_structural_innovations = length(ReducedForm.Q);
number_of_particles = ParticleOptions.number_of_particles;

if isempty(ReducedForm.H)
    ReducedForm.H = 0;
end

const_lik = log(2*pi)*number_of_observed_variables +log(det(ReducedForm.H)) ;
lik  = NaN(sample_size,1);

if rcond(ReducedForm.H) < 1e-12
    LIK = NaN;
    return
end

StateVectorVarianceSquareRoot = chol(ReducedForm.StateVectorVariance)';
state_variance_rank = size(StateVectorVarianceSquareRoot,2);

Q_lower_triangular_cholesky = chol(ReducedForm.Q)';

if ParticleOptions.pruning
    error('sequential_importance_particle_filter_cpm: pruning not supported in this test harness.')
end

if size(U.init,1)~=state_variance_rank || size(U.init,2)~=number_of_particles
    error('sequential_importance_particle_filter_cpm: U.init has the wrong size.')
end
if size(U.eps,1)~=number_of_structural_innovations || size(U.eps,2)~=number_of_particles || size(U.eps,3)~=sample_size
    error('sequential_importance_particle_filter_cpm: U.eps has the wrong size.')
end
if numel(U.resample)~=sample_size
    error('sequential_importance_particle_filter_cpm: U.resample has the wrong size.')
end

weights = ones(1,number_of_particles)/number_of_particles ;
StateVectors = bsxfun(@plus,StateVectorVarianceSquareRoot*U.init,ReducedForm.StateVectorMean);
StateVectors_ = StateVectors;
mf0_ = ReducedForm.mf0;

for t=1:sample_size
    epsilon = Q_lower_triangular_cholesky*U.eps(:,:,t);
    tmp=iterate_law_of_motion(StateVectors,epsilon,ReducedForm,M_,options_,ReducedForm.use_k_order_solver,ParticleOptions.pruning);

    PredictionError = bsxfun(@minus,Y(:,t),tmp(ReducedForm.mf1,:));

    lnw = -.5*(const_lik+sum(PredictionError.*(ReducedForm.H\PredictionError),1));

    dfac = max(lnw);
    wtilde = weights.*exp(lnw-dfac);
    lik(t) = log(sum(wtilde))+dfac;
    weights = wtilde/sum(wtilde);
    if (ParticleOptions.resampling.status.generic && neff(weights)<ParticleOptions.resampling.threshold*sample_size) || ParticleOptions.resampling.status.systematic
        u_resample_t = normcdf(U.resample(t));
        StateVectors = resample(tmp(ReducedForm.mf0,:)',weights',ParticleOptions,u_resample_t)';
        weights = ones(1,number_of_particles)/number_of_particles;
    elseif ParticleOptions.resampling.status.none
        StateVectors = tmp(ReducedForm.mf0,:);
    end
end

LIK = -sum(lik(start:end));
