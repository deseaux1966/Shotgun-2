function [LIK,lik] = sequential_importance_particle_filter(ReducedForm,Y,start,ParticleOptions,ThreadsOptions, options_, M_)
% [LIK,lik] = sequential_importance_particle_filter(ReducedForm,Y,start,ParticleOptions,ThreadsOptions, options_, M_)
% Evaluates the likelihood of a nonlinear model with a particle filter employing a sequential importance sampling approach
%
% INPUTS
%  - ReducedForm            [structure] decision rules
%  - Y                      [double]    dataset
%  - start                  [integer]   first observation for likelihood evaluation
%  - ParticleOptions        [structure] filter options
%  - ThreadsOptions         [structure] options for threading of mex files
%  - options_               [structure] describing the options
%  - M_                     [structure] describing the model
%
% OUTPUTS
% - LIK                [double]    scalar, likelihood
% - lik                [double]    (T-s+1)×1 vector, density of observations in each period.
%
% References:
% Implementation is e.g. described in Godsill/Doucet/West (2004): "Monte Carlo Smoothing for Nonlinear Time Series",
% Journal of the American Statistical Association, March 2004, 99(465)

% Copyright © 2011-2026 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <https://www.gnu.org/licenses/>.
%
% ---------------------------------------------------------------------
% PROJECT-LOCAL PATCH (Shotgun-2, mmt8bOccBin3_covid_bT_baa_softzlb*):
% stock Dynare 7.2 calls
%     set_dynare_seed_local_options([],false,'default')
% here, which unconditionally resets MATLAB's global RNG to a FIXED seed
% (seed=0, per set_dynare_seed_local_options.m's 'default' branch) on
% EVERY likelihood evaluation. Confirmed directly: three different
% set_dynare_seed() calls made before invoking the (unpatched) likelihood
% produced bit-identical results to 12 decimal places. This is a
% deliberate common-random-numbers design that helps mode-finding
% (see mmt8bOccBin3_covid_bT_baa_softzlb_mode.mod) but is invalid for
% Metropolis-Hastings estimation: pseudo-marginal MCMC (Andrieu & Roberts,
% 2009, Annals of Statistics 37(2), 697-725) requires the likelihood
% estimator to be unbiased AND independently redrawn at every iteration.
%
% Fix below: seed a fresh mt19937ar stream from a persistent,
% monotonically-incrementing call counter instead of a fixed constant, so
% every call gets a distinct (and, across the many calls in an MH chain,
% effectively independent) random stream. This is deterministic given the
% counter's value -- useful for debugging -- but never repeats within a
% session, unlike stock Dynare's fixed seed=0. The counter is intentionally
% NOT reset by this function; call `clear sequential_importance_particle_filter`
% (or `clear functions`) to reset it, e.g. between independent diagnostic runs.
%
% This file lives in Shotgun-2/dynare_patches, added to the MATLAB path
% ahead of the installed C:\dynare\7.2\matlab\nonlinear-filters copy so it
% shadows the original; the Dynare installation itself is untouched.
% ---------------------------------------------------------------------

% Set default value for start
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

% Initialization of the likelihood.
const_lik = log(2*pi)*number_of_observed_variables +log(det(ReducedForm.H)) ;
lik  = NaN(sample_size,1);

% filter out singular measurement error case
if rcond(ReducedForm.H) < 1e-12
    LIK = NaN;
    return
end

StateVectorVarianceSquareRoot = chol(ReducedForm.StateVectorVariance)';%reduced_rank_cholesky(ReducedForm.StateVectorVariance)';
state_variance_rank = size(StateVectorVarianceSquareRoot,2); % Get the rank of StateVectorVarianceSquareRoot

% Factorize the covariance matrix of the structural innovations
Q_lower_triangular_cholesky = chol(ReducedForm.Q)';

% PATCHED: fresh seed per call (see header note above), not a fixed default.
persistent sipf_patch_call_counter
if isempty(sipf_patch_call_counter)
    sipf_patch_call_counter = 0;
end
sipf_patch_call_counter = sipf_patch_call_counter + 1;
set_dynare_seed_local_options([],false,'mt19937ar',sipf_patch_call_counter);

% Initialization of the weights across particles.
weights = ones(1,number_of_particles)/number_of_particles ;
StateVectors = bsxfun(@plus,StateVectorVarianceSquareRoot*randn(state_variance_rank,number_of_particles),ReducedForm.StateVectorMean);
if ParticleOptions.pruning && ~(options_.order ==1)
    if options_.order == 2
        StateVectors_ = StateVectors;
        mf0_ = ReducedForm.mf0;
    elseif options_.order == 3
        StateVectors_ = repmat(StateVectors,3,1);
        mf0_ = repmat(ReducedForm.mf0,1,3);
        mask2 = number_of_state_variables+1:2*number_of_state_variables;
        mask3 = 2*number_of_state_variables+1:3*number_of_state_variables;
        mf0_(mask2) = mf0_(mask2)+size(ReducedForm.ghx,1);
        mf0_(mask3) = mf0_(mask3)+2*size(ReducedForm.ghx,1);
    else
        error('Pruning is not available for orders > 3');
    end
else
    StateVectors_ = StateVectors;
    mf0_ = ReducedForm.mf0;
end

% Loop over observations
for t=1:sample_size
    epsilon = Q_lower_triangular_cholesky*randn(number_of_structural_innovations,number_of_particles);
    if ParticleOptions.pruning
        [tmp, tmp_]=iterate_law_of_motion(StateVectors,epsilon,ReducedForm,M_,options_,ReducedForm.use_k_order_solver,ParticleOptions.pruning,StateVectors_);
    else
        tmp=iterate_law_of_motion(StateVectors,epsilon,ReducedForm,M_,options_,ReducedForm.use_k_order_solver,ParticleOptions.pruning);
    end

    PredictionError = bsxfun(@minus,Y(:,t),tmp(ReducedForm.mf1,:));

    lnw = -.5*(const_lik+sum(PredictionError.*(ReducedForm.H\PredictionError),1));

    dfac = max(lnw);
    wtilde = weights.*exp(lnw-dfac);
    lik(t) = log(sum(wtilde))+dfac;
    weights = wtilde/sum(wtilde);
    if (ParticleOptions.resampling.status.generic && neff(weights)<ParticleOptions.resampling.threshold*sample_size) || ParticleOptions.resampling.status.systematic
        if ParticleOptions.pruning
            temp = resample([tmp(ReducedForm.mf0,:)' tmp_(mf0_,:)'],weights',ParticleOptions);
            StateVectors = temp(:,1:number_of_state_variables)';
            StateVectors_ = temp(:,number_of_state_variables+1:end)';
        else
            StateVectors = resample(tmp(ReducedForm.mf0,:)',weights',ParticleOptions)';
        end
        weights = ones(1,number_of_particles)/number_of_particles;
    elseif ParticleOptions.resampling.status.none
        StateVectors = tmp(ReducedForm.mf0,:);
        if ParticleOptions.pruning
            StateVectors_ = tmp_(mf0_,:);
        end
    end
end

LIK = -sum(lik(start:end));
