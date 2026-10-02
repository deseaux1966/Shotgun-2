function resampled_output = resample(particles,weights,ParticleOptions,u_external)
% resampled_output = resample(particles,weights,ParticleOptions,u_external)
% Resamples particles.
% if particles = 0, returns the resampling index (except for smooth resampling)
% Otherwise, returns the resampled particles set.
%
% INPUTS
%  - particles              [double]    n*1 vector of particles
%  - weights                [double]    n*1 vector of particles' weights.
%  - ParticleOptions        [structure] filter options
%  - u_external              [double]   (optional) externally supplied uniform
%                                       random draw(s) to use in place of an
%                                       internal rand() call -- see note below.
%
% OUTPUTS
%  - resampled_output       [double]    vector of resample particles
%
% This function is called by: sequential_importance_particle_filter
% This function calls: residual_resampling, traditional_resampling

% Copyright © 2011-2025 Dynare Team
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
% PROJECT-LOCAL PATCH (Shotgun-2): adds an optional u_external argument so
% the correlated-pseudo-marginal MCMC experiment
% (sequential_importance_particle_filter_cpm.m, in this same directory)
% can supply the resampling draw itself (as normcdf(z) of a standard
% normal z that is part of the same correlated random vector driving
% particle propagation), rather than consuming an independent internal
% rand() call every period. With no 4th argument, behaves exactly as
% stock Dynare (internal rand()/rand(size(weights))) -- backward compatible.
% ---------------------------------------------------------------------

defaultmethod = 1; % For residual based method set this variable equal to 0.

if nargin<4
    u_external = [];
end

if defaultmethod
    if ParticleOptions.resampling.method.kitagawa
        if isempty(u_external)
            resampled_output = traditional_resampling(particles,weights,rand);
        else
            resampled_output = traditional_resampling(particles,weights,u_external);
        end
    elseif ParticleOptions.resampling.method.stratified
        if isempty(u_external)
            resampled_output = traditional_resampling(particles,weights,rand(size(weights)));
        else
            resampled_output = traditional_resampling(particles,weights,u_external);
        end
    elseif ParticleOptions.resampling.method.smooth
        if particles==0
            error('Particle = 0 is incompatible with the resampling_method=smooth method!')
        end
        resampled_output = multivariate_smooth_resampling(particles,weights);
    else
        error('Unknown sampling method!')
    end
else
    if ParticleOptions.resampling.method.kitagawa
        if isempty(u_external)
            resampled_output = residual_resampling(particles,weights,rand);
        else
            resampled_output = residual_resampling(particles,weights,u_external);
        end
    elseif ParticleOptions.resampling.method.stratified
        if isempty(u_external)
            resampled_output = residual_resampling(particles,weights,rand(size(weights)));
        else
            resampled_output = residual_resampling(particles,weights,u_external);
        end
    else
        error('Unknown sampling method!')
    end
end
