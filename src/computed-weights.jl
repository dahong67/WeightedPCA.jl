## Weighted PCA with computed weights

# Abstract type for computed weights
"""
    ComputedWeights

Abstract supertype for weights that are computed from properties of the data.
"""
abstract type ComputedWeights end

# Uniform weights
"""
    UniformWeights <: ComputedWeights

Uniform weighting, i.e., `w[l] = 1`.
Corresponds to conventional (unweighted) PCA.
"""
struct UniformWeights <: ComputedWeights end

_wpca(Y, i, ::UniformWeights) =
    _wpca(Y, i, [one(eltype(eltype(Y))) for _ in 1:length(Y)])

# Inverse noise variance weights
"""
    InverseVarianceWeights <: ComputedWeights

Inverse noise variance weighting, i.e., `w[l] = 1/v[l]`.

# Constructors
+ `InverseVarianceWeights(v=noisevar)` for known noise variances `noisevar`
+ `InverseVarianceWeights()` for unknown noise variances; noise variances will be estimated from data
"""
struct InverseVarianceWeights{Tv} <: ComputedWeights
    v::Tv    # vector of noise variances
end
InverseVarianceWeights(; v=NoiseNormEstimator()) =
    InverseVarianceWeights(v)

_wpca(Y, i, method::InverseVarianceWeights{<:AbstractVector{<:Real}}) =
    _wpca(Y, i, inv.(method.v))
_wpca(Y, i, method::InverseVarianceWeights{<:AbstractNoiseVarEstimator}) =
    _wpca(Y, i, InverseVarianceWeights(estimatev(Y, method.v)))

# Optimal weights
"""
    OptimalWeights <: ComputedWeights

Optimal weighting, i.e., `w[l] = 1/v[l] * 1/(1+v[l]/λ)`.

# Constructors
+ `OptimalWeights(v=noisevar, λ=signalvar)` for known noise variances `noisevar` and signal variance `signalvar`
+ `OptimalWeights(λ=signalvar)` for known signal variance `signalvar`; noise variances will be estimated from data
+ `OptimalWeights(v=noisevar)` for known noise variances `noisevar`; signal variance will be estimated from data
+ `OptimalWeights()` for unknown noise and signal variances; noise and signal variances will be estimated from data
"""
struct OptimalWeights{Tv,Tλ} <: ComputedWeights
    v::Tv    # vector of noise variances
    λ::Tλ    # signal variance
end
OptimalWeights(; v=NoiseNormEstimator(), λ=InvNoiseWeightedShrinkageEstimator()) =
    OptimalWeights(v, λ)

_wpca(Y, i::Number, method::OptimalWeights{<:AbstractVector{<:Real},<:Real}) =
    _wpca(Y, i, inv.(method.v) .* inv.(one(method.λ) .+ method.v ./ method.λ))
_wpca(Y, i::Number, method::OptimalWeights{<:AbstractNoiseVarEstimator,<:Real}) =
    _wpca(Y, i, OptimalWeights(estimatev(Y, method.v), method.λ))
_wpca(Y, i::Number, method::OptimalWeights{<:AbstractVector{<:Real},<:AbstractSignalVarEstimator}) =
    _wpca(Y, i, OptimalWeights(method.v, estimateλ(Y, i, method.v, method.λ)))
_wpca(Y, i, method::OptimalWeights{<:AbstractNoiseVarEstimator,<:AbstractSignalVarEstimator}) =
    _wpca(Y, i, OptimalWeights(estimatev(Y, method.v), method.λ))
