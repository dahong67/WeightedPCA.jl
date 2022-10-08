## Variance estimators

# Noise variance estimators
"""
    AbstractNoiseVarEstimator

Abstract supertype for noise variance estimators.
Should implement `estimatev`.
"""
abstract type AbstractNoiseVarEstimator end

"""
    estimatev(Y,method::AbstractNoiseVarEstimator) -> Vector{<:Real}

Estimate the vector of noise variances from data `Y`
using the estimator `method`.
"""
function estimatev end

"""
    NoiseNormEstimator <: AbstractNoiseVarEstimator

Norm-based estimate of noise variance: `vh[l] = norm(Y[l])^2/length(Y[l])`
"""
struct NoiseNormEstimator <: AbstractNoiseVarEstimator end
estimatev(Y, ::NoiseNormEstimator) = [norm(Yl)^2 / length(Yl) for Yl in Y]

# Signal variance estimators
"""
    AbstractSignalVarEstimator

Abstract supertype for signal variance estimators.
Should implement `estimateλ`.
"""
abstract type AbstractSignalVarEstimator end

"""
    estimateλ(Y,v,method::AbstractSignalVarEstimator)

Estimate the vector of signal variances from data `Y`
using the estimator `method`.
"""
function estimateλ end

"""
    estimateλ(Y,i,v,method::AbstractSignalVarEstimator)

Estimate the `i`th signal variance from data `Y`
using the estimator `method`.
"""
function estimateλ(Y, i::Integer, v, method::AbstractSignalVarEstimator)
    λ = estimateλ(Y, v, method)
    i <= length(λ) || error("could not estimate λi, component may be too weak")
    return λ[i]
end

"""
    InvNoiseWeightedShrinkageEstimator <: AbstractSignalVarEstimator

Inverse noise variance weighted shrinkage-based estimate of signal variances:
`λh[i] = Ξ(λhinv[i])`
+ `Ξ(λ) = -(vb+vb/c-λ)/2 + sqrt((vb+vb/c-λ)^2-4*vb^2/c)/2` is the shrinkage
+ `vb = ( Σ_l p[l]/v[l] )^(-1)` where `p[l] = n[l]/n`
+ `c = n/d` is the data aspect ratio
+ `λhinv = eigvals( Σ_l (1/v[l])/(n[1]/v[1]+⋯+n[L]/v[L]) Y[l]*Y[l]' )`
  are inverse noise variance weighted eigenvalues
"""
struct InvNoiseWeightedShrinkageEstimator <: AbstractSignalVarEstimator end
function estimateλ(Y, v, ::InvNoiseWeightedShrinkageEstimator)
    # Shrinkage function
    Ξ = (λ, v, c) -> -(v + v / c - λ) / 2 + sqrt((v + v / c - λ)^2 - 4 * v^2 / c) / 2

    # Compute inverse noise variance weighted eigenvalues
    n, L = size.(Y, 2), length(Y)
    w = [(1 / v[l]) / sum(n[lp] / v[lp] for lp in 1:L) for l in 1:L]
    Yw = reduce(hcat, sqrt(w[l]) * Y[l] for l in 1:L)
    λhinv = svdvals(Yw) .^ 2

    # Compute number of components above (inv weighted) phase transition
    d = (only ∘ unique)(size.(Y, 1))
    p, c = n ./ sum(n), sum(n) / d
    vb = inv(sum(p[l] / v[l] for l in 1:L))
    k = count(>(vb * (1 + 1 / sqrt(c))^2), λhinv)

    return [Ξ(λhinv[i], vb, c) for i in 1:k]
end
