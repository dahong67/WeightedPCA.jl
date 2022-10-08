"Weighted PCA module. Provides weighted principal component analysis (PCA) for data with samples of heterogeneous quality (heteroscedastic noise)."
module WeightedPCA

using LinearAlgebra: norm, svd, svdvals

export wpca
export UniformWeights, InverseVarianceWeights, OptimalWeights

# Main function
"""
    wpca(Y,i,weights=UniformWeights())

Compute `i`th principal component of data `Y` via weighted PCA using `weights`,
i.e., output is the `i`th eigenvector of the weighted sample covariance
`Σ_l w[l] Y[l]*Y[l]'`.
Data `Y` is a list of matrices (each column is a sample).

# Choices for `weights`
+ `UniformWeights()` : uniform weights, i.e., `w[l] = 1` [default]
+ `InverseVarianceWeights([v])` : inverse noise variance weights, i.e., `w[l] = 1/v[l]`
+ `OptimalWeights([v,λ])` : optimal weights for signal with variance `λ`, i.e., `w[l] = 1/v[l] * 1/(1+v[l]/λ)`

The `weights` can also be manually set by passing in an `AbstractVector{<:Real}`.

See also: [`UniformWeights`](@ref), [`InverseVarianceWeights`](@ref), [`OptimalWeights`](@ref).
"""
wpca(Y, i::Integer, weights=UniformWeights()) = _wpca(Y, i, weights)

function _wpca(Y, i, w::AbstractVector{<:Real})
    axes(Y) == axes(w) || throw(DimensionMismatch("`axes(Y)` must match `axes(w)`"))
    Yw = reduce(hcat, sqrt(wl) * Yl for (wl, Yl) in zip(w, Y))
    return svd(Yw).U[:, i]
end

# Variance estimators and weighted PCA with computed weights
include("variance-estimators.jl")
include("computed-weights.jl")

end
