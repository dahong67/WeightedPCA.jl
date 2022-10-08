using WeightedPCA
using Test
using WeightedPCA: estimatev, NoiseNormEstimator
using WeightedPCA: estimateλ, InvNoiseWeightedShrinkageEstimator
using LinearAlgebra
using StableRNGs

include("ref-imp.jl") # reference implementations

@testset "Noise variance estimators" begin
    rng = StableRNG(0)

    # Data setup
    c, v = [4, 8], [1, 3]
    λ1 = 1

    d = 20
    n = c .* d
    k = 1
    L = length.((c, v)) |> only ∘ unique

    # Generate data
    u1 = normalize(randn(rng, d))
    F = reshape(sqrt(λ1) * u1, :, 1)
    Y = [F * randn(rng, k, n[l]) + sqrt(v[l]) * randn(rng, d, n[l]) for l in 1:L]

    @testset "NoiseNormEstimator" begin
        @test estimatev(Y, NoiseNormEstimator()) == RefImp.vest(Y)
    end
end

@testset "Signal variance estimators" begin
    rng = StableRNG(0)

    # Data setup
    c, v = [4, 8], [1, 3]
    λ1 = 1

    d = 20
    n = c .* d
    k = 1
    L = length.((c, v)) |> only ∘ unique

    # Generate data
    u1 = normalize(randn(rng, d))
    F = reshape(sqrt(λ1) * u1, :, 1)
    Y = [F * randn(rng, k, n[l]) + sqrt(v[l]) * randn(rng, d, n[l]) for l in 1:L]

    @testset "InvNoiseWeightedShrinkageEstimator" begin
        # true v
        λr = RefImp.λest(Y; v=v)
        @test estimateλ(Y, v, InvNoiseWeightedShrinkageEstimator()) == λr
        for i in axes(λr, 1)
            @test estimateλ(Y, i, v, InvNoiseWeightedShrinkageEstimator()) == λr[i]
        end
        @test_throws ErrorException estimateλ(Y, k + 1, v, InvNoiseWeightedShrinkageEstimator())

        # estimated v
        λr = RefImp.λest(Y)
        @test estimateλ(Y, estimatev(Y, NoiseNormEstimator()), InvNoiseWeightedShrinkageEstimator()) == λr
        for i in axes(λr, 1)
            @test estimateλ(Y, i, estimatev(Y, NoiseNormEstimator()), InvNoiseWeightedShrinkageEstimator()) == λr[i]
        end
        @test_throws ErrorException estimateλ(Y, k + 1, estimatev(Y, NoiseNormEstimator()), InvNoiseWeightedShrinkageEstimator())
    end
end

@testset "Weighted PCA: manually set weights" begin
    rng = StableRNG(0)

    # Data setup
    c, v = [4, 8], [1, 3]
    λ1 = 1

    d = 20
    n = c .* d
    k = 1
    L = length.((c, v)) |> only ∘ unique

    # Generate data
    u1 = normalize(randn(rng, d))
    F = reshape(sqrt(λ1) * u1, :, 1)
    Y = [F * randn(rng, k, n[l]) + sqrt(v[l]) * randn(rng, d, n[l]) for l in 1:L]

    @testset "i=$i, w=$w" for i in 1:3, w in [[1, 2], [2, 3]]
        @test wpca(Y, i, w) == let w = w
            Yw = reduce(hcat, sqrt(w[l]) * Y[l] for l in 1:L)
            svd(Yw).U[:, i]
        end
    end
end

@testset "Weighted PCA: uniform weights" begin
    rng = StableRNG(0)

    # Data setup
    c, v = [4, 8], [1, 3]
    λ1 = 1

    d = 20
    n = c .* d
    k = 1
    L = length.((c, v)) |> only ∘ unique

    # Generate data
    u1 = normalize(randn(rng, d))
    F = reshape(sqrt(λ1) * u1, :, 1)
    Y = [F * randn(rng, k, n[l]) + sqrt(v[l]) * randn(rng, d, n[l]) for l in 1:L]

    @testset "i=$i" for i in 1:3
        @test wpca(Y, i, UniformWeights()) == svd(reduce(hcat, Y)).U[:, i]
    end
end

@testset "Weighted PCA: inverse noise variance weights" begin
    rng = StableRNG(0)

    # Data setup
    c, v = [4, 8], [1, 3]
    λ1 = 1

    d = 20
    n = c .* d
    k = 1
    L = length.((c, v)) |> only ∘ unique

    # Generate data
    u1 = normalize(randn(rng, d))
    F = reshape(sqrt(λ1) * u1, :, 1)
    Y = [F * randn(rng, k, n[l]) + sqrt(v[l]) * randn(rng, d, n[l]) for l in 1:L]

    @testset "i=$i" for i in 1:3
        # true v
        @test wpca(Y, i, InverseVarianceWeights(v=v)) == let v = v
            w = inv.(v)
            Yw = reduce(hcat, sqrt(w[l]) * Y[l] for l in 1:L)
            svd(Yw).U[:, i]
        end

        # estimated v
        @test wpca(Y, i, InverseVarianceWeights()) == let v = RefImp.vest(Y)
            w = inv.(v)
            Yw = reduce(hcat, sqrt(w[l]) * Y[l] for l in 1:L)
            svd(Yw).U[:, i]
        end
    end
end

@testset "Weighted PCA: optimal weights" begin
    rng = StableRNG(0)

    # Data setup
    c, v = [4, 8], [1, 3]
    λ1 = 1.5

    d = 20
    n = c .* d
    k = 1
    L = length.((c, v)) |> only ∘ unique

    # Generate data
    u1 = normalize(randn(rng, d))
    F = reshape(sqrt(λ1) * u1, :, 1)
    Y = [F * randn(rng, k, n[l]) + sqrt(v[l]) * randn(rng, d, n[l]) for l in 1:L]

    # true v, true λi
    @test wpca(Y, 1, OptimalWeights(v=v, λ=λ1)) ≈ let v = v, λi = λ1
        w = inv.(v .* (λi .+ v))
        Yw = reduce(hcat, sqrt(w[l]) * Y[l] for l in 1:L)
        svd(Yw).U[:, 1]
    end

    # true v, estimated λi
    @test wpca(Y, 1, OptimalWeights(v=v)) ≈ let v = v, λi = RefImp.λest(Y; v=v)[1]
        w = inv.(v .* (λi .+ v))
        Yw = reduce(hcat, sqrt(w[l]) * Y[l] for l in 1:L)
        svd(Yw).U[:, 1]
    end

    # estimated v, true λi
    @test wpca(Y, 1, OptimalWeights(λ=λ1)) ≈ let v = RefImp.vest(Y), λi = λ1
        w = inv.(v .* (λi .+ v))
        Yw = reduce(hcat, sqrt(w[l]) * Y[l] for l in 1:L)
        svd(Yw).U[:, 1]
    end

    # estimated v, estimated λi
    @test wpca(Y, 1, OptimalWeights()) ≈ let v = RefImp.vest(Y), λi = RefImp.λest(Y)[1]
        w = inv.(v .* (λi .+ v))
        Yw = reduce(hcat, sqrt(w[l]) * Y[l] for l in 1:L)
        svd(Yw).U[:, 1]
    end
end
