## Reference implementations

module RefImp

using LinearAlgebra

function vest(Y)
    return map(Y) do Yl
        d, nl = size(Yl)
        return norm(Yl)^2 / (d * nl)
    end
end

function λest_inv(Y, v)
    n, L = size.(Y, 2), length(Y)
    w = [(1 / v[l]) / sum(n[lp] / v[lp] for lp in 1:L) for l in 1:L]

    Yw = reduce(hcat, sqrt(w[l]) * Y[l] for l in 1:L)
    return svdvals(Yw) .^ 2
end

Ξ(λ, v; c) = -(v + v / c - λ) / 2 + sqrt((v + v / c - λ)^2 - 4 * v^2 / c) / 2

function λest(Y; v=vest(Y))
    d, n, L = (only ∘ unique)(size.(Y, 1)), size.(Y, 2), length(Y)
    p = n ./ sum(n)

    λinv = λest_inv(Y, v)
    c = sum(n) / d
    vb = inv(sum(p[l] / v[l] for l in 1:L))

    k = count(>(vb * (1 + 1 / sqrt(c))^2), λinv)
    return [Ξ(λinv[i], vb; c) for i in 1:k]
end

end
