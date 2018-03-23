function estimate(entity::FundamentalMatrix, matches...)
    ℳ, ℳʹ = matches
    N = length(ℳ)
    if (N != length(ℳʹ))
          throw(ArgumentError("There should be an equal number of points for each view."))
    end
    (ℳ,𝐓) = hartley_normalization(ℳ)
    (ℳʹ,𝐓ʹ) = hartley_normalization!(ℳʹ)
    𝐀::Matrix{Float64} = moments(FundamentalMatrix(), ℳ, ℳʹ)
    (λ::Float64, f::Vector{Float64}) = smallest_eigenpair(𝐀)
    𝐅::Matrix{Float64} = reshape(f,(3,3))
    # Enforce the rank-2 constraint.
    U,S,V = svd(𝐅)
    S[end] = 0.0
    𝐅 = U*diagm(S)*V'
end
