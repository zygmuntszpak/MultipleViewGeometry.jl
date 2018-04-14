function estimate(entity::FundamentalMatrix, method::DirectLinearTransform, matches...)
    ℳ, ℳʹ = matches
    N = length(ℳ)
    if (N != length(ℳʹ))
          throw(ArgumentError("There should be an equal number of points for each view."))
    end
    (ℳ,𝐓) = transform(HomogeneousCoordinates(),CanonicalToHartley(),ℳ)
    (ℳʹ,𝐓ʹ) = transform(HomogeneousCoordinates(),CanonicalToHartley(),ℳʹ)
    𝐀::Matrix{Float64} = moments(FundamentalMatrix(), ℳ, ℳʹ)
    (λ::Float64, f::Vector{Float64}) = smallest_eigenpair(𝐀)
    𝐅::Matrix{Float64} = reshape(f,(3,3))
    enforce_ranktwo!(𝐅)
    # Transform estimate back to the original (unnormalised) coordinate system.
    𝐅 = 𝐓ʹ'*𝐅*𝐓
end

# TODO fix numerical instability
function estimate(entity::FundamentalMatrix, method::Taubin, matches...)
    ℳ, ℳʹ = matches
    N = length(ℳ)
    if (N != length(ℳʹ))
          throw(ArgumentError("There should be an equal number of points for each view."))
    end
    (ℳ,𝐓) = hartley_normalization(ℳ)
    (ℳʹ,𝐓ʹ) = hartley_normalization(ℳʹ)
    𝐀::Matrix{Float64} = moments(FundamentalMatrix(), ℳ, ℳʹ)
    𝐁::Matrix{Float64} = mean_covariance(FundamentalMatrix(), ℳ, ℳʹ)
    (λ::Float64, f::Vector{Float64}) = smallest_eigenpair(𝐀,𝐁)
    𝐅::Matrix{Float64} = reshape(f,(3,3))
    enforce_ranktwo!(𝐅)
    # Transform estimate back to the original (unnormalised) coordinate system.
    𝐅 = 𝐓ʹ'*𝐅*𝐓
end

function estimate(entity::FundamentalMatrix, method::FundamentalNumericalScheme, Λ::Vector{T1},  matches...) where T1 <: Matrix
    ℳ, ℳʹ = matches
    N = length(ℳ)
    if (N != length(ℳʹ))
          throw(ArgumentError("There should be an equal number of points for each view."))
    end
    if (N != length(Λ))
          throw(ArgumentError("There should be a covariance matrix for each point correspondence."))
    end
    (ℳ,𝐓) = transform(HomogeneousCoordinates(),CanonicalToHartley(),ℳ)
    (ℳʹ,𝐓ʹ) = transform(HomogeneousCoordinates(),CanonicalToHartley(),ℳʹ)
    # Map seed to the normalized coordinate system.
    𝛉 = (inv(𝐓') ⊗ inv(𝐓ʹ')) * method.𝛉₀
    # Map covariance matrices to the normalized coordinate system.
    Λʹ = transform(CovarianceMatrices(), CanonicalToHartley(), Λ , tuple(𝐓,𝐓ʹ))
    for i = 1:method.max_iter
        𝐗 = X(AML(),FundamentalMatrix(), 𝛉, Λʹ, ℳ, ℳʹ)
        (λ::Float64, 𝛉⁺::Vector{Float64}) = smallest_eigenpair(𝐗/N)
        𝛉 = reshape(𝛉⁺,length(𝛉⁺),1)
    end
    𝐅::Matrix{Float64} = reshape(𝛉,(3,3))
    enforce_ranktwo!(𝐅)
    # Transform estimate back to the original (unnormalised) coordinate system.
    𝐅 = 𝐓ʹ'*𝐅*𝐓
end


#𝛉 = reshape(𝛉₀,length(𝛉₀),1)
# function z(entity::FundamentalMatrix, 𝐦::Vector{Float64}, 𝐦ʹ::Vector{Float64})
# 𝐮 = 𝐦 ⊗ 𝐦ʹ
# 𝐮[1:end-1]
# end

function mean_covariance(entity::ProjectiveEntity, matches...)
    ℳ, ℳʹ = matches
    N = length(ℳ)
    if (N != length(ℳʹ))
          throw(ArgumentError("There should be an equal number of points for each view."))
    end
    𝐁 = fill(0.0,(9,9))
    𝚲 = eye(4)
    for correspondence in zip(ℳ, ℳʹ)
        m , mʹ = correspondence
        𝐦  = 𝑛(collect(Float64,m.coords))
        𝐦ʹ = 𝑛(collect(Float64,mʹ.coords))
        ∂ₓ𝐮 = ∂ₓu(entity, 𝐦 , 𝐦ʹ)
        𝐁 = 𝐁 + ∂ₓ𝐮 * 𝚲 * ∂ₓ𝐮'
    end
    𝐁/N
end


function enforce_ranktwo!(𝐅::Matrix{Float64})::Matrix{Float64}
    # Enforce the rank-2 constraint.
    U,S,V = svd(𝐅)
    S[end] = 0.0
    𝐅 = U*diagm(S)*V'
end
