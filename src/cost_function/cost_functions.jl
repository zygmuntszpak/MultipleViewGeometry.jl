function cost(c::CostFunction, entity::FundamentalMatrix, 𝛉::Matrix,Λ::Vector{T1}, matches...) where T1 <: Matrix
    ℳ, ℳʹ = matches
    N = length(ℳ)
    if (N != length(ℳʹ))
          throw(ArgumentError("There should be an equal number of points for each view."))
    end
    if (N != length(Λ))
          throw(ArgumentError("There should be a covariance matrix for each point correspondence."))
    end
    Jₐₘₗ = fill(0.0,(1,1))
    for correspondence in zip(ℳ, ℳʹ,Λ)
        m , mʹ, 𝚲 = correspondence
        𝐦  = 𝑛(collect(Float64,m.coords))
        𝐦ʹ = 𝑛(collect(Float64,mʹ.coords))
        𝐮ₓ =  uₓ(entity, 𝐦 , 𝐦ʹ)
        ∂ₓ𝐮 = ∂ₓu(entity, 𝐦 , 𝐦ʹ)
        𝐁 =  ∂ₓ𝐮 * 𝚲 * ∂ₓ𝐮'
        𝚺⁻¹ = inv(𝛉' * 𝐁 * 𝛉)
        Jₐₘₗ = Jₐₘₗ + 𝛉' * 𝐮ₓ * 𝚺⁻¹ * 𝐮ₓ' * 𝛉
    end
    Jₐₘₗ[1]
end

function ∂cost(c::CostFunction, entity::FundamentalMatrix, 𝛉::Matrix,Λ::Vector{T1}, matches...) where T1 <: Matrix
𝐗 = X(c, entity, 𝛉, Λ, matches)
2*𝐗*𝛉
end


function X(c::CostFunction, entity::FundamentalMatrix, 𝛉::Matrix, Λ::Vector{T1},  matches...) where T1 <: Matrix
    ℳ, ℳʹ = matches
    N = length(ℳ)
    if (N != length(ℳʹ))
          throw(ArgumentError("There should be an equal number of points for each view."))
    end
    if (N != length(Λ))
          throw(ArgumentError("There should be a covariance matrix for each point correspondence."))
    end
    _X(c, entity, 𝛉, Λ, ℳ, ℳʹ)
end

function _X(c::CostFunction, entity::ProjectiveEntity, 𝛉::Matrix, Λ::Vector{T1}, 𝒟...) where T1 <: Matrix
    l = length(𝛉)
    𝐈ₗ = eye(l)
    𝐍 = fill(0.0,(l,l))
    𝐌 = fill(0.0,(l,l))
    n = 1
    for dataₙ in zip(𝒟...)
        𝒟ₙ = dataₙ
        𝚲ₙ = Λ[n]
        𝐔ₙ = uₓ(entity,𝒟ₙ)
        ∂ₓ𝐮ₙ = ∂ₓu(entity, 𝒟ₙ)
        𝐁ₙ =  ∂ₓ𝐮ₙ * 𝚲ₙ * ∂ₓ𝐮ₙ'
        𝚺ₙ = Σₙ(entity,𝛉, 𝐁ₙ)
        𝚺ₙ⁻¹ = inv(𝚺ₙ)
        𝛈ₙ = 𝚺ₙ⁻¹ * 𝐔ₙ' * 𝛉
        𝐍 = 𝐍 + (𝛈ₙ' ⊗ 𝐈ₗ) * 𝐁ₙ * (𝛈ₙ ⊗ 𝐈ₗ)
        𝐌 = 𝐌 + 𝐔ₙ * 𝚺ₙ⁻¹ * 𝐔ₙ'
        n = n + 1
    end
    𝐗 = 𝐌 - 𝐍
end

function Σₙ(entity::FundamentalMatrix, 𝛉::Matrix, 𝐁ₙ::Matrix)
𝛉' * 𝐁ₙ * 𝛉
end
