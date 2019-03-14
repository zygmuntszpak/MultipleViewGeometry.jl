
function cost(c::CostFunction, entity::FundamentalMatrix, 𝛉::AbstractArray, 𝒞::Tuple{AbstractArray, Vararg{AbstractArray}}, 𝒟::Tuple{AbstractArray, Vararg{AbstractArray}})
    ℳ, ℳʹ = 𝒟
    Λ₁, Λ₂ = 𝒞
    Jₐₘₗ = 0.0
    N = length(𝒟[1])
    𝚲ₙ = @MMatrix zeros(4,4)
    𝐞₁ = @SVector [1.0, 0.0, 0.0]
    𝐞₂ = @SVector [0.0, 1.0, 0.0]
    index = SVector(1,2)
    @inbounds for n = 1:N
        𝚲ₙ[1:2,1:2] .=  Λ₁[n][index,index]
        𝚲ₙ[3:4,3:4] .=  Λ₂[n][index,index]
        𝐦 = hom(ℳ[n])
        𝐦ʹ= hom(ℳʹ[n])
        𝐔ₙ = (𝐦 ⊗ 𝐦ʹ)
        ∂ₓ𝐮ₙ =  [(𝐞₁ ⊗ 𝐦ʹ) (𝐞₂ ⊗ 𝐦ʹ) (𝐦 ⊗ 𝐞₁) (𝐦 ⊗ 𝐞₂)]
        𝐁ₙ =  ∂ₓ𝐮ₙ * 𝚲ₙ * ∂ₓ𝐮ₙ'
        𝚺ₙ = 𝛉' * 𝐁ₙ * 𝛉
        𝚺ₙ⁻¹ = inv(𝚺ₙ)
        Jₐₘₗ +=  𝛉' * 𝐔ₙ * 𝚺ₙ⁻¹ * 𝐔ₙ' * 𝛉
    end
    Jₐₘₗ
end

# function datum(c::CostFunction, entity::FundamentalMatrix,𝒞::Tuple{AbstractArray, Vararg{AbstractArray}}, 𝒟::Tuple{AbstractArray, Vararg{AbstractArray}}, n::Integer, 𝚲ₙ::AbstractArray)
#     Λ, Λʹ = collect(𝒞)
#     ℳ, ℳʹ = collect(𝒟)
#     #𝚲 = Λ[n]
#     #𝚲ʹ = Λʹ[n]
#     #@typeof 𝚲
#     dim = 2
#     #dim, _ = size(𝚲)
#     dim = dim - 1
#     #𝚲ₙ =  @SMatrix [𝚲[1:dim,1:dim] zeros(dim,dim); zeros(dim,dim) 𝚲ʹ[1:dim,1:dim]] #TODO SMatrix?
#     #𝚲ₙ = sparse(zeros(4,4))
#     #𝚲ₙ = eye(4)
#     #𝒟ₙ = (ℳ[n], ℳʹ[n])
#     #𝒟ₙ
# end

function ∂cost(c::CostFunction, entity::FundamentalMatrix, 𝛉::Matrix, Λ::Vector{T1}, matches...) where T1 <: Matrix
    𝐗 = X(c, entity, 𝛉, Λ, matches)
    2*𝐗*𝛉
end


function covariance_matrix(c::CostFunction, s::HessianApproximation, entity::FundamentalMatrix, 𝛉::AbstractArray, 𝒞::Tuple{AbstractArray, Vararg{AbstractArray}}, 𝒟::Tuple{AbstractArray, Vararg{AbstractArray}})
    𝛉 = SVector{9}(𝛉 / norm(𝛉))
    ℳ, ℳʹ = 𝒟
    Λ₀, Λ₀ʹ = 𝒞
    N = length(ℳ)
    if (N != length(ℳʹ))
          throw(ArgumentError("There should be an equal number of points for each view."))
    end
    if (N != length(Λ₀) || N != length(Λ₀ʹ) )
          throw(ArgumentError("There should be a covariance matrix for each point correspondence."))
    end
    # Map corresponding points to the normalized coordinate system.
    𝒪, 𝐓 = transform(HomogeneousCoordinates(),CanonicalToHartley(),ℳ)
    𝒪ʹ, 𝐓ʹ = transform(HomogeneousCoordinates(),CanonicalToHartley(),ℳʹ)
    # Map estimate to the normalized coordinate system.
    𝛉₁ = (inv(𝐓') ⊗ inv(𝐓ʹ')) * 𝛉
    𝛉₁ =  𝛉₁ / norm(𝛉₁)
    # Map covariance matrices to the normalized coordinate system.
    Λ₁ = transform(CovarianceMatrices(), CanonicalToHartley(), Λ₀ , 𝐓)
    Λ₁ʹ = transform(CovarianceMatrices(), CanonicalToHartley(), Λ₀ʹ , 𝐓ʹ)

    𝐇 = _H(c, entity, 𝛉₁, (Λ₁,Λ₁ʹ), (𝒪 , 𝒪ʹ)) * 0.5 # Magic half

    # Rank-8 constrained Moore-Pensore pseudo inverse.
    d = length(𝛉)
    U,S,V = svd(𝐇)
    S = SizedArray{Tuple{9}}(S)
    for i = 1:d-1
        S[i] = 1/S[i]
    end
    S[d] = 0.0
    𝐇⁻¹ = U*diagm(S)*V'

    𝐏 = (1/norm(𝛉₁)) * ( Matrix{Float64}(I, 9, 9) - ((𝛉₁*𝛉₁') / norm(𝛉₁)^2) )
    𝚲 = 𝐏 * 𝐇⁻¹ * 𝐏


    # Derivative of the determinant of 𝚯 = reshape(𝛉₁,(3,3)).
    φ₁ = 𝛉₁[5]*𝛉₁[9] - 𝛉₁[8]*𝛉₁[6]
    φ₂ = -(𝛉₁[4]*𝛉₁[5] - 𝛉₁[7]*𝛉₁[6])
    φ₃ = 𝛉₁[4]*𝛉₁[8] - 𝛉₁[7]*𝛉₁[5]
    φ₄ = -(𝛉₁[2]*𝛉₁[9] - 𝛉₁[8]*𝛉₁[3])
    φ₅ = 𝛉₁[1]*𝛉₁[9] - 𝛉₁[7]*𝛉₁[3]
    φ₆ = -(𝛉₁[1]*𝛉₁[8] - 𝛉₁[7]*𝛉₁[2])
    φ₇ = 𝛉₁[2]*𝛉₁[6] - 𝛉₁[5]*𝛉₁[3]
    φ₈ = -(𝛉₁[1]*𝛉₁[6] - 𝛉₁[4]*𝛉₁[3])
    φ₉ = 𝛉₁[1]*𝛉₁[5] - 𝛉₁[4]*𝛉₁[2]
    ∂𝛟 = [φ₁; φ₂; φ₃; φ₄; φ₅; φ₆; φ₇; φ₈; φ₉]

    A = [ Matrix{Float64}(I, 9, 9) ; zeros(1,9)]
    B = [ Matrix{Float64}(I, 9, 9) ∂𝛟; ∂𝛟' 0]
    𝚲 = inv(B)*A*𝚲*A'*inv(B)
    𝚲 = 𝚲[1:9,1:9]

    𝛉₀ = (𝐓' ⊗ 𝐓ʹ') * 𝛉₁
    # Jacobian of the unit normalisation transformation: 𝛉 / norm(𝛉)
    ∂𝛉= (1/norm(𝛉₀)) * ( Matrix{Float64}(I, 9, 9) - ((𝛉₀*𝛉₀') / norm(𝛉₀)^2) )
    F = ∂𝛉*(𝐓' ⊗ 𝐓ʹ')
    F * 𝚲 * F'
end

function covariance_matrix(c::CostFunction, s::CanonicalApproximation, entity::FundamentalMatrix, 𝛉::AbstractArray, 𝒞::Tuple{AbstractArray, Vararg{AbstractArray}}, 𝒟::Tuple{AbstractArray, Vararg{AbstractArray}})
    𝛉 = SVector{9}(𝛉 / norm(𝛉))
    ℳ, ℳʹ = 𝒟
    Λ₀, Λ₀ʹ = 𝒞
    N = length(ℳ)
    if (N != length(ℳʹ))
          throw(ArgumentError("There should be an equal number of points for each view."))
    end
    if (N != length(Λ₀) || N != length(Λ₀ʹ) )
          throw(ArgumentError("There should be a covariance matrix for each point correspondence."))
    end
    # Map corresponding points to the normalized coordinate system.
    𝒪, 𝐓 = transform(HomogeneousCoordinates(),CanonicalToHartley(),ℳ)
    𝒪ʹ, 𝐓ʹ = transform(HomogeneousCoordinates(),CanonicalToHartley(),ℳʹ)
    # Map estimate to the normalized coordinate system.
    𝛉₁ = (inv(𝐓') ⊗ inv(𝐓ʹ')) * 𝛉
    # Map covariance matrices to the normalized coordinate system.
    Λ₁ = transform(CovarianceMatrices(), CanonicalToHartley(), Λ₀ , 𝐓)
    Λ₁ʹ = transform(CovarianceMatrices(), CanonicalToHartley(), Λ₀ʹ , 𝐓ʹ)

    𝚲  = _covariance_matrix(AML(),FundamentalMatrix(), 𝛉₁, (Λ₁,Λ₁ʹ), (𝒪 , 𝒪ʹ))

    𝛉₁ =  𝛉₁ / norm(𝛉₁)

    # Derivative of the determinant of 𝚯 = reshape(𝛉₁,(3,3)).
    φ₁ = 𝛉₁[5]*𝛉₁[9] - 𝛉₁[8]*𝛉₁[6]
    φ₂ = -(𝛉₁[4]*𝛉₁[5] - 𝛉₁[7]*𝛉₁[6])
    φ₃ = 𝛉₁[4]*𝛉₁[8] - 𝛉₁[7]*𝛉₁[5]
    φ₄ = -(𝛉₁[2]*𝛉₁[9] - 𝛉₁[8]*𝛉₁[3])
    φ₅ = 𝛉₁[1]*𝛉₁[9] - 𝛉₁[7]*𝛉₁[3]
    φ₆ = -(𝛉₁[1]*𝛉₁[8] - 𝛉₁[7]*𝛉₁[2])
    φ₇ = 𝛉₁[2]*𝛉₁[6] - 𝛉₁[5]*𝛉₁[3]
    φ₈ = -(𝛉₁[1]*𝛉₁[6] - 𝛉₁[4]*𝛉₁[3])
    φ₉ = 𝛉₁[1]*𝛉₁[5] - 𝛉₁[4]*𝛉₁[2]
    ∂𝛟 = [φ₁; φ₂; φ₃; φ₄; φ₅; φ₆; φ₇; φ₈; φ₉]

    A = [Matrix{Float64}(I, 9, 9) ; zeros(1,9)]
    B = [Matrix{Float64}(I, 9, 9) ∂𝛟; ∂𝛟' 0]
    𝚲 = inv(B)*A*𝚲*A'*inv(B)
    𝚲 = 𝚲[1:9,1:9]

    𝛉₀ = (𝐓' ⊗ 𝐓ʹ') * 𝛉₁

    # Jacobian of the unit normalisation transformation: 𝛉 / norm(𝛉)
    ∂𝛉= (1/norm(𝛉₀)) * (Matrix{Float64}(I, 9, 9) - ((𝛉₀*𝛉₀') / norm(𝛉₀)^2) )
    F = ∂𝛉*(𝐓' ⊗ 𝐓ʹ')
    F * 𝚲 * F'
end

function _covariance_matrix(c::CostFunction, entity::FundamentalMatrix, 𝛉::AbstractArray, 𝒞::Tuple{AbstractArray, Vararg{AbstractArray}}, 𝒟::Tuple{AbstractArray, Vararg{AbstractArray}})
    𝛉 = 𝛉 / norm(𝛉)
    ℳ, ℳʹ = 𝒟
    Λ₁, Λ₂ = 𝒞
    N = length(𝒟[1])
    𝚲ₙ = @MMatrix zeros(4,4)
    𝐞₁ = @SMatrix [1.0; 0.0; 0.0]
    𝐞₂ = @SMatrix [0.0; 1.0; 0.0]
    index = SVector(1,2)
    𝐌 = fill(0.0,(9,9))
    for n = 1:N
        𝚲ₙ[1:2,1:2] .=  Λ₁[n][index,index]
        𝚲ₙ[3:4,3:4] .=  Λ₂[n][index,index]
        𝐦 = ℳ[n]
        𝐦ʹ= ℳʹ[n]
        𝐔ₙ = (𝐦 ⊗ 𝐦ʹ)
        𝐀 = 𝐔ₙ*𝐔ₙ'
        ∂ₓ𝐮ₙ =  [(𝐞₁ ⊗ 𝐦ʹ) (𝐞₂ ⊗ 𝐦ʹ) (𝐦 ⊗ 𝐞₁) (𝐦 ⊗ 𝐞₂)]
        𝐁ₙ =  ∂ₓ𝐮ₙ * 𝚲ₙ * ∂ₓ𝐮ₙ'
        𝐌 = 𝐌 + 𝐀/(𝛉'*𝐁ₙ*𝛉)
    end
    d = length(𝛉)
    𝐏 = Matrix{Float64}(I, d, d) - norm(𝛉)^-2 * (𝛉*𝛉')
    U,S,V = svd(𝐌)
    S = SizedArray{Tuple{9}}(S)
    for i = 1:d-1
        S[i] = 1/S[i]
    end
    S[d] = 0.0
    𝐌⁻¹ = U*diagm(S)*V'
    𝐏 * 𝐌⁻¹ * 𝐏
end


function X(c::CostFunction, entity::ProjectiveEntity, 𝛉::AbstractArray,𝒞::Tuple{AbstractArray, Vararg{AbstractArray}}, 𝒟::Tuple{AbstractArray, Vararg{AbstractArray}})
    ℳ, ℳʹ = 𝒟
    Λ, Λʹ = 𝒞
    N = length(ℳ)
    if (N != length(ℳʹ))
          throw(ArgumentError("There should be an equal number of points for each view."))
    end

    if (N != length(Λ) || N != length(Λʹ))
          throw(ArgumentError("There should be a covariance matrix for each point correspondence."))
    end

    _X(c, entity, 𝛉, 𝒞, 𝒟)

end


# function _X(c::CostFunction, entity::ProjectiveEntity, 𝛉::AbstractArray,𝒞::Tuple{AbstractArray, Vararg{AbstractArray}}, 𝒟::Tuple{AbstractArray, Vararg{AbstractArray}})
#     l = length(𝛉)
#     𝐈ₗ = eye(l)
#     𝐍 = fill(0.0,(l,l))
#     𝐌 = fill(0.0,(l,l))
#     N = length(𝒟[1])
#     for n = 1:N
#         𝒟ₙ, 𝚲ₙ = datum(c,entity,𝒞, 𝒟, n)
#         𝐔ₙ = uₓ(entity,𝒟ₙ)
#         ∂ₓ𝐮ₙ = ∂ₓu(entity, 𝒟ₙ)
#         𝐁ₙ =  ∂ₓ𝐮ₙ * 𝚲ₙ * ∂ₓ𝐮ₙ'
#         𝚺ₙ = Σₙ(entity,𝛉, 𝐁ₙ)
#         𝚺ₙ⁻¹ = inv(𝚺ₙ)
#         𝛈ₙ = 𝚺ₙ⁻¹ * 𝐔ₙ' * 𝛉
#         𝐍 = 𝐍 + (𝛈ₙ' ⊗ 𝐈ₗ) * 𝐁ₙ * (𝛈ₙ ⊗ 𝐈ₗ)
#         𝐌 = 𝐌 + 𝐔ₙ * 𝚺ₙ⁻¹ * 𝐔ₙ'
#     end
#     𝐗 = 𝐌 - 𝐍
# end

function _X(c::CostFunction, entity::FundamentalMatrix, 𝛉::AbstractArray,𝒞::Tuple{AbstractArray, Vararg{AbstractArray}}, 𝒟::Tuple{AbstractArray, Vararg{AbstractArray}})
    l = length(𝛉)
    𝐈ₗ = SMatrix{l,l}(1.0I)
    𝐍 = @SMatrix zeros(l,l)
    𝐌 = @SMatrix zeros(l,l)
    N = length(𝒟[1])
    ℳ, ℳʹ = 𝒟
    Λ₁, Λ₂ = 𝒞
    𝚲ₙ = @MMatrix zeros(4,4)
    𝐞₁ = @SMatrix [1.0; 0.0; 0.0]
    𝐞₂ = @SMatrix [0.0; 1.0; 0.0]
    @inbounds for n = 1:N
        index = SVector(1,2)
        𝚲ₙ[1:2,1:2] .=  Λ₁[n][index,index]
        𝚲ₙ[3:4,3:4] .=  Λ₂[n][index,index]
        𝐦 = hom(ℳ[n])
        𝐦ʹ= hom(ℳʹ[n])
        𝐔ₙ = (𝐦 ⊗ 𝐦ʹ)
        ∂ₓ𝐮ₙ =  [(𝐞₁ ⊗ 𝐦ʹ) (𝐞₂ ⊗ 𝐦ʹ) (𝐦 ⊗ 𝐞₁) (𝐦 ⊗ 𝐞₂)]
        𝐁ₙ =  ∂ₓ𝐮ₙ * 𝚲ₙ * ∂ₓ𝐮ₙ'
        𝚺ₙ = 𝛉' * 𝐁ₙ * 𝛉
        𝚺ₙ⁻¹ = inv(𝚺ₙ)
        𝛈ₙ = 𝚺ₙ⁻¹ * 𝐔ₙ' * 𝛉
        𝐍 = 𝐍 + ((𝛈ₙ' ⊗ 𝐈ₗ) * 𝐁ₙ * (𝛈ₙ ⊗ 𝐈ₗ))
        𝐌 = 𝐌 + (𝐔ₙ * 𝚺ₙ⁻¹ * 𝐔ₙ')
    end
    𝐗 = 𝐌 - 𝐍
end

function Σₙ(entity::FundamentalMatrix, 𝛉::AbstractArray, 𝐁ₙ::AbstractArray)
    𝛉' * 𝐁ₙ * 𝛉
end

function H(c::CostFunction, entity::FundamentalMatrix,  𝛉::AbstractArray,𝒞::Tuple{AbstractArray, Vararg{AbstractArray}}, 𝒟::Tuple{AbstractArray, Vararg{AbstractArray}})
    ℳ, ℳʹ = 𝒟
    Λ, Λʹ = 𝒞
    N = length(ℳ)
    if (N != length(ℳʹ))
          throw(ArgumentError("There should be an equal number of points for each view."))
    end

    if (N != length(Λ) || N != length(Λʹ))
          throw(ArgumentError("There should be a covariance matrix for each point correspondence."))
    end
    _H(c, entity, 𝛉, 𝒞, 𝒟)
end


function _H(c::CostFunction, entity::ProjectiveEntity, 𝛉::AbstractArray, 𝒞::Tuple{AbstractArray, Vararg{AbstractArray}}, 𝒟::Tuple{AbstractArray, Vararg{AbstractArray}})
    𝐗 = X(c, entity, 𝛉, 𝒞, 𝒟)
    𝐓 = T(c, entity, 𝛉, 𝒞, 𝒟)

    𝐇 = 2*(𝐗-𝐓)
end


function T(c::CostFunction, entity::ProjectiveEntity, 𝛉::AbstractArray, 𝒞::Tuple{AbstractArray, Vararg{AbstractArray}}, 𝒟::Tuple{AbstractArray, Vararg{AbstractArray}})
    l = length(𝛉)
    𝐈ₗ = SMatrix{l,l}(1.0I)
    𝐈ₘ =  Iₘ(entity)
    𝐍 = @SMatrix zeros(l,l)
    𝐌 = @SMatrix zeros(l,l)
    𝐓 = @SMatrix zeros(l,l)
    N = length(𝒟[1])
    ℳ, ℳʹ = 𝒟
    Λ₁, Λ₂ = 𝒞
    𝚲ₙ = @MMatrix zeros(4,4)
    𝐞₁ = @SMatrix [1.0; 0.0; 0.0]
    𝐞₂ = @SMatrix [0.0; 1.0; 0.0]
    for n = 1: N
        index = SVector(1,2)
        𝚲ₙ[1:2,1:2] .=  Λ₁[n][index,index]
        𝚲ₙ[3:4,3:4] .=  Λ₂[n][index,index]
        𝐦 = hom(ℳ[n])
        𝐦ʹ= hom(ℳʹ[n])
        𝐔ₙ = (𝐦 ⊗ 𝐦ʹ)
        ∂ₓ𝐮ₙ =  [(𝐞₁ ⊗ 𝐦ʹ) (𝐞₂ ⊗ 𝐦ʹ) (𝐦 ⊗ 𝐞₁) (𝐦 ⊗ 𝐞₂)]
        𝐁ₙ = ∂ₓ𝐮ₙ * 𝚲ₙ * ∂ₓ𝐮ₙ'
        𝚺ₙ = 𝛉' * 𝐁ₙ * 𝛉
        𝚺ₙ⁻¹ = inv(𝚺ₙ)
        𝐓₁ = @SMatrix zeros(Float64,l,l)
        𝐓₂ = @SMatrix zeros(Float64,l,l)
        𝐓₃ = @SMatrix zeros(Float64,l,l)
        𝐓₄ = @SMatrix zeros(Float64,l,l)
        𝐓₅ = @SMatrix zeros(Float64,l,l)
        # The additional parentheses around some of the terms are needed as
        # a workaround to a bug where Base.afoldl allocates memory unnecessarily.
        # https://github.com/JuliaArrays/StaticArrays.jl/issues/537
        for k = 1:l
            𝐞ₖ = 𝐈ₗ[:,k]
            ∂𝐞ₖ𝚺ₙ = (𝐈ₘ ⊗ 𝐞ₖ') * 𝐁ₙ * (𝐈ₘ ⊗ 𝛉) + (𝐈ₘ ⊗ 𝛉') * 𝐁ₙ * (𝐈ₘ ⊗ 𝐞ₖ)
            𝐓₁ = 𝐓₁ + (((𝐔ₙ * 𝚺ₙ⁻¹) * (∂𝐞ₖ𝚺ₙ)) * 𝚺ₙ⁻¹) * 𝐔ₙ' * 𝛉 * 𝐞ₖ'
            𝐓₂ = 𝐓₂ + (𝐞ₖ' * 𝐔ₙ * 𝚺ₙ⁻¹ ⊗ 𝐈ₗ) * 𝐁ₙ * (𝚺ₙ⁻¹ * 𝐔ₙ' * 𝛉 ⊗ 𝐈ₗ) * 𝛉 * 𝐞ₖ'
            𝐓₄ = 𝐓₄ + (𝛉' * 𝐔ₙ * 𝚺ₙ⁻¹ * (∂𝐞ₖ𝚺ₙ) * 𝚺ₙ⁻¹ ⊗ 𝐈ₗ) * 𝐁ₙ * (𝚺ₙ⁻¹ * 𝐔ₙ' * 𝛉 ⊗ 𝐈ₗ) * 𝛉 * 𝐞ₖ'
            𝐓₅ = 𝐓₅ + (𝛉' * 𝐔ₙ * 𝚺ₙ⁻¹ ⊗ 𝐈ₗ) * 𝐁ₙ * (𝚺ₙ⁻¹ * (∂𝐞ₖ𝚺ₙ) * 𝚺ₙ⁻¹ * 𝐔ₙ' * 𝛉 ⊗ 𝐈ₗ) * 𝛉 * 𝐞ₖ'
        end
        𝐓₃ =  (𝛉' * 𝐔ₙ * 𝚺ₙ⁻¹ ⊗ 𝐈ₗ) * 𝐁ₙ * (𝐈ₘ ⊗ 𝛉) * 𝚺ₙ⁻¹ * 𝐔ₙ'
        𝐓 = 𝐓 + 𝐓₁ + 𝐓₂ + 𝐓₃ - 𝐓₄ - 𝐓₅
    end
    𝐓
end

@inline function Iₘ(entity::FundamentalMatrix)
     SMatrix{1,1}(1.0I)
end
