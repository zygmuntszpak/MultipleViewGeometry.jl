function cost(c::CostFunction, entity::HomographyMatrix, 𝛉::AbstractArray, 𝒞::Tuple{AbstractArray, Vararg{AbstractArray}}, 𝒟::Tuple{AbstractArray, Vararg{AbstractArray}})
    ℳ, ℳʹ = 𝒟
    Λ₁, Λ₂ = 𝒞
    Jₐₘₗ = 0.0
    N = length(𝒟[1])
    𝚲ₙ = @MMatrix zeros(4,4)
    𝐞₁ = @SVector [1.0, 0.0, 0.0]
    𝐞₂ = @SVector [0.0, 1.0, 0.0]
    𝐞₁ₓ = vec2antisym(𝐞₁)
    𝐞₂ₓ = vec2antisym(𝐞₂)
    𝐈₃₂ = @SMatrix [1.0  0.0 ; 0.0 1.0 ; 0.0 0.0]
    𝐈₂ = @SMatrix  [1.0  0.0 ; 0.0 1.0]
    index = SVector(1,2)
    @inbounds for n = 1:N
        𝚲ₙ[1:2,1:2] .=  Λ₁[n][index,index]
        𝚲ₙ[3:4,3:4] .=  Λ₂[n][index,index]
        𝐦 = hom(ℳ[n])
        𝐦ʹ= hom(ℳʹ[n])
        𝐦ʹₓ = vec2antisym(𝐦ʹ)
        𝐔ₙ = (-𝐦 ⊗ 𝐦ʹₓ)
        𝐕ₙ = 𝐔ₙ * 𝐈₃₂
        ∂ₓ𝐯ₙ = -hcat(vec((𝐞₁ ⊗ 𝐦ʹₓ)*𝐈₃₂), vec((𝐞₂ ⊗ 𝐦ʹₓ)*𝐈₃₂), vec((𝐦 ⊗ 𝐞₁ₓ)*𝐈₃₂), vec((𝐦 ⊗ 𝐞₂ₓ)*𝐈₃₂))
        𝐁ₙ =  ∂ₓ𝐯ₙ * 𝚲ₙ * ∂ₓ𝐯ₙ'
        𝚺ₙ = (𝐈₂ ⊗ 𝛉') * 𝐁ₙ * (𝐈₂ ⊗ 𝛉)
        𝚺ₙ⁻¹ = inv(𝚺ₙ)
        Jₐₘₗ +=  𝛉' * 𝐕ₙ * 𝚺ₙ⁻¹ * 𝐕ₙ' * 𝛉
    end
    Jₐₘₗ
end

function _X(c::CostFunction, entity::HomographyMatrix, 𝛉::AbstractArray,𝒞::Tuple{AbstractArray, Vararg{AbstractArray}}, 𝒟::Tuple{AbstractArray, Vararg{AbstractArray}})
    𝐈₉ = SMatrix{9,9}(1.0I)
    𝐈₃₂ = @SMatrix [1.0  0.0 ; 0.0 1.0 ; 0.0 0.0]
    𝐈₂ = @SMatrix  [1.0  0.0 ; 0.0 1.0]
    𝐍 = @SMatrix zeros(9,9)
    𝐌 = @SMatrix zeros(9,9)
    N = length(𝒟[1])
    ℳ, ℳʹ = 𝒟
    Λ₁, Λ₂ = 𝒞
    𝚲ₙ = @MMatrix zeros(4,4)
    𝐞₁ = @SMatrix [1.0; 0.0; 0.0]
    𝐞₂ = @SMatrix [0.0; 1.0; 0.0]
    𝐞₁ₓ = vec2antisym(𝐞₁)
    𝐞₂ₓ = vec2antisym(𝐞₂)
    @inbounds for n = 1:N
        index = SVector(1,2)
        𝚲ₙ[1:2,1:2] .=  Λ₁[n][index,index]
        𝚲ₙ[3:4,3:4] .=  Λ₂[n][index,index]
        𝐦 = hom(ℳ[n])
        𝐦ʹ= hom(ℳʹ[n])
        𝐦ʹₓ = vec2antisym(𝐦ʹ)
        𝐔ₙ = -𝐦 ⊗ 𝐦ʹₓ
        𝐕ₙ = 𝐔ₙ * 𝐈₃₂
        ∂ₓ𝐯ₙ = -hcat(vec((𝐞₁ ⊗ 𝐦ʹₓ)*𝐈₃₂), vec((𝐞₂ ⊗ 𝐦ʹₓ)*𝐈₃₂), vec((𝐦 ⊗ 𝐞₁ₓ)*𝐈₃₂), vec((𝐦 ⊗ 𝐞₂ₓ)*𝐈₃₂))
        𝐁ₙ =  ∂ₓ𝐯ₙ * 𝚲ₙ * ∂ₓ𝐯ₙ'
        𝚺ₙ = (𝐈₂ ⊗ 𝛉') * 𝐁ₙ * (𝐈₂ ⊗ 𝛉)
        𝚺ₙ⁻¹ = inv(𝚺ₙ)
        𝛈ₙ = 𝚺ₙ⁻¹ * 𝐕ₙ' * 𝛉
        𝐍 = 𝐍 + ((𝛈ₙ' ⊗ 𝐈₉) * 𝐁ₙ * (𝛈ₙ ⊗ 𝐈₉))
        𝐌 = 𝐌 + (𝐕ₙ * 𝚺ₙ⁻¹ * 𝐕ₙ')
    end
    𝐗 = 𝐌 - 𝐍
end


function covariance_matrix(c::CostFunction, s::CanonicalApproximation, entity::HomographyMatrix, 𝛉::AbstractArray, 𝒞::Tuple{AbstractArray, Vararg{AbstractArray}}, 𝒟::Tuple{AbstractArray, Vararg{AbstractArray}})
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
    𝛉₁ = (inv(𝐓') ⊗ 𝐓ʹ) * 𝛉
    𝛉₁ =  𝛉₁ / norm(𝛉₁)
    # Map covariance matrices to the normalized coordinate system.
    Λ₁ = transform(CovarianceMatrices(), CanonicalToHartley(), Λ₀ , 𝐓)
    Λ₁ʹ = transform(CovarianceMatrices(), CanonicalToHartley(), Λ₀ʹ , 𝐓ʹ)

    𝚲  = _covariance_matrix(AML(),HomographyMatrix(), 𝛉₁, (Λ₁,Λ₁ʹ), (𝒪 , 𝒪ʹ))

    𝛉₀ = (𝐓' ⊗ inv(𝐓ʹ)) * 𝛉₁
    𝛉₀ = 𝛉₀ / norm(𝛉₀)

    # Jacobian of the unit normalisation transformation: 𝛉 / norm(𝛉)
    ∂𝛉= (1/norm(𝛉₀)) * (Matrix{Float64}(I, 9, 9) - ((𝛉₀*𝛉₀') / norm(𝛉₀)^2) )
    F = ∂𝛉*((𝐓' ⊗ inv(𝐓ʹ)))
    F * 𝚲 * F'
end

function _covariance_matrix(c::CostFunction, entity::HomographyMatrix, 𝛉::AbstractArray, 𝒞::Tuple{AbstractArray, Vararg{AbstractArray}}, 𝒟::Tuple{AbstractArray, Vararg{AbstractArray}})
    𝛉 = 𝛉 / norm(𝛉)
    ℳ, ℳʹ = 𝒟
    Λ₁, Λ₂ = 𝒞
    N = length(𝒟[1])
    𝐈₉ = SMatrix{9,9}(1.0I)
    𝐈₃₂ = @SMatrix [1.0  0.0 ; 0.0 1.0 ; 0.0 0.0]
    𝐈₂ = @SMatrix  [1.0  0.0 ; 0.0 1.0]
    𝚲ₙ = @MMatrix zeros(4,4)
    𝐞₁ = @SMatrix [1.0; 0.0; 0.0]
    𝐞₂ = @SMatrix [0.0; 1.0; 0.0]
    𝐞₁ₓ = vec2antisym(𝐞₁)
    𝐞₂ₓ = vec2antisym(𝐞₂)
    index = SVector(1,2)
    𝐌 = fill(0.0,(9,9))
    for n = 1:N
        𝚲ₙ[1:2,1:2] .=  Λ₁[n][index,index]
        𝚲ₙ[3:4,3:4] .=  Λ₂[n][index,index]
        𝐦 = ℳ[n]
        𝐦ʹ= ℳʹ[n]
        𝐦ʹₓ = vec2antisym(𝐦ʹ)
        𝐔ₙ = -𝐦 ⊗ 𝐦ʹₓ
        𝐕ₙ = 𝐔ₙ * 𝐈₃₂
        ∂ₓ𝐯ₙ = -hcat(vec((𝐞₁ ⊗ 𝐦ʹₓ)*𝐈₃₂), vec((𝐞₂ ⊗ 𝐦ʹₓ)*𝐈₃₂), vec((𝐦 ⊗ 𝐞₁ₓ)*𝐈₃₂), vec((𝐦 ⊗ 𝐞₂ₓ)*𝐈₃₂))
        𝐁ₙ =  ∂ₓ𝐯ₙ * 𝚲ₙ * ∂ₓ𝐯ₙ'
        𝚺ₙ = (𝐈₂ ⊗ 𝛉') * 𝐁ₙ * (𝐈₂ ⊗ 𝛉)
        𝚺ₙ⁻¹ = inv(𝚺ₙ)
        𝐌 = 𝐌 + (𝐕ₙ * 𝚺ₙ⁻¹ * 𝐕ₙ')
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
