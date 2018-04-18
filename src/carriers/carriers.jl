# Carrier vector for fundamental matrix estimation.
@inline function ∂ₓu(entity::FundamentalMatrix, 𝒟)
    m, mʹ = collect(𝒟)
    # 𝐦  = 𝑛(collect(Float64,m.coords))
    # 𝐦ʹ = 𝑛(collect(Float64,mʹ.coords))
    𝐦  = 𝑛(m)
    𝐦ʹ = 𝑛(mʹ)
    ∂ₓu(entity, 𝐦 , 𝐦ʹ)
end

@inline function ∂ₓu(entity::FundamentalMatrix, 𝐦::AbstractVector, 𝐦ʹ::AbstractVector)
    𝐞₁ = [1.0 0.0 0.0]'
    𝐞₂ = [0.0 1.0 0.0]'
    [(𝐞₁ ⊗ 𝐦ʹ) (𝐞₂ ⊗ 𝐦ʹ) (𝐦 ⊗ 𝐞₁) (𝐦 ⊗ 𝐞₂)]
end


@inline function uₓ(entity::FundamentalMatrix, 𝒟)
    m, mʹ = collect(𝒟)
    # 𝐦  = 𝑛(collect(Float64,m.coords))
    # 𝐦ʹ = 𝑛(collect(Float64,mʹ.coords))
    𝐦  = 𝑛(m)
    𝐦ʹ = 𝑛(mʹ)
    uₓ(entity, 𝐦 , 𝐦ʹ)
end

@inline function uₓ(entity::FundamentalMatrix, 𝐦::AbstractVector, 𝐦ʹ::AbstractVector)
    𝐦 ⊗ 𝐦ʹ
end
