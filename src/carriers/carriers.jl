# Carrier vector for fundamental matrix estimation.
function ∂ₓu(entity::FundamentalMatrix, 𝒟)
    m, mʹ = 𝒟
    𝐦  = 𝑛(collect(Float64,m.coords))
    𝐦ʹ = 𝑛(collect(Float64,mʹ.coords))
    ∂ₓu(entity, 𝐦 , 𝐦ʹ)
end

function ∂ₓu(entity::FundamentalMatrix, 𝐦::Vector{Float64}, 𝐦ʹ::Vector{Float64})
    𝐞₁ = [1.0 0.0 0.0]'
    𝐞₂ = [0.0 1.0 0.0]'
    [(𝐞₁ ⊗ 𝐦ʹ) (𝐞₂ ⊗ 𝐦ʹ) (𝐦 ⊗ 𝐞₁) (𝐦 ⊗ 𝐞₂)]
end


function uₓ(entity::FundamentalMatrix, 𝒟)
    m, mʹ = 𝒟
    𝐦  = 𝑛(collect(Float64,m.coords))
    𝐦ʹ = 𝑛(collect(Float64,mʹ.coords))
    uₓ(entity, 𝐦 , 𝐦ʹ)
end

function uₓ(entity::FundamentalMatrix, 𝐦::Vector{Float64}, 𝐦ʹ::Vector{Float64})
    𝐦 ⊗ 𝐦ʹ
end
