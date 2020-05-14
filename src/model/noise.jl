# Homogeneous and isotropic noise (the same covariance matrix for all points)
function apply_noise(ℳ::AbstractVector, 𝚲::AbstractMatrix)
    F = cholesky(Symmetric(𝚲))
    𝐕 = F.L
    dim = length(first(ℳ))
    𝒪 = map(ℳ) do 𝐦
        𝐫 = @SVector randn(dim)
        Δ = 𝐕*𝐫
        𝐦 + Δ
    end
end

function apply_noise(ℳ::AbstractVector, 𝒞::AbstractVector)
    𝒪 = similar(ℳ)
    dim = length(first(ℳ))
    for i in eachindex(ℳ)
        𝚲 = 𝒞[i]
        𝐦 = ℳ[i]
        F = cholesky(Symmetric(𝚲))
        𝐕 = F.L
        𝐫 = @SVector randn(dim)
        Δ = 𝐕*𝐫
        𝒪[i] = 𝐦 + Δ
    end
    return 𝒪
end
