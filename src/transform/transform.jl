function transform(entity::HomogeneousCoordinates, coordinate_system::CanonicalToHartley, ℳ::AbstractArray{T}) where T<:HomogeneousPoint
    (ℳ,𝐓) = hartley_normalization(ℳ)
end

function transform(entity::CovarianceMatrices, coordinate_system::CanonicalToHartley, Λ::Vector{T},𝒯) where T<:Matrix
Λ₂ = deepcopy(Λ)
dim, _ = size(Λ₂[1])
N = length(𝒯)
blocksize = Int8(round(dim/N))

k = [(i-1)* blocksize + 1  for i = 1:N+1]
    map!(Λ₂ , Λ₂) do 𝚲
        for n = 1:N
            i = k[n]
            j = k[n+1]-1
            𝚲ₙ = vcat(hcat(𝚲[i:j,i:j],zeros(blocksize,1)),zeros(1,blocksize+1))
            𝚲ₙ = 𝒯[n] * 𝚲ₙ * 𝒯[n]'
            𝚲[i:j,i:j] = 𝚲ₙ[1:blocksize,1:blocksize]
        end
        𝚲
    end
Λ₂
end
