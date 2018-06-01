function transform(entity::HomogeneousCoordinates, coordinate_system::CanonicalToHartley, ℳ::Vector{<:AbstractArray})
    𝒪, 𝐓 = hartley_normalization(ℳ)
end

function transform(entity::CovarianceMatrices, coordinate_system::CanonicalToHartley, Λ::Vector{<:AbstractArray}, 𝐓::AbstractArray)

    Λ₂ = map(Λ) do 𝚲ₙ
         (𝐓 * 𝚲ₙ * 𝐓')
    end

end
