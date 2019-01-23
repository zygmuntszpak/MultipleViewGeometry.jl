
function project(e::Pinhole, 𝐏::AbstractArray, 𝒳::Vector{<:AbstractArray})

    if size(𝐏) != (3,4)
        throw(ArgumentError("Expect 3 x 4 projection matrix."))
    end
    ℳ = map(𝒳) do 𝐗
        𝐦 = hom⁻¹(𝐏 * hom(𝐗))
    end
end
