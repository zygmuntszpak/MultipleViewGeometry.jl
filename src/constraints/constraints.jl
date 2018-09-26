abstract type Constraint end

mutable struct EpipolarConstraint <: Constraint
end

# Triangulation from Two Views Revisited: Hartley-Sturm vs. Optimal Correction
function satisfy(entity::FundamentalMatrix, constraint::EpipolarConstraint, 𝐅::AbstractArray, 𝒟::Tuple{AbstractArray, Vararg{AbstractArray}})
    ℳ, ℳʹ = 𝒟
    𝒪 = similar(ℳ)
    𝒪ʹ = similar(ℳʹ)

    N = length(ℳ)
    𝐏₂ = SMatrix{3,3,Float64,3^2}([1 0 0; 0 1 0; 0 0 0])

    I = 10
    for n = 1:N
        𝐦 = ℳ[n]
        𝐦ʹ = ℳʹ[n]
        𝐦ₕ = init_correction_view_1(𝐅, 𝐦, 𝐦ʹ, 𝐏₂)
        𝐦ₕʹ= init_correction_view_2(𝐅, 𝐦, 𝐦ʹ, 𝐏₂)
        for i = 1:I
            𝐦ₜ =  𝐦 - 𝐦ₕ
            𝐦ₜʹ =  𝐦ʹ - 𝐦ₕʹ
            𝐦ₕ = update_correction_view_1(𝐅, 𝐦, 𝐦ₕ, 𝐦ₜ, 𝐦ʹ, 𝐦ₕʹ, 𝐦ₜʹ, 𝐏₂)
            𝐦ₕʹ = update_correction_view_2(𝐅, 𝐦, 𝐦ₕ, 𝐦ₜ, 𝐦ʹ, 𝐦ₕʹ, 𝐦ₜʹ, 𝐏₂)
        end
        𝒪[n] = 𝐦ₕ
        𝒪ʹ[n] = 𝐦ₕʹ
    end
    𝒪 ,𝒪ʹ
end

function init_correction_view_1(𝐅::AbstractArray, 𝐦::AbstractVector, 𝐦ʹ::AbstractVector, 𝐏₂::AbstractArray)
    𝐦 -  dot(𝐦,𝐅*𝐦ʹ)*𝐏₂*𝐅*𝐦ʹ / ( dot(𝐅*𝐦ʹ,𝐏₂*𝐅*𝐦ʹ) + dot(𝐅'*𝐦, 𝐏₂*𝐅'*𝐦) )
end

function init_correction_view_2(𝐅::AbstractArray, 𝐦::AbstractVector, 𝐦ʹ::AbstractVector, 𝐏₂::AbstractArray)
    𝐦ʹ -  dot(𝐦,𝐅*𝐦ʹ)*𝐏₂*𝐅'*𝐦 / ( dot(𝐅*𝐦ʹ,𝐏₂*𝐅*𝐦ʹ) + dot(𝐅'*𝐦, 𝐏₂*𝐅'*𝐦) )
end

function update_correction_view_1(𝐅::AbstractArray, 𝐦::AbstractVector, 𝐦ₕ::AbstractVector, 𝐦ₜ::AbstractVector,  𝐦ʹ::AbstractVector, 𝐦ₕʹ::AbstractVector, 𝐦ₜʹ::AbstractVector,  𝐏₂::AbstractArray)
    𝐦 -  ( ( dot(𝐦ₕ,𝐅*𝐦ₕʹ) + dot(𝐅*𝐦ₕʹ, 𝐦ₜ) + dot(𝐅'*𝐦ₕ, 𝐦ₜʹ)  ) * 𝐏₂*𝐅*𝐦ₕʹ)   / (dot(𝐅*𝐦ₕʹ, 𝐏₂*𝐅*𝐦ₕʹ) + dot(𝐅'*𝐦ₕ, 𝐏₂*𝐅'*𝐦ₕ) )
end

function update_correction_view_2(𝐅::AbstractArray, 𝐦::AbstractVector, 𝐦ₕ::AbstractVector, 𝐦ₜ::AbstractVector,  𝐦ʹ::AbstractVector, 𝐦ₕʹ::AbstractVector, 𝐦ₜʹ::AbstractVector,  𝐏₂::AbstractArray)
    𝐦ʹ -  ( ( dot(𝐦ₕ,𝐅*𝐦ₕʹ) + dot(𝐅*𝐦ₕʹ, 𝐦ₜ) + dot(𝐅'*𝐦ₕ, 𝐦ₜʹ)  ) * 𝐏₂*𝐅'*𝐦ₕ)   / (dot(𝐅*𝐦ₕʹ, 𝐏₂*𝐅*𝐦ₕʹ) + dot(𝐅'*𝐦ₕ, 𝐏₂*𝐅'*𝐦ₕ) )
end
