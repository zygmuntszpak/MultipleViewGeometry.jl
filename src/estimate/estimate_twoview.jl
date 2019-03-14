include("fundamental_matrix.jl")
include("homography_matrix.jl")

# 𝒟::Tuple{AbstractArray, Vararg{AbstractArray}}
function estimate(entity::ProjectiveEntity, method::ManualEstimation,  𝒟::Tuple{AbstractArray, Vararg{AbstractArray}})
    method.𝚹
end

# # estimate(entity::HomographyMatrix, method::DirectLinearTransform, 𝓓::Tuple{Vector{T₁} where T₁ <: AbstractArray, Vector{T₂} where T₂ <: AbstractArray})
# function estimate(entity::ProjectiveEntity, method::ManualEstimation, 𝓓::Tuple{Vector{T₁} where T₁ <: AbstractArray, Vector{T₂} where T₂ <: AbstractArray})
#     # map((ℳ, ℳʹ) -> estimate(entity, method, (𝓜, 𝓜ʹ)), 𝓓)
#     @show "Here"
#     𝓜, 𝓜ʹ =  𝓓
#     𝓡 = Vector{typeof{entity}}(undef,length(𝓓))
#     for k = 1:length(𝓓)
#         𝓡[k] = estimate(entity, method, (𝓜[k], 𝓜ʹ[k]))
#     end
#     𝓡
# end
