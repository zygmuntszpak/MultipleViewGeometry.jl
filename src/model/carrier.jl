function construct_carrier(task::HomographyEstimationTask, 𝐦::AbstractVector, 𝐦ʹ::AbstractVector)
    𝐔 = -𝐦 ⊗ vec2antisym(𝐦ʹ)
    return 𝐔
end

function construct_carrier(task::FundamentalMatrixEstimationTask, 𝐦::AbstractVector, 𝐦ʹ::AbstractVector)
    𝐔 = 𝐦 ⊗ 𝐦ʹ
    return 𝐔
end
