
function fit_fundamental_matrix(pts₁::AbstractVector, pts₂::AbstractVector, method::AbstractProjectiveOptimizationScheme)
    observations = Observations(tuple(pts₁, pts₂))
    fit_fundamental_matrix(observations, method)
end

function fit_fundamental_matrix(observations::AbstractObservations, method::DirectLinearTransform)
    task = FundamentalMatrixEstimationTask()
    objective = AlgebraicLeastSquares(task)
    # method has a field `apply_normalization` which determines whether or not
    # to apply Hartley data normalization to the observations
    𝛉, λ = method(objective, observations) # TODO change order to: observations, objective
    #TODO decide whether to enforce in normalized or unnormalized space
    𝐅 = SMatrix{3,3,Float64,9}(reshape(𝛉,(3,3)))
    𝐅 = enforce_ranktwo(𝐅)
    return FundamentalMatrix(𝐅)
end

function enforce_ranktwo(𝐅::AbstractArray)
    # Enforce the rank-2 constraint.
    F = svd(𝐅)
    𝐒 =  Diagonal(SVector(F.S[1], F.S[2], 0.0))
    𝐅 = F.U * 𝐒 * F.Vt
    return SMatrix{3,3,Float64,9}(𝐅 / norm(𝐅))
end

# FundamentalMatrixTask
function (coordinate_transformation::HartleyNormalizeDataContext)(task::FundamentalMatrixEstimationTask, direction::ToNormalizedSpace, 𝛉::AbstractVector)
    𝒯 = matrices(coordinate_transformation)
    𝐓 = 𝒯[1]
    𝐓ʹ = 𝒯[2]
    𝛉′ = (inv(𝐓') ⊗ inv(𝐓ʹ')) * 𝛉
    return 𝛉′ / norm(𝛉′)
end

function (coordinate_transformation::HartleyNormalizeDataContext)(task::FundamentalMatrixEstimationTask, direction::FromNormalizedSpace, 𝛉′::AbstractVector)
    𝒯 = matrices(coordinate_transformation)
    𝐓 = 𝒯[1]
    𝐓ʹ = 𝒯[2]
    𝛉 = (𝐓' ⊗ 𝐓ʹ') * 𝛉′
    return 𝛉 / norm(𝛉)
end
