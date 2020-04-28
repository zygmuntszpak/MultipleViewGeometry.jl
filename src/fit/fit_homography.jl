
function fit_homography(pts₁::AbstractVector, pts₂::AbstractVector, method::AbstractProjectiveOptimizationScheme)
    observations = Observations(tuple(pts₁, pts₂))
    fit_homography(observations, method)
end

function fit_homography(observations::AbstractObservations, method::DirectLinearTransform)
    task = HomographyEstimationTask()
    objective = AlgebraicLeastSquares(task)
    # method has a field `apply_normalization` which determines whether or not
    # to apply Hartley data normalization to the observations
    𝛉, λ = method(objective, observations) # TODO change order to: observations, objective
    𝐇 = reshape(𝛉,(3,3))
    𝐇 = SMatrix{3,3,Float64,9}(𝐇 / norm(𝐇))
    return HomographyMatrix(𝐇)
end

# HomographyMatrixTask
function (coordinate_transformation::HartleyNormalizeDataContext)(task::HomographyEstimationTask, direction::ToNormalizedSpace, 𝛉::AbstractVector)
    𝒯 = matrices(coordinate_transformation)
    𝐓 = 𝒯[1]
    𝐓ʹ = 𝒯[2]
    𝛉′ = (inv(𝐓') ⊗ 𝐓ʹ) * 𝛉
    return 𝛉′ / norm(𝛉′)
end

function (coordinate_transformation::HartleyNormalizeDataContext)(task::HomographyEstimationTask, direction::FromNormalizedSpace, 𝛉′::AbstractVector)
    𝒯 = matrices(coordinate_transformation)
    𝐓 = 𝒯[1]
    𝐓ʹ = 𝒯[2]
    𝛉 = (𝐓' ⊗ inv(𝐓ʹ)) * 𝛉′
    return 𝛉 / norm(𝛉)
end




#fit(data, covariance_matrix; objective = AlgebraicLeastSquares(HomographyMatrixTask()))
