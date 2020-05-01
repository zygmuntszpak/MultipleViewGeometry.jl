
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

function fit_homography(observations::AbstractObservations, method₀::LevenbergMarquardt)
    task = HomographyEstimationTask()
    @unpack data = observations
    # Determine total number of 2D points.
    N = length(first(data))
    # Initialise the residual vector so that it need not be recreated for
    # each iteration of the LevenbergMarquardt optimization step.
    residuals = zeros(Float64, 2*N)
    objective = SumOfSquares(task, VectorValuedObjective(task, residuals))

    # TODO Instantiate the Jacobian matrix so that it need not be recreated
    # for each step of the optimization loop.
    jacobian_matrix = JacobianMatrix(task, objective, zeros(Float64, 4*N, 9+2*N))
    # 𝐉 = zeros(4*N,9+2*N)
    method @set method₀ =

    # method has a field `apply_normalization` which determines whether or not
    # to apply Hartley data normalization to the observations
    𝛉, λ = method(objective, observations) # TODO change order to: observations, objective
    𝐇 = reshape(𝛉,(3,3))
    𝐇 = SMatrix{3,3,Float64,9}(𝐇 / norm(𝐇))
    return HomographyMatrix(𝐇)
end

function (objective::VectorValuedObjective{T})(𝛉::AbstractVector, observations::AbstractObservations) where T <: HomographyEstimationTask
    @unpack residuals = objective
    @unpack data = observations
    ℳ = data[1]
    ℳ′ = data[2]

    # Nine parameters for the projection matrix, and 2 parameters per 2D point.
    N = Int((length(𝛉) - 9) / 2)
    index₁ = SVector(1,2)
    index₂ = SVector(3,4)
    reprojections = reshape(residuals,(4,N))
    𝛉v = @view 𝛉[1:9]
    𝐇 = SMatrix{3,3,Float64,9}(reshape(𝛉v,(3,3)))
    i = 10
    for n = 1:N
        # Extract the observed image points in the first and second view.
        𝐦 = ℳ[n]
        𝐦′ = ℳ′[n]
        # Extract 2D point and convert to homogeneous coordinates
        𝐦ᶿ = hom(SVector{2,Float64}(𝛉[i],𝛉[i+1]))
        reprojections[index₁,n] = hom⁻¹(𝐦ᶿ) - 𝐦
        reprojections[index₂,n] = hom⁻¹(𝐇 * 𝐦ᶿ) - 𝐦′
        i = i + 2
    end
    return residuals
end

function (jac::JacobianMatrix{T₁, T₂, T₃})(𝛉::AbstractVector, observations::AbstractObservations) where T₁ <: HomographyEstimationTask, T₂ <: SOS, T₃ <: AbstractMatrix
    @unpack jacobian = jac
    # Nine parameters for the homography matrix, and 2 parameters per 2D point.
    N = Int((length(𝛉) - 9) / 2)
    index₁ = SVector(1,2)
    index₂ = SVector(3,4)
    𝛉v = @view 𝛉[1:9]
    𝐇 = SMatrix{3,3,Float64,9}(reshape(𝛉v,(3,3)))

    # Create a view of the jacobian matrix 𝐉 and reshape it so that
    # it will be more convenient to index into the appropriate entries
    # whilst looping over all of the data points.
    𝐉 = reshape(reinterpret(Float64, jacobian), 4, N, 9+2*N)
    𝐀 = SMatrix{2,3,Float64,6}(1,0,0,1,0,0)
    𝐈₃ = SMatrix{3,3}(1.0I)
    i = 10
    for n = 1:N
        # Extract 3D point and convert to homogeneous coordinates.
        𝐦 = hom(SVector{2,Float64}(𝛉[i], 𝛉[i+1]))

        # Derivative of residual in first and second image w.r.t 2D point in the
        # first image.
        ∂𝐫₁_d𝐦 = 𝐀 * 𝐈₃
        ∂𝐫₂_d𝐦 = 𝐀 * ∂hom⁻¹(𝐇 * 𝐦) * 𝐇

        # Derivative of residual in second image w.r.t homography martix.
        # ∂𝐫₁_d𝐇 is the zero vector.
        ∂𝐫₂_d𝐇 = 𝐀 * ∂hom⁻¹(𝐇  * 𝐦) * (𝐦' ⊗ 𝐈₃)
        𝐉v[index₂,n,1:9] = ∂𝐫₂_d𝐇
        𝐉v[index₁,n,i:i+1] = ∂𝐫₁_d𝐦[:,1:2]
        𝐉v[index₂,n,i:i+1] = ∂𝐫₂_d𝐦[:,1:2]
        i = i + 2
    end
    return 𝐉
end


#fit(data, covariance_matrix; objective = AlgebraicLeastSquares(HomographyMatrixTask()))
