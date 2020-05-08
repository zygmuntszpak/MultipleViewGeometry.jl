# TODO take a camera model as input
function fit_camera_matrix(𝐀::AbstractArray, 𝐤::AbstractVector, ℰ::AbstractVector, 𝒳::AbstractVector, 𝓜::AbstractVector, method::LevenbergMarquardt)
    task = CameraCalibrationTask()
    # The total number of views.
    M = length(ℰ)
    # Camera intrinsics, lens distortion and extrinsics.
    𝛈 = compose_parameter_vector(𝐀, 𝐤, ℰ)
    # The total number of image points across all views.
    N = sum(map(x-> length(x), 𝓜))
    # The projections of the 3D points onto each image, and the actual 3D points.
    observations = Observations(tuple(tuple(𝓜, 𝒳)))

    # Initialise the residual vector so that it need not be recreated for
    # each iteration of the LevenbergMarquardt optimization step.
    residuals = zeros(Float64, 2*N*M)
    objective = SumOfSquares(task, VectorValuedObjective(task, residuals))


    # TODO Instantiate the Jacobian matrix so that it need not be recreated
    # for each step of the optimization loop.
    jacobian_matrix = JacobianMatrix(task, objective, observations, zeros(Float64, 2*N, length(𝛈)))

    method = @set method.seed = ManualEstimation(𝛈) # TODO check/resolve type instability

    # TODO There needs to be an explicit test for the veracity of the Jacobian matrix.
    𝐉 = jacobian_matrix(𝛈)
    @unpack vector_valued_objective = objective
    g = x-> vector_valued_objective(x, observations)
    z = g(𝛈)
    𝐉₂ = FiniteDiff.finite_difference_jacobian(g, 𝛈)
    println("The first")
    display(𝐉)
    # println("The second")
    display(𝐉₂)
    # println("The end")
    display(norm(𝐉 - 𝐉₂))

    𝛈, λ = method(objective, observations, jacobian_matrix) # TODO change order to: observations, objective
    @show λ
    # 𝐇 = reshape(𝛈[1:9],(3,3))
    # 𝐇 = SMatrix{3,3,Float64,9}(𝐇 / norm(𝐇))
    # return HomographyMatrix(𝐇)
    return nothing
end

function (objective::VectorValuedObjective{T})(𝛉::AbstractVector, observations::AbstractObservations) where T <: CameraCalibrationTask
    #@unpack residuals = objective # TODO See "Find workaround"
    @unpack data = observations
    𝓜, 𝒳  = first(data) # TODO reconsider how data is stored in Observations struct.
    M = length(𝓜)
    N = length(𝒳)
    index₁ = SVector(1,2)
    # Assume that all 3D points are projected onto all images.
    residuals = zeros(2*N*M)   # TODO Find workaround. At the moment we recreate the residual vector so that we can use finite_difference_jacobian for verification
    projection_error = reshape(residuals,(2,N*M))
    # Camera intrinsic parameters.
    𝐢 = SVector{5, Float64}(𝛉[1], 𝛉[2], 𝛉[3], 𝛉[4], 𝛉[5])
    # Lens distortion parameters.
    𝐤 = SVector{2, Float64}(𝛉[6], 𝛉[7])
    # Keep track of the extrinsics for each view.
    i = 8
    # Keep track of the residual for each data point in each view.
    j = 0
    for m = 1:M
        # The set of image points in the mth view.
        ℳ = 𝓜[m]
        # Extrinsic parameters for the mth view.
        𝛚 = SVector{6, Float64}(𝛉[i], 𝛉[i+1], 𝛉[i+2], 𝛉[i+3], 𝛉[i+4], 𝛉[i+5])
        for n = 1:N
            𝐦 = ℳ[n]
            𝐗 = 𝒳[n]
            j = j + 1
            projection_error[index₁, j] = 𝐦 - project_with_lens(𝐗, 𝐢, 𝐤, 𝛚)
        end
        i = i + 6
    end
    return residuals
end

function (jacobian_functor::JacobianMatrix{T₁, T₂, T₃})(𝛉::AbstractVector) where {T₁ <: CameraCalibrationTask, T₂ <: SumOfSquares, T₃ <: AbstractObservations, T₄ <: AbstractMatrix}
    @unpack jacobian = jacobian_functor
    @unpack observations = jacobian_functor
    @unpack data = observations

    𝓜, 𝒳  = first(data) # TODO reconsider how data is stored in Observations struct. Perhaps introduce a Correspondence type
    # The total number of views.
    M = length(𝓜)
    # The total number of 3D points.
    N = length(𝒳)
    index₁ = SVector(1,2)

    # Create a view of the jacobian matrix 𝐉 and reshape it so that
    # it will be more convenient to index into the appropriate entries
    # whilst looping over all of the data points.
    𝐉 = reshape(reinterpret(Float64, jacobian), (2, M*N, length(𝛉)))

    # Camera intrinsic parameters.
    𝐢 = SVector{5, Float64}(𝛉[1], 𝛉[2], 𝛉[3], 𝛉[4], 𝛉[5])
    # Lens distortion parameters.
    𝐤 = SVector{2, Float64}(𝛉[6], 𝛉[7])
    # Keep track of the extrinsics for each view.
    i = 8
    # Keep track of the residual for each data point in each view.
    j = 0
    for m = 1:M
        # Extrinsic parameters for the mth view.
        𝛚 = SVector{6, Float64}(𝛉[i], 𝛉[i+1], 𝛉[i+2], 𝛉[i+3], 𝛉[i+4], 𝛉[i+5])
        for n = 1:N
            𝐗 = 𝒳[n]
            j = j + 1
            𝐉₀ = -∂project_with_lens(𝐗, 𝐢, 𝐤, 𝛚)
            # 𝐢, 𝐤
            𝐉[index₁, j, 1:7] .= 𝐉₀[:,1:7]
            # 𝛚
            𝐉[index₁, j, i:i+5] = 𝐉₀[:,8:end]
        end
        i = i + 6
    end
    return jacobian
end


"""
   compose_parameter_vector(𝐀::AbstractArray, 𝛋::AbstractVector, ℰ::AbstractVector)

Given the camera instrincs 𝐀, lens distortion 𝛋, and extrinsic
view parameters ℰ, returns a length-(7 + M6) parameter vector 𝛈.
"""
function compose_parameter_vector(𝐀::AbstractArray, 𝐤::AbstractVector, ℰ::AbstractVector)
    α = 𝐀[1,1]
    γ = 𝐀[1,2]
    uc = 𝐀[1,3]
    β = 𝐀[2,2]
    vc = 𝐀[2,3]
    k₁ = 𝐤[1]
    k₂ = 𝐤[2]
    𝐚 = [α, β, γ, uc, vc, k₁,k₂]
    𝛈 = 𝐚
    M = length(ℰ)
    for m = 1:M
        𝐄ₘ = ℰ[m]
        𝐑 = RotMatrix{3}(𝐄ₘ[:, 1:3])
        # "stereographic projection" of a normalized quaternion
        𝐫 = SPQuat(𝐑)
        𝐭 = 𝐄ₘ[:, 4]
        𝐰 = [𝐫.x, 𝐫.y, 𝐫.z, 𝐭[1], 𝐭[2], 𝐭[3]]
        𝛈 = vcat(𝛈, 𝐰)
    end
    return 𝛈
end


"""
    P̃(𝐗::AbstractVector, 𝐰::AbstractVector)

Maps the 3D point 𝐗 to image coordinates (i.e. sensor coordinates  in
the `OpticalSystem` coordinate system,  before camera intrinsics are applied)
using the view parameters 𝛚 = [𝐫, 𝐭], where 𝐫 are modified Rodrigues parameters
and 𝐭 is a translation.
"""
function P̃(𝐗::AbstractVector, 𝛚::AbstractVector)
    𝐫 = SPQuat(𝛚[1], 𝛚[2], 𝛚[3])
    𝐭 = SVector(𝛚[4], 𝛚[5], 𝛚[6])
    𝐑 = RotMatrix(𝐫)
    𝐖 = hcat(𝐑, 𝐭)
    #=
        f(𝛚; 𝐗) =  hom⁻¹(𝐖 * hom(𝐗)) where 𝐖 = reshape(𝐰(𝛚), (3,4)).
        This is equivalent to  f(𝛚; 𝐗) = hom⁻¹((𝐗' ⊗ 𝐈₃) * 𝐰(𝛚)).
    =#
    𝐦 = hom⁻¹(𝐖 * hom(𝐗))
    return 𝐦
end

# TODO make this a functor
function project_with_lens(𝐗::AbstractVector, 𝛈::AbstractVector)
    # Camera intrinsic parameters.
    𝐢 = SVector{5, Float64}(𝛈[1:5]...)
    # Lens distortion parameters.
    𝐤 = SVector{2, Float64}(𝛈[6:7]...)
    # Extrinsics parameters (modified Rodrigues rotation and translation)
    𝛚 = SVector{6, Float64}(𝛈[8:13]...)
    # Project assuming the identity matrix for camera intrinsics
    𝐱₀ = P̃(𝐗, 𝛚)
    # Apply lens distortion.
    𝐱₁ = distort(𝐱₀, 𝐤)
    # Apply the affine transformation associated with the intrinsic camera parameters.
    𝐲 = apply_intrinsics(𝐱₁, 𝐢)
    return 𝐲
end

# TODO make this a functor
function project_with_lens(𝐗::AbstractVector, 𝐢::AbstractVector, 𝐤::AbstractVector, 𝛚::AbstractVector)
    # Project assuming the identity matrix for camera intrinsics
    𝐱₀ = P̃(𝐗, 𝛚)
    # Apply lens distortion.
    𝐱₁ = distort(𝐱₀, 𝐤)
    # Apply the affine transformation associated with the intrinsic camera parameters.
    𝐲 = apply_intrinsics(𝐱₁, 𝐢)
    return 𝐲
end

function ∂project_with_lens(𝐗::AbstractVector, 𝛈::AbstractVector)
    # Camera intrinsic parameters.
    𝐢 = SVector{5, Float64}(𝛈[1:5]...)
    # Lens distortion parameters.
    𝐤 = SVector{2, Float64}(𝛈[6:7]...)
    # Extrinsics parameters (modified Rodrigues rotation and translation)
    𝛚 = SVector{6, Float64}(𝛈[8:13]...)

    # Project assuming the identity matrix for camera intrinsics.
    𝐱₀ = P̃(𝐗, 𝛚)
    # Apply lens distortion.
    𝐱₁ = distort(𝐱₀, 𝐤)
    # Apply affine transformation associated with the intrinsic camera parameters.
    𝐲 = apply_intrinsics(𝐱₁, 𝐢)

    𝐉₁ = ∂𝐀₀_𝐢(𝐱₁, 𝐢)
    𝐉₂ = ∂𝐀₀_𝐱(𝐱₁, 𝐢) * ∂ₖdistort(𝐱₀, 𝐤)
    𝐉₃ = ∂𝐀₀_𝐱(𝐱₁, 𝐢) * ∂ₓdistort(𝐱₀, 𝐤) * ∂P̃_𝛚(𝐗, 𝛚)
    𝐉  = hcat(𝐉₁, 𝐉₂, 𝐉₃)

    return 𝐉
end
function ∂project_with_lens(𝐗::AbstractVector, 𝐢::AbstractVector, 𝐤::AbstractVector, 𝛚::AbstractVector)
    # Project assuming the identity matrix for camera intrinsics.
    𝐱₀ = P̃(𝐗, 𝛚)
    # Apply lens distortion.
    𝐱₁ = distort(𝐱₀, 𝐤)
    # Apply affine transformation associated with the intrinsic camera parameters.
    𝐲 = apply_intrinsics(𝐱₁, 𝐢)

    𝐉₁ = ∂𝐀₀_𝐢(𝐱₁, 𝐢)
    𝐉₂ = ∂𝐀₀_𝐱(𝐱₁, 𝐢) * ∂ₖdistort(𝐱₀, 𝐤)
    𝐉₃ = ∂𝐀₀_𝐱(𝐱₁, 𝐢) * ∂ₓdistort(𝐱₀, 𝐤) * ∂P̃_𝛚(𝐗, 𝛚)
    𝐉  = hcat(𝐉₁, 𝐉₂, 𝐉₃)

    return 𝐉
end

function apply_intrinsics(𝐱::AbstractVector, 𝐢::AbstractVector)
    𝐦 = hom(𝐱)
    α = 𝐢[1]
    γ = 𝐢[2]
    β = 𝐢[3]
    uc = 𝐢[4]
    vc = 𝐢[5]
    𝐀₀ = SMatrix{2,3,Float64,6}(α, 0 , γ , β, uc, vc)
    𝐮 = 𝐀₀ * 𝐦
    return 𝐮
end

function ∂𝐀₀_𝐱(𝐱::AbstractVector, 𝐢::AbstractVector)
    α = 𝐢[1]
    γ = 𝐢[2]
    β = 𝐢[3]
    𝐉 = SMatrix{2,2,Float64,4}(α, 0 , γ , β)
    return 𝐉
end

function ∂𝐀₀_𝛈(𝐱::AbstractVector, 𝛈::AbstractVector)
    # Camera intrinsic parameters.
    𝐢 = SVector{5, Float64}(𝛈[1:5]...)
    𝐉₁ = ∂𝐀₀_𝐢(𝐱, 𝐢)
    𝐉₂ = zeros(2,2)
    𝐉₃ = zeros(2,6)
    𝐉  = hcat(𝐉₁, 𝐉₂, 𝐉₃)
    return 𝐉
end

function ∂𝐀₀_𝐢(𝐱::AbstractVector, 𝐢::AbstractVector)
    𝐦 = hom(𝐱)
    α = 𝐢[1]
    γ = 𝐢[2]
    β = 𝐢[3]
    uc = 𝐢[4]
    vc = 𝐢[5]
    𝐃 = SMatrix{6,5,Float64, 30}(1, 0, 0, 0, 0, 0,
                                 0, 0, 1, 0, 0, 0,
                                 0, 0, 0, 1, 0, 0,
                                 0, 0, 0, 0, 1, 0,
                                 0, 0, 0, 0, 0, 1)
    𝐈₂ = SMatrix{2,2, Float64}(I(2))
    𝐉 = (𝐦' ⊗ 𝐈₂) * 𝐃
    return 𝐉
end

function ∂𝐰_𝛚(𝛚::AbstractVector)
    𝐫 = SPQuat(𝛚[1], 𝛚[2], 𝛚[3])
    𝐭 = SVector(𝛚[4], 𝛚[5], 𝛚[6])
    𝐑 = RotMatrix(𝐫)
    𝐖 = hcat(𝐑, 𝐭)
    𝐈₃ = SMatrix{3,3, Float64}(I(3))
    ∂R_𝛚 = vcat(Rotations.jacobian(RotMatrix, 𝐫), zeros(3,3))
    ∂t_𝛚 = vcat(zeros(9,3), 𝐈₃)
    𝐉 = hcat(∂R_𝛚, ∂t_𝛚)
    return 𝐉
end

function ∂P̃_𝛈(𝐗::AbstractVector, 𝛈::AbstractVector)
    # Camera intrinsic parameters.
    𝐢 = SVector{5, Float64}(𝛈[1:5]...)
    # Lens distortion parameters.
    𝐤 = SVector{2, Float64}(𝛈[6:7]...)
    # Extrinsics parameters (modified Rodrigues rotation and translation)
    𝛚 = SVector{6, Float64}(𝛈[8:13]...)

    # TODO ∂𝐗

    𝐉₁ = zeros(2,5)
    𝐉₂ = zeros(2,2)
    𝐉₃ = ∂P̃_𝛚(𝐗, 𝛚)
    𝐉  = hcat(𝐉₁, 𝐉₂, 𝐉₃)
    return 𝐉
end

function  ∂P̃_𝛚(𝐗::AbstractVector, 𝛚::AbstractVector)
    𝐫 = SPQuat(𝛚[1], 𝛚[2], 𝛚[3])
    𝐭 = SVector(𝛚[4], 𝛚[5], 𝛚[6])
    𝐑 = RotMatrix(𝐫)
    𝐖 = hcat(𝐑, 𝐭)
    𝐰 = vec(𝐖)
    𝐈₃ = SMatrix{3,3, Float64}(I(3))
    𝐃 = SMatrix{2,3,Float64,6}(1, 0, 0, 1, 0, 0)
    #=
        f(𝛚; 𝐗) =  hom⁻¹(𝐖 * hom(𝐗)) where 𝐖 = reshape(𝐰(𝛚), (3,4)).
        This is equivalent to  f(𝛚; 𝐗) = (hom(𝐗)' ⊗ 𝐈₃) * 𝐰(𝛚).
    =#
    𝐌 = hom(𝐗)
    𝐉 = 𝐃 * ∂hom⁻¹((𝐌' ⊗ 𝐈₃) * 𝐰) * (𝐌' ⊗ 𝐈₃) * ∂𝐰_𝛚(𝛚)
    return 𝐉
end


# TODO Make this a functor
function distort(𝐱::AbstractVector, 𝐤::AbstractVector)
    D = radial_deviation
    𝐲 = 𝐱 + 𝐱 * D(norm(𝐱), 𝐤)
    return 𝐲
end

function ∂ₖdistort(𝐱::AbstractVector, 𝐤::AbstractVector)
    𝐝ᵤ = 𝐱 * ∂ₖradial_deviation(norm(𝐱), 𝐤)
    return 𝐝ᵤ
end

function ∂ₓdistort(𝐱::AbstractVector, 𝐤::AbstractVector)
    𝐈 = SMatrix{2,2}(I)
    𝐝ᵤ =  𝐈 + 𝐈 * radial_deviation(norm(𝐱), 𝐤) + 𝐱 * ∂ᵣradial_deviation(norm(𝐱), 𝐤) *  ∂norm(𝐱)
    return 𝐝ᵤ
end

# TODO Make this a functor
function radial_deviation(r::Number, 𝐤::AbstractVector)
    return 𝐤[1] * r^2 + 𝐤[2] * r^4
end

function ∂ᵣradial_deviation(r::Number, 𝐤::AbstractVector)
    return 2*𝐤[1]*r + 4*𝐤[2]*r^3
end

function ∂ₖradial_deviation(r::Number, 𝐤::AbstractVector)
    ∂₁ = r^2
    ∂₂ = r^4
    return transpose(SVector(∂₁, ∂₂))
end

function ∂norm(𝐱)
    return 𝐱' / norm(𝐱)
end
