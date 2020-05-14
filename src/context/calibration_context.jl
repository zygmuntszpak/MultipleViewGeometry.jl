struct CalibrateCamera <: AbstractContext end

# TODO world should be restricted to a "CalibrationWorld" type
function (calibrate::CalibrateCamera)(world::AbstractWorld, cameras::Vector{<:AbstractCamera})
    @unpack points = world
    @unpack coordinate_system = world

    aquire = AquireImage()
    # Determine projections of the 3D points in each camera view.
    𝓜 = [aquire(world, camera) for camera in cameras]
    return calibrate(world, 𝓜 )


    # # Drop the z-coordinates from the points on the calibration plane since they
    # # are zero anyway. The aim will be to compute a homography between the
    # # points in ℳ′ and the projections of the 3D points in each camera view.
    # ℳ′ = [Point(p[1], p[2]) for p in points]
    #
    # aquire = AquireImage()
    # # Determine projections of the 3D points in each camera view.
    # 𝓜 = [aquire(world, camera) for camera in cameras]
    #
    # # Estimate a homography matrix between the points on the calibration plane
    # # and the points in each image.
    # ℋ = [fit_homography(ℳ′, ℳ, DirectLinearTransform()) for ℳ in 𝓜]
    #
    # # TODO Refine the homography matrix estimate by minimising the gold-standard
    # # reprojection error.
    #
    # # Determine intrinsics parameters from the set of homographies.
    # 𝐀 = get_camera_intrinsics(ℋ; use_analytical_method = false)
    #
    # # TODO There is still an unresolved ambiguity which ought to be enforced
    # # which has to do with whether the recovered extrinsic parameters locate
    # # the camera in front of or behind the calibration grid.
    # ℰ = get_camera_extrinsics(𝐀, ℋ)
    #
    # # Determine the lens distortion parameters.
    # 𝐤 = get_lens_distortion(𝐀, ℰ, points, 𝓜)
    #
    # # Refine all estimates by minimising the reprojection error.
    # cameras = fit_sole_camera_rig(𝐀, 𝐤, ℰ, points, 𝓜, LevenbergMarquardt())
    #
    # return cameras
end

function (calibrate::CalibrateCamera)(world::AbstractWorld, 𝓜::AbstractVector{<: AbstractVector})
    @unpack points = world
    #@unpack coordinate_system = world

    # Drop the z-coordinates from the points on the calibration plane since they
    # are zero anyway. The aim will be to compute a homography between the
    # points in ℳ′ and the projections of the 3D points in each camera view.
    ℳ′ = [Point(p[1], p[2]) for p in points]

    # Estimate a homography matrix between the points on the calibration plane
    # and the points in each image.
    ℋ = [fit_homography(ℳ′, ℳ, DirectLinearTransform()) for ℳ in 𝓜]

    # TODO Refine the homography matrix estimate by minimising the gold-standard
    # reprojection error.

    # Determine intrinsics parameters from the set of homographies.
    𝐀 = get_camera_intrinsics(ℋ; use_analytical_method = false)

    # TODO There is still an unresolved ambiguity which ought to be enforced
    # which has to do with whether the recovered extrinsic parameters locate
    # the camera in front of or behind the calibration grid.
    ℰ = get_camera_extrinsics(𝐀, ℋ)

    # Determine the lens distortion parameters.
    𝐤 = get_lens_distortion(𝐀, ℰ, points, 𝓜)

    # Refine all estimates by minimising the reprojection error.
    cameras = fit_sole_camera_rig(𝐀, 𝐤, ℰ, points, 𝓜, LevenbergMarquardt())

    return cameras
end

"""
    get_camera_intrinsics(ℋ::Vector{<: HomographyMatrix})

Takes a sequence of homography matrices and returns the common intrinsic
camera matrix 𝐀.
"""
function get_camera_intrinsics(ℋ::Vector{<: HomographyMatrix}; use_analytical_method::Bool = true)
    M = length(ℋ)
    𝐕 = zeros(2*M, 6)
    i = 1
    for m = 1:M
        𝐇 = matrix(ℋ[m])
        𝐕[i,:] .= helper(1, 2, 𝐇)
        𝐕[i + 1,:] .= helper(1, 1, 𝐇) - helper(2, 2, 𝐇)
        i = i + 2
    end
    # Solve 𝐕 * 𝐛 = 0
    F = svd(𝐕)
    𝐛 = F.Vt[end,:]
    B₀, B₁, B₂, B₃, B₄, B₅ = 𝐛
    𝐀 = use_analytical_method ? obtain_intrinsics_analytically(𝐛) : obtain_intrinsics_numerically(𝐛)
    return 𝐀
end

function obtain_intrinsics_analytically(𝐛::AbstractVector)
    B₀, B₁, B₂, B₃, B₄, B₅ = 𝐛
    ω = B₀*B₂*B₅ - B₁^2*B₅ - B₀*B₄^2 + 2*B₁*B₃*B₄ - B₂*B₃^2
    d = B₀*B₂ - B₁^2
    α = sqrt(max(ω / (d*B₀), 0.0))
    β = sqrt(max((ω / d^2) * B₀, 0.0))
    γ = sqrt(max(ω /  (d^2*B₀), 0.0)) * B₁
    uc = (B₁*B₄ - B₂*B₃) / d
    vc = (B₁*B₃ - B₀*B₄) / d
    # Intrinsic parameter matrix
    𝐀 = SMatrix{3,3,Float64,9}(α, 0, 0, γ, β, 0, uc, vc, 1)
    return 𝐀
end

function obtain_intrinsics_numerically(𝐛::AbstractVector)
    B₀, B₁, B₂, B₃, B₄, B₅ = 𝐛
    𝐁 = [B₀ B₁ B₃;
         B₁ B₂ B₄;
         B₃ B₄ B₅]
    # Make sure 𝐁 is positive definite
    if B₀ < 0 || B₂ < 0 || B₅ < 0
        𝐁 = -𝐁
    end
    𝐋 = cholesky(Symmetric(𝐁)).L
    # Intrinsic parameter matrix
    𝐀 = SMatrix{3,3,Float64,9}(inv(𝐋)') * 𝐋[3,3]
    return 𝐀
end

function helper(s::Integer, t::Integer, 𝐇::AbstractArray)
    H₁ₛ = 𝐇[1,s]
    H₁ₜ = 𝐇[1,t]
    H₂ₛ = 𝐇[2,s]
    H₂ₜ = 𝐇[2,t]
    H₃ₛ = 𝐇[3,s]
    H₃ₜ = 𝐇[3,t]
    𝐯 = [H₁ₛ * H₁ₜ,
         H₁ₛ * H₂ₜ + H₂ₛ * H₁ₜ,
         H₂ₛ * H₂ₜ,
         H₃ₛ * H₁ₜ + H₁ₛ * H₃ₜ,
         H₃ₛ * H₂ₜ + H₂ₛ * H₃ₜ,
         H₃ₛ * H₃ₜ]
    return 𝐯
end

"""
    get_camera_extrinsics(𝐀::AbstractArray, ℋ::Vector{<: HomographyMatrix})

Takes a matrix 𝐀 representing the intrinsic camera parameters together
with a sequence of homography matrices and returns a sequence of extrinsic
view parameters ℰ = (𝐄₁, ..., 𝐄ₘ) with 𝐄ᵢ = [𝐑ᵢ, 𝐭ᵢ].
"""
function get_camera_extrinsics(𝐀::AbstractArray, ℋ::Vector{<: HomographyMatrix})
    ℰ = [estimate_view_transform(𝐀, H) for H in ℋ]
    return ℰ
end

function estimate_view_transform(𝐀::AbstractArray, H::HomographyMatrix)
    𝐇 = matrix(H)
    𝐇 = 𝐇 / 𝐇[3,3] # TODO follow up on this convention
    𝐡₁ = 𝐇[:,1]
    𝐡₂ = 𝐇[:,2]
    𝐡₃ = 𝐇[:,3]
    𝐀⁻¹ = inv(𝐀)
    κ = 1 / norm(𝐀⁻¹ * 𝐡₁)
    𝐫₁ = κ * (𝐀⁻¹ * 𝐡₁)
    𝐫₂ = κ * (𝐀⁻¹ * 𝐡₂)
    𝐫₃ = cross(𝐫₁, 𝐫₂)
    𝐑₀ = hcat(𝐫₁, 𝐫₂, 𝐫₃)
    𝐑 = make_true_rotation_matrix(𝐑₀)
    𝐭 = SVector{3, Float64}(κ * (𝐀⁻¹ * 𝐡₃))
    return hcat(𝐑, 𝐭)
end

function make_true_rotation_matrix(𝐑₀::AbstractArray)
    F = svd(𝐑₀)
    # Zhang suggests 𝐑 = 𝐔 * 𝐕', but including 𝐒 accomodates for reflections.
    𝐒 = Diagonal([1.0 1.0 det(F.U * F.Vt)])
    𝐑 = SMatrix{3,3,Float64,9}(F.U * 𝐒 * F.Vt)
    return 𝐑
end

"""
    get_lens_distortion(𝐀::AbstractArray, ℰ::AbstractArray, ℳ′::AbstractArray, 𝓜::AbstractArray)

Takes a matrix 𝐀 representing the intrinsic camera parameters; the estimated
extrinsics parameters ℰ = (𝐄₁, ..., 𝐄ₘ) with 𝐄ᵢ = [𝐑ᵢ, 𝐭ᵢ]; the target model
points 𝒳 and the observed sensor point 𝓜 = [ℳ₁, ..., ℳₘ] with
ℳᵢ = [𝐦₁, ... ,𝐦ₙ] being the points for view i. Returns the vector 𝐤
of estimated lens distortion coefficients.
"""
function get_lens_distortion(𝐀::AbstractArray, ℰ::AbstractArray, 𝒳::AbstractArray, 𝓜::AbstractArray)
    # The number of views.
    M = length(ℰ)
    # The number of model points.
    N = length(𝒳)
    # The projection center (in sensor coordinates).
    uc = 𝐀[1,3]
    vc = 𝐀[2,3]

    𝐃 = zeros(2*M*N, 2)
    𝐝 = zeros(2*M*N, 1)
    l = 0
    i = 1
    for m = 1:M
        𝐄ₘ = ℰ[m]
        #=
         Construct a 'canonical projection' matrix by assuming a focal length
         of 1, with principal point at (0,0) such that the resulting intrinsic
         matrix is identity.
        =#
        𝐏₀  = 𝐄ₘ
        #=
            Construct the actual projection matrix based on the given intrinsics
            and extrinsics.
        =#
        𝐏  = 𝐀 * 𝐄ₘ
        ℳ = 𝓜[m]
        for n = 1:N
            𝐗ₙ = 𝒳[n]
            # Canonical projection.
            𝐦₀ = hom⁻¹(𝐏₀ * hom(𝐗ₙ))
            # Radius in the canonical projection coordinates.
            r = norm(𝐦₀)
            # Projection with actual camera.
            𝐦₁ = hom⁻¹(𝐏 * hom(𝐗ₙ))
            u₁, v₁ = 𝐦₁
            δu = u₁ - uc
            δv = v₁ - vc
            𝐃[i, :] .= [δu*r^2, δu*r^4]
            𝐃[i + 1, :] .= [δv*r^2, δv*r^4]
            # observed image point
            𝐦₂ = ℳ[n]
            u₂, v₂ = 𝐦₂
            𝐝[i] = u₂ - u₁
            𝐝[i + 1] = v₂ - v₁
            i = i + 2
        end
    end
    # Solve the linear system of equations to obtain lens parameter vector.
    𝐤 = 𝐃 \ 𝐝
    return SVector{2,Float64}(𝐤...)
end

"""
   refine_parameters(𝐀::AbstractArray, 𝛋::AbstractVector, ℰ::AbstractVector, 𝒳::AbstractVector, 𝓜::AbstractVector)

Given a seed for the camera instrincs 𝐀, lens distortion 𝛋, extrinsic
view parameters ℰ, 3D points on the calibration target 𝒳 and observed
image points 𝓜, returns refined estimates for the camera intrinsics,
distortion parameters and camera view parameters, respectively.
"""
function refine_parameters(𝐀::AbstractArray, 𝐤::AbstractVector, ℰ::AbstractVector, 𝒳::AbstractVector, 𝓜::AbstractVector)
    𝛈 = compose_parameter_vector(𝐀, 𝐤, ℰ)
    return 𝛈
end
