struct CalibrateCamera <: AbstractContext end

# TODO world should be restricted to a "CalibrationWorld" type
function (calibrate::CalibrateCamera)(world::AbstractWorld, cameras::Vector{<:AbstractCamera})
    @unpack points = world
    # Drop the z-coordinates from the points on the calibration plane since they
    # are zero anyway. The aim will be to compute a homography between the
    # points in ℳ′ and the projections of the 3D points in each camera view.
    ℳ′ = [Point(p[1], p[2]) for p in points]

    aquire = AquireImage()
    # Determine projections of the 3D points in each camera view.
    𝓜 = [aquire(world, camera) for camera in cameras]

    # Estimate a homography matrix between the points on the calibration plane
    # and the points in each image.
    ℋ = [fit_homography(ℳ′, ℳ, DirectLinearTransform()) for ℳ in 𝓜]

    # TODO Refine the homography matrix estimate by minimising the gold-standard
    # reprojection error.
    𝐀₁ = get_camera_intrinsics(ℋ)
    display(𝐀₁)
    𝐀₂ = get_camera_intrinsics(ℋ; use_analytical_method = false)
    display(𝐀₂)
    return nothing
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

function get_camera_extrinsics(𝐀::AbstractArray, ℋ::Vector{<: HomographyMatrix})

end
