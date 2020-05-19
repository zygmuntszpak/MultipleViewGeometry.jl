#=
    Create an instance of an AbstractRectificationAlgorithm
    AbstractRectification
=#
abstract type AbstractRectificationAlgorithm end
struct FusielloCalibratedRectification <: AbstractRectificationAlgorithm end

"""
    Rectify(camera₁::AbstractCamera, camera₂::AbstractCamera, algorithm::AbstractRectificationAlgorithm, 𝐇₁::AbstractMatrix, 𝐇₂::AbstractMatrix)


The idea behind rectification is to define a new pair of cameras by rotating the
old camera around their optical centers until focal planes become coplanar thereby
containing the baseline [X]. This ensures that the epipoles are at infinity,
hence epipolar lines are parallel. To have horizontal eipolar lines, the baseline
must be parallel to the new x-axis of both cameras. In addition, to have a proper
rectification, point correspondences must have the same vertical coordinate.
This is obtained by requiring that the new cameras have the same intrinsic
parameters.
"""
struct Rectify{T₁ <: AbstractCamera, T₂ <: AbstractCamera, T₃ <: AbstractRectificationAlgorithm, T₄ <: AbstractMatrix, T₅ <: AbstractMatrix} <: AbstractContext
    camera₁::T₁
    camera₂::T₂
    algorithm::T₃
    𝐇₁::T₄
    𝐇₂::T₅
end

"""
    Rectify(camera₁::AbstractCamera, camera₂::AbstractCamera, algorithm::AbstractRectificationAlgorithm)

TODO

"""
# function Rectify(camera₁::AbstractCamera, camera₂::AbstractCamera, algorithm::AbstractRectificationAlgorithm)
#     rectified_camera₁, rectified_camera₂, 𝐇₁, 𝐇₂ = algorithm(camera₁, camera₂)
#     return Rectify(rectified_camera₁, rectified_camera₂, algorithm, 𝐇₁, 𝐇₂)
# end
function Rectify(camera₁::AbstractCamera, camera₂::AbstractCamera, algorithm::AbstractRectificationAlgorithm)
 rectified_camera₁, rectified_camera₂, 𝐇₁, 𝐇₂ = algorithm(camera₁, camera₂)
 return Rectify(rectified_camera₁, rectified_camera₂, algorithm, 𝐇₁, 𝐇₂)
end


#Rectify(camera₁::AbstractCamera, camera₂::AbstractCamera, algorithm::AbstractRectificationAlgorithm) = Rectify(algorithm(camera₁, camera₂)..., algorithm)

# function Rectify(camera₁::AbstractCamera, camera₂::AbstractCamera, algorithm::AbstractRectificationAlgorithm)
# end

# TODO Add functor method for rectifying images and 2D points

function (algorithm::FusielloCalibratedRectification)(camera₁::AbstractCamera, camera₂::AbstractCamera)
    𝐏₁ = matrix(Projection(camera₁))
    𝐏₂ = matrix(Projection(camera₂))

    intrinsics₁ = intrinsics(camera₁)
    intrinsics₂ = intrinsics(camera₂)
    focal_length = (intrinsics₁.focal_length + intrinsics₂.focal_length) / 2
    skewedness = 0.0
    scale_x = (intrinsics₁.scale_x + intrinsics₂.scale_x) / 2
    scale_y = (intrinsics₁.scale_y + intrinsics₂.scale_y) / 2
    # TODO determine the width and height based on the minimum rectangle that
    # contains the rectified image coordinates.
    width =  (intrinsics₁.width + intrinsics₂.width) / 2
    height = (intrinsics₁.height + intrinsics₂.height) / 2
    principal_point = (intrinsics₁.principal_point + intrinsics₂.principal_point) / 2
    common_intrinsics = IntrinsicParameters(focal_length = focal_length,
                                     skewedness = skewedness,
                                     scale_x = scale_x,
                                     scale_y = scale_y,
                                     width = width,
                                     height = height,
                                     principal_point = principal_point,
                                     coordinate_system = OpticalSystem())
    extrinsics₁ = extrinsics(camera₁)
    extrinsics₂ = extrinsics(camera₂)
    𝐨₁ = origin(extrinsics₁)
    𝐨₂ = origin(extrinsics₂)
    𝐞₁, 𝐞₂, 𝐞₃ = basis_vectors(extrinsics₁)
    # The new x-axis points in the direction of the baseline
    𝐞₁′ = Vec(𝐨₁ - 𝐨₂) / norm(Vec(𝐨₁ - 𝐨₂))
    # The new y-axis is orthogonal to the new x-axis and the old z-axis
    𝐞₂′ = cross(𝐞₃, 𝐞₁′) / norm(cross(𝐞₃, 𝐞₁′))
    # The new z-axis is orthogonal to the baseline and the y-axis
    𝐞₃′ = cross(𝐞₁′, 𝐞₂′) / norm(cross(𝐞₁′, 𝐞₂′))
    rectified_extrinsics₁ = ExtrinsicParameters(coordinate_system = CartesianSystem(𝐨₁, 𝐞₁′, 𝐞₂′, 𝐞₃′))
    rectified_extrinsics₂ = ExtrinsicParameters(coordinate_system = CartesianSystem(𝐨₂, 𝐞₁′, 𝐞₂′, 𝐞₃′))

    rectified_camera₁ = deepcopy(camera₁)
    rectified_camera₁ = @set camera₁.model.intrinsics = common_intrinsics
    rectified_camera₁ = @set camera₁.model.extrinsics = rectified_extrinsics₁

    rectified_camera₂ = deepcopy(camera₂)
    rectified_camera₂ = @set camera₂.model.intrinsics = common_intrinsics
    rectified_camera₂ = @set camera₂.model.extrinsics = rectified_extrinsics₂

    𝐏₁′ = matrix(Projection(rectified_camera₁))
    𝐏₂′ = matrix(Projection(rectified_camera₂))

    𝐇₁ = SMatrix{3,3, Float64}(𝐏₁′[1:3, 1:3] * inv(𝐏₁[1:3, 1:3]))
    𝐇₂ = SMatrix{3,3, Float64}(𝐏₂′[1:3, 1:3] * inv(𝐏₂[1:3, 1:3]))
    # 𝐇₁ = SMatrix{3,3, Float64}(𝐏₁[1:3, 1:3] * inv(𝐏₁′[1:3, 1:3]))
    # 𝐇₂ = SMatrix{3,3, Float64}(𝐏₂[1:3, 1:3] * inv(𝐏₂′[1:3, 1:3]))

    return rectified_camera₁, rectified_camera₂, 𝐇₁, 𝐇₂
end
