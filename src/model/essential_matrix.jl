struct EssentialMatrix{T₁ <: AbstractMatrix} <: ProjectiveEntity
    𝐄::T₁
end

function matrix(entity::EssentialMatrix)
    entity.𝐄
end

EssentialMatrix(camera₁::AbstractCamera, camera₂::AbstractCamera) = EssentialMatrix(camera₁, camera₂, CartesianSystem(Point(0.0, 0.0, 0.0),Vec(1.0, 0.0, 0.0), Vec(0.0, 1.0, 0.0), Vec(0.0, 0.0, 1.0)))
EssentialMatrix(camera₁::AbstractCamera, camera₂::AbstractCamera, reference_system::AbstractCoordinateSystem) = EssentialMatrix(construct_essential_matrix(camera₁, camera₂, reference_system))


#EssentialMatrix(P₁::Projection, P₂::Projection) = EssentialMatrix(construct_essential_matrix(P₁, P₂))
# TODO remove reference_system
function construct_essential_matrix(camera₁::AbstractCamera, camera₂::AbstractCamera,  reference_system::AbstractCoordinateSystem)
    image_type₁ = image_type(camera₁)
    image_system₁ = coordinate_system(image_type₁)
    image_type₂ = image_type(camera₂)
    image_system₂ = coordinate_system(image_type₂)
    intrinsics₁ = intrinsics(camera₁)
    intrinsics₂ = intrinsics(camera₂)
    𝐊₁ = matrix(intrinsics₁, image_system₁)
    𝐊₂ = matrix(intrinsics₂, image_system₂)
    𝐅 = matrix(FundamentalMatrix(camera₁, camera₂))
    𝐄 = 𝐊₂'*𝐅*𝐊₁
    return 𝐄 / norm(𝐄)
end
