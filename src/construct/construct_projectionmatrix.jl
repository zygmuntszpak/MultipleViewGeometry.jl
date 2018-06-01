function construct( e::ProjectionMatrix,
                   𝐊::AbstractArray,
                   𝐑::AbstractArray,
                    𝐭::AbstractArray )

    if size(𝐊) != (3,3) || size(𝐑) != (3,3)
        throw(ArgumentError("Expect 3 x 3 calibration and rotation matrices."))
    end
    if length(𝐭) != 3
        throw(ArgumentError("Expect length-3 translation vectors."))
    end
    𝐏 = 𝐊*[𝐑 -𝐑*𝐭]
    SMatrix{3,4,Float64,3*4}(𝐏)
end

function construct( e::ProjectionMatrix, 𝐅::AbstractArray)
    𝐞 = epipole(𝐅')
    𝐏₁ = eye(3,4)
    𝐏₂ = [vec2antisym(𝐞) * 𝐅  𝐞]

    SMatrix{3,4,Float64,3*4}(𝐏₁), SMatrix{3,4,Float64,3*4}(𝐏₂)

end
