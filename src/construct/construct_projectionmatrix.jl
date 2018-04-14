function construct( e::ProjectionMatrix,
                   𝐊::AbstractArray{T,2},
                   𝐑::AbstractArray{T,2},
                    𝐭::AbstractArray{T,1} ) where T<:Real

    if size(𝐊) != (3,3) || size(𝐑) != (3,3)
        throw(ArgumentError("Expect 3 x 3 calibration and rotation matrices."))
    end
    if length(𝐭) != 3
        throw(ArgumentError("Expect length-3 translation vectors."))
    end
    𝐏 = 𝐊*[𝐑 -𝐑*𝐭]
end

function construct( e::ProjectionMatrix, 𝐅::AbstractArray{T,2}) where T<:Real
    𝐞 = collect(epipole(𝐅').coords)
    eye(3,4), [vec2antisym(𝐞) * 𝐅  𝐞]
end
