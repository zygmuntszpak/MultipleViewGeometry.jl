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

function construct( e::ProjectionMatrix, 𝐄::AbstractArray, 𝒟::Tuple{AbstractArray, Vararg{AbstractArray}})
    𝐖 = SMatrix{3,3,Float64,3*3}([0 -1 0; 1 0 0; 0 0 1])
    𝐙 = SMatrix{3,3,Float64,3*3}([0 1 0; -1 0 0; 0 0 0])
    𝐔,𝐒,𝐕 = svd(𝐄)
    𝐭 = 𝐔[:,3]
    𝐏₁ = SMatrix{3,4,Float64,3*4}(eye(3,4))
    𝐏₂₁ = SMatrix{3,4,Float64,3*4}([𝐔*𝐖*𝐕'  𝐭])
    𝐏₂₂ = SMatrix{3,4,Float64,3*4}([𝐔*𝐖'*𝐕' 𝐭])
    𝐏₂₃ = SMatrix{3,4,Float64,3*4}([𝐔*𝐖*𝐕' -𝐭])
    𝐏₂₄ = SMatrix{3,4,Float64,3*4}([𝐔*𝐖'*𝐕' -𝐭])

    𝒳₁ = triangulate(DirectLinearTransform(), 𝐏₁, 𝐏₂₁, 𝒟)
    𝒳₂ = triangulate(DirectLinearTransform(), 𝐏₁, 𝐏₂₂, 𝒟)
    𝒳₃ = triangulate(DirectLinearTransform(), 𝐏₁, 𝐏₂₃, 𝒟)
    𝒳₄ = triangulate(DirectLinearTransform(), 𝐏₁, 𝐏₂₄, 𝒟)

    # Determine which projection matrix in the second view triangulated
    # the majority of points in front of the cameras.
    ℳ₁ = map(𝒳₁) do 𝐗
        𝐦 = 𝐏₂₁ * 𝐗
        𝐦[3] > 0
    end

    ℳ₂ = map(𝒳₂) do 𝐗
        𝐦 = 𝐏₂₂ * 𝐗
        𝐦[3] > 0
    end

    ℳ₃ = map(𝒳₃) do 𝐗
        𝐦 = 𝐏₂₃ * 𝐗
        𝐦[3] > 0
    end

    ℳ₄ = map(𝒳₄) do 𝐗
        𝐦 = 𝐏₂₄ * 𝐗
        𝐦[3] > 0
    end

    total, index = findmax((sum(ℳ₁), sum(ℳ₂), sum(ℳ₃), sum(ℳ₄)))

    if index == 1
        return 𝐏₁,  𝐏₂₁
    elseif index == 2
        return 𝐏₁,  𝐏₂₂
    elseif index == 3
        return 𝐏₁,  𝐏₂₃
    else
        return 𝐏₁,  𝐏₂₄
    end
end
