function construct( e::FundamentalMatrix,
                   𝐊₁::AbstractArray,
                   𝐑₁::AbstractArray,
                   𝐭₁::AbstractArray,
                   𝐊₂::AbstractArray,
                   𝐑₂::AbstractArray,
                   𝐭₂::AbstractArray)

    if size(𝐊₁) != (3,3) || size(𝐊₂) != (3,3) ||
       size(𝐑₁) != (3,3) || size(𝐑₂) != (3,3)
        throw(ArgumentError("Expect 3 x 3 calibration and rotation matrices."))
    end
    if length(𝐭₁) != 3 || length(𝐭₂) != 3
        throw(ArgumentError("Expect length-3 translation vectors."))
    end
    𝐅 = vec2antisym(𝐊₂*𝐑₂*(𝐭₁ .- 𝐭₂))*𝐊₂*𝐑₂/𝐑₁/𝐊₁
    SMatrix{3,3,Float64,3*3}(𝐅)
end

function construct( e::FundamentalMatrix, 𝐏₁::AbstractArray, 𝐏₂::AbstractArray)
    if (size(𝐏₁) != (3,4)) || (size(𝐏₂) != (3,4))
        throw(ArgumentError("Expect 3 x 4 projection matrices."))
    end
    𝐜₁ = SVector{4,Float64}(nullspace(Array(𝐏₁)))
    𝐞₂ = 𝐏₂*𝐜₁
    𝐅 = vec2antisym(𝐞₂)*𝐏₂*pinv(Array(𝐏₁))
    SMatrix{3,3,Float64,3*3}(𝐅)
end
