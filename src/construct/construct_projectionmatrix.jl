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
    # TODO: Reconcile this change in convention with the rest of the code.
    #𝐏 = 𝐊*[𝐑 -𝐑*𝐭]
    𝐏 = 𝐊*[𝐑' -𝐑'*𝐭]
    SMatrix{3,4,Float64,3*4}(𝐏)
end

function construct(e::ProjectionMatrix, 𝐅::AbstractArray)
    𝐞 = epipole(𝐅')
    𝐏₁ = Matrix{Float64}(I,3,4)
    𝐏₂ = [vec2antisym(𝐞) * 𝐅  𝐞]

    SMatrix{3,4,Float64,3*4}(𝐏₁), SMatrix{3,4,Float64,3*4}(𝐏₂)

end

function construct(e::ProjectionMatrix, 𝐄::AbstractArray, 𝒟::Tuple{AbstractArray, Vararg{AbstractArray}})
    𝐖 = SMatrix{3,3,Float64,3*3}([0 -1 0; 1 0 0; 0 0 1])
    𝐙 = SMatrix{3,3,Float64,3*3}([0 1 0; -1 0 0; 0 0 0])
    𝐔,𝐒,𝐕 = svd(𝐄)
    𝐭 = 𝐔[:,3]
    𝐏₁ = SMatrix{3,4,Float64,3*4}(1.0I)
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

function construct(e::ProjectionMatrices, ℋ::Tuple{AbstractArray, Vararg{AbstractArray}})
    construct(LatentVariables(HomographyMatrices()), ℋ)
end

function construct(lv::LatentVariables, ℋ::Tuple{AbstractArray, Vararg{AbstractArray}})
        N = length(ℋ)
        if N < 2
            throw(ArgumentError("Please supply at least two homography matrices."))
        end
        𝐗₁ = ℋ[1]
        𝛍 = zeros(N)
        𝐉 = Array{Float64}(undef,(3,(N-1)*6))
        i₁ = range(1, step = 6, length = N - 1)
        i₂ = range(6, step = 6, length = N - 1)
        for n = 2:N
            𝐗ₙ = ℋ[n]
            e₁, e₂ = find_nearest_eigenvalues(eigvals(Array(𝐗₁), Array(𝐗ₙ)))
            𝐘 = hcat(e₁ * 𝐗ₙ - 𝐗₁, e₂ * 𝐗ₙ - 𝐗₁)
            μ = (e₁ + e₂) / 2
            𝛍[n] = real(μ)
            𝐉[:,i₁[n-1]:i₂[n-1]] .= 𝐘
        end
        𝛈 = initialisation_procedure(𝐉, 𝛍, ℋ)
        𝐚  = @view 𝛈[1:9]
        𝐀 = reshape(𝐚, (3,3))
        𝐛 = @view 𝛈[10:12]
        𝐏₁ = SMatrix{3,4}(1.0I)
        𝐏₂ = SMatrix{3,4}(hcat(𝐀, 𝐛))
        𝐏₁, 𝐏₂
end

function initialisation_procedure(𝐉::AbstractArray, 𝛍::AbstractArray, ℋ::Tuple{AbstractArray, Vararg{AbstractArray}})
    N = length(ℋ)
    if N < 2
        throw(ArgumentError("Please supply at least two homography matrices."))
    end
    F = svd(𝐉)
    𝛈 = zeros(9 + 3 + (N*3) + N)
    𝐛 = real(F.U[:,1])
    𝐗₁ = ℋ[1]
    𝐀 = 𝐗₁
    𝐯₁ = SVector(0,0,0)
    wₙ = 1
    # pack 𝛈 = [𝐚,𝐛, 𝐯₁,...,𝐯ₙ, w₁, ..., wₙ]
    𝛈[1:9] .= vec(𝐀)
    𝛈[10:12] .= 𝐛
    for (n,i) in enumerate(range(13, step = 3, length = N))
        if n == 1
            𝛈[i:i+2] .= 𝐯₁
        else
            𝐗ₙ = ℋ[n]
            𝛈[i:i+2] .= 𝐯₁ +  (𝛍[n] * 𝐗ₙ - 𝐗₁)' * 𝐛 / (norm(𝐛)^2)
        end
    end
    𝛈[end-(N-1):end] .= wₙ
    𝛈
end

function find_nearest_eigenvalues(e::AbstractArray)
    dist = SVector(abs(e[1]-e[2]), abs(e[1]-e[3]), abs(e[2]-e[3]))
    minval, index = findmin(dist)
    if index == 3
        i₁ = 2
        i₂ = 3
    elseif index == 2
        i₁ = 1
        i₂ = 3
    else
        i₁ = 1
        i₂ = 2
    end
    e[i₁], e[i₂]
end

function unpack(e::LatentVariables, 𝛈::AbstractArray)
    N = div(length(𝛈) - 12,  4)
    𝐚  = @view 𝛈[1:9]
    𝐀 = reshape(𝐚, (3,3))
    𝐛 = @view 𝛈[10:12]
    𝐰 = @view 𝛈[end-(N-1):end]
    r = range(13, step = 3, length = N+1)
    @show first(r), last(r), N
    𝐯 = reshape(view(𝛈,first(r):last(r)-1), (3,N))
    𝒫 = Array{Array{Float64,2},1}(undef,(N,))
    for n = 1:N
        𝒫[n] = 𝐰[n]*𝐀 + 𝐛*𝐯[:,n]'
    end
    𝒫
end
