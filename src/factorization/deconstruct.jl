
function deconstruct(lv::LatentVariables, ℋ::Tuple{AbstractArray, Vararg{AbstractArray}})
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

# function unpack(e::LatentVariables, 𝛈::AbstractArray)
#     N = div(length(𝛈) - 12,  4)
#     𝐚  = @view 𝛈[1:9]
#     𝐀 = reshape(𝐚, (3,3))
#     𝐛 = @view 𝛈[10:12]
#     𝐰 = @view 𝛈[end-(N-1):end]
#     r = range(13, step = 3, length = N+1)
#     𝐯 = reshape(view(𝛈,first(r):last(r)-1), (3,N))
#     ℋ = Array{Array{Float64,2},1}(undef,(N,))
#     for n = 1:N
#         ℋ[n] = 𝐰[n]*𝐀 + 𝐛*𝐯[:,n]'
#     end
#     ℋ
# end
