function estimate(entity::FundamentalMatrix, method::DirectLinearTransform, 𝒟::Tuple{AbstractArray, Vararg{AbstractArray}})
    ℳ, ℳʹ =  𝒟
    N = length(ℳ)
    if (N != length(ℳʹ))
          throw(ArgumentError("There should be an equal number of points for each view."))
    end
    𝒪, 𝐓 = transform(HomogeneousCoordinates(),CanonicalToHartley(),ℳ)
    𝒪ʹ, 𝐓ʹ  = transform(HomogeneousCoordinates(),CanonicalToHartley(),ℳʹ)
    𝐀 = moments(FundamentalMatrix(), (𝒪, 𝒪ʹ))
    λ, f = smallest_eigenpair(Symmetric(𝐀))
    𝐅 = reshape(f,(3,3))
    𝐅 = enforce_ranktwo!(Array(𝐅))
    𝐅 = 𝐅 / norm(𝐅)
    # Transform estimate back to the original (unnormalised) coordinate system.
    𝐓ʹ'*𝐅*𝐓
end

# TODO fix numerical instability
# function estimate(entity::FundamentalMatrix, method::Taubin, matches...)
#     ℳ, ℳʹ = matches
#     N = length(ℳ)
#     if (N != length(ℳʹ))
#           throw(ArgumentError("There should be an equal number of points for each view."))
#     end
#     (ℳ,𝐓) = hartley_normalization(ℳ)
#     (ℳʹ,𝐓ʹ) = hartley_normalization(ℳʹ)
#     𝐀::Matrix{Float64} = moments(FundamentalMatrix(), ℳ, ℳʹ)
#     𝐁::Matrix{Float64} = mean_covariance(FundamentalMatrix(), ℳ, ℳʹ)
#     (λ::Float64, f::Vector{Float64}) = smallest_eigenpair(𝐀,𝐁)
#     𝐅::Matrix{Float64} = reshape(f,(3,3))
#     enforce_ranktwo!(𝐅)
#     # Transform estimate back to the original (unnormalised) coordinate system.
#     𝐅 = 𝐓ʹ'*𝐅*𝐓
# end

function estimate(entity::FundamentalMatrix, method::FundamentalNumericalScheme,  𝒞::Tuple{AbstractArray, Vararg{AbstractArray}}, 𝒟::Tuple{AbstractArray, Vararg{AbstractArray}})
    ℳ, ℳʹ = 𝒟
    Λ₀, Λ₀ʹ = 𝒞
    N = length(ℳ)
    if (N != length(ℳʹ))
          throw(ArgumentError("There should be an equal number of points for each view."))
    end
    if (N != length(Λ₀) || N != length(Λ₀ʹ) )
          throw(ArgumentError("There should be a covariance matrix for each point correspondence."))
    end
    # Map corresponding points to the normalized coordinate system.
    𝒪, 𝐓 = transform(HomogeneousCoordinates(),CanonicalToHartley(),ℳ)
    𝒪ʹ, 𝐓ʹ = transform(HomogeneousCoordinates(),CanonicalToHartley(),ℳʹ)
    # Map seed to the normalized coordinate system.
    𝛉 = (inv(𝐓') ⊗ inv(𝐓ʹ')) * method.𝛉₀
    # Map covariance matrices to the normalized coordinate system.
    Λ₁ = transform(CovarianceMatrices(), CanonicalToHartley(), Λ₀ , 𝐓)
    Λ₁ʹ = transform(CovarianceMatrices(), CanonicalToHartley(), Λ₀ʹ , 𝐓ʹ)
    for i = 1:method.max_iter
        𝐗 = X(AML(),FundamentalMatrix(), 𝛉, (Λ₁,Λ₁ʹ), (𝒪, 𝒪ʹ))
        λ, 𝛉⁺ = smallest_eigenpair(Symmetric(𝐗/N))
        𝛉 = reshape(𝛉⁺,length(𝛉⁺),1)
    end
    𝐅 = reshape(𝛉,(3,3))
    𝐅 = enforce_ranktwo!(Array(𝐅))
    # Transform estimate back to the original (unnormalised) coordinate system.
    𝐅 = 𝐓ʹ'*𝐅*𝐓
end

function estimate(entity::FundamentalMatrix, method::BundleAdjustment,  𝒟::Tuple{AbstractArray, Vararg{AbstractArray}})
    ℳ, ℳʹ = 𝒟
    N = length(ℳ)
    if (N != length(ℳʹ))
          throw(ArgumentError("There should be an equal number of points for each view."))
    end
    𝐅 = reshape(method.𝛉₀,(3,3))
    𝒳 = triangulate(DirectLinearTransform(),𝐅,(ℳ,ℳʹ))

    𝐏₁, 𝐏₂ = construct(ProjectionMatrix(),𝐅)

    # Construct a length-(12+3*N) vector consisting of the projection matrix associated
    # with the second view (the first twelve dimensions), as well as N three-dimensional points
    # (the remaining dimensions).
    𝛉 = pack(FundamentalMatrix(), 𝐏₂, 𝒳)

    index₁ = SVector(1,2)
    index₂ = SVector(3,4)
    pts = Matrix{Float64}(undef,4,N)
    for n = 1:N
        pts[index₁,n] = ℳ[n][index₁]
        pts[index₂,n] = ℳʹ[n][index₁]
    end


    #fit = curve_fit(model_fundamental, jacobian_model, 𝐏₁, reinterpret(Float64,pts,(4*N,1)), 𝛉; show_trace = false)
    #fit = curve_fit(model_fundamental, jacobian_model, 𝐏₁, temp, 𝛉; show_trace = false)
    fit = curve_fit(model_fundamental, jacobian_model,  𝐏₁, reshape(reinterpret(Float64,vec(pts)),(4*N,)) , 𝛉; show_trace = false)
    #reshape(reinterpret(T, vec(a)), dims)
    #reinterpret(::Type{T}, a::Array{S}, dims::NTuple{N, Int}) where {T, S, N}
    #@show typeof(reshape(reinterpret(Float64,vec(pts)),(4*N,)))
    #@show typeof(reinterpret(Float64,pts,(4*N,)))
    #fit = curve_fit(model_fundamental, jacobian_model, 𝐏₁, reshape(reinterpret(Float64,pts),(4*N,)), 𝛉; show_trace = false)
    𝐏₂ = reshape(fit.param[1:12],(3,4))
    𝐅 = construct(FundamentalMatrix(), 𝐏₁, 𝐏₂)
    𝐅, fit
end

function model_fundamental(𝐏₁,𝛉)
    # Twelve parameters for the projection matrix, and 3 parameters per 3D point.
    N = Int((length(𝛉) - 12) / 3)
    index₁ = SVector(1,2)
    index₂ = SVector(3,4)
    reprojections = Matrix{Float64}(undef,4,N)
    𝛉v = @view 𝛉[1:12]
    𝐏₂ = SMatrix{3,4,Float64,12}(reshape(𝛉v,(3,4)))
    i = 13
    for n = 1:N
        # Extract 3D point and convert to homogeneous coordinates
        v = @view 𝛉[i:i+2]
        M = hom(SVector{3,Float64}(𝛉[i:i+2]))
        reprojections[index₁,n] = hom⁻¹(𝐏₁ * M)
        reprojections[index₂,n] = hom⁻¹(𝐏₂ * M)
        i = i + 3
    end
    reshape(reinterpret(Float64,vec(reprojections)),(4*N,))
end

function jacobian_model(𝐏₁,𝛉)
    # Twelve parameters for the projection matrix, and 3 parameters per 3D point.
    N = Int((length(𝛉) - 12) / 3)
    index₁ = SVector(1,2)
    index₂ = SVector(3,4)
    reprojections = Matrix{Float64}(undef,4,N)
    𝛉v = @view 𝛉[1:12]
    𝐏₂ = SMatrix{3,4,Float64,12}(reshape(𝛉v,(3,4)))
    𝐉 = zeros(4*N,12+3*N)
    # Create a view of the jacobian matrix 𝐉 and reshape it so that
    # it will be more convenient to index into the appropriate entires
    # whilst looping over all of the data points.
    𝐉v = reshape(reinterpret(Float64,𝐉), 4, N, 12+3*N)
    𝐀 = SMatrix{2,3,Float64,6}(1,0,0,1,0,0)
    𝐈₃ = SMatrix{3,3}(1.0I)
    i = 13
    for n = 1:N
        # Extract 3D point and convert to homogeneous coordinates
        𝐌 = hom(SVector{3,Float64}(𝛉[i:i+2]))

        # Derivative of residual in first and second image w.r.t 3D point.
        ∂𝐫₁_d𝐌 = -𝐀 * ∂hom⁻¹(𝐏₁ * 𝐌) * 𝐏₁
        ∂𝐫₂_d𝐌 = -𝐀 * ∂hom⁻¹(𝐏₂ * 𝐌) * 𝐏₂

        # Derivative of residual in second image w.r.t projection martix
        # ∂𝐫₁_d𝐏₁ is the zero vector.
        ∂𝐫₂_d𝐏₂ = 𝐀 * ∂hom⁻¹(𝐏₂ * 𝐌) * (𝐌' ⊗ 𝐈₃)

        𝐉v[index₂,n,1:12] = ∂𝐫₂_d𝐏₂
        𝐉v[index₁,n,i:i+2] = ∂𝐫₁_d𝐌[:,1:3]
        𝐉v[index₂,n,i:i+2] = ∂𝐫₂_d𝐌[:,1:3]
        i = i + 3
    end
    𝐉
end

#𝛉 = reshape(𝛉₀,length(𝛉₀),1)
# function z(entity::FundamentalMatrix, 𝐦::Vector{Float64}, 𝐦ʹ::Vector{Float64})
# 𝐮 = 𝐦 ⊗ 𝐦ʹ
# 𝐮[1:end-1]
# end

function mean_covariance(entity::ProjectiveEntity, matches...)
    ℳ, ℳʹ = matches
    N = length(ℳ)
    if (N != length(ℳʹ))
          throw(ArgumentError("There should be an equal number of points for each view."))
    end
    𝐁 = fill(0.0,(9,9))
    𝚲 = eye(4)
    for correspondence in zip(ℳ, ℳʹ)
        m , mʹ = correspondence
        𝐦  = 𝑛(collect(Float64,m.coords))
        𝐦ʹ = 𝑛(collect(Float64,mʹ.coords))
        ∂ₓ𝐮 = ∂ₓu(entity, 𝐦 , 𝐦ʹ)
        𝐁 = 𝐁 + ∂ₓ𝐮 * 𝚲 * ∂ₓ𝐮'
    end
    𝐁/N
end


function enforce_ranktwo!(𝐅::AbstractArray)
    # Enforce the rank-2 constraint.
    U,S,V = svd(𝐅)
    S[end] = 0.0
    𝐅 = U*Matrix(Diagonal(S))*V'
end

# Construct a parameter vector consisting of a projection matrix and 3D points
function pack(entity::FundamentalMatrix, 𝐏₂::AbstractArray, 𝒳::AbstractArray, )
    N = length(𝒳)
    𝛉 = Vector{Float64}(undef,12+N*3)
    𝛉[1:12] = Array(𝐏₂[:])
    i = 13
    for n = 1:N
        𝛉[i:i+2] = 𝒳[n][1:3]
        i = i + 3
    end
    𝛉
end
