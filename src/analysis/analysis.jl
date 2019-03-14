
function assess(c::ReprojectionError, entity::HomographyMatrix, ℋ::AbstractArray, 𝓓::Tuple{Vector{Vector{T₁}} where T₁ <: AbstractArray, Vector{Vector{T₂}} where T₂ <: AbstractArray})
    𝓜, 𝓜ʹ =  𝓓
    results = Vector{Float64}(undef,length(𝓓))
    for k = 1:length(𝓓)
        results[k] = assess(c, entity, ℋ[k], (𝓜[k], 𝓜ʹ[k]))
    end
    results
end

function assess(c::ReprojectionError, entity::HomographyMatrix, 𝐇::AbstractArray, 𝒟::Tuple{AbstractArray, AbstractArray})
    ℳ, ℳʹ = 𝒟
    N = length(𝒟[1])
    if (N != length(ℳʹ))
          throw(ArgumentError("There should be an equal number of points for each view."))
    end
    # Construct a length-(2*N) vector consisting of N two-dimensional points in the
    # first view.
    𝛉 = Vector{Float64}(undef, N*2)
    i = 1
    for n = 1:N
        𝛉[i:i+1] = @view ℳ[n][1:2]
        i = i + 2
    end

    index₁ = SVector(1,2)
    index₂ = SVector(3,4)
    pts = Matrix{Float64}(undef,4,N)
    for n = 1:N
        pts[index₁,n] = ℳ[n][index₁]
        pts[index₂,n] = ℳʹ[n][index₁]
    end
    #fit = curve_fit(model_homography,  𝐇, reshape(reinterpret(Float64,vec(pts)),(4*N,)) , 𝛉; show_trace = false, maxIter = 2)
    fit = curve_fit(model_homography, jacobian_model_homography, 𝐇, reshape(reinterpret(Float64,vec(pts)),(4*N,)) , 𝛉;  show_trace = false)
    # TODO Investigate NaN for initial values of Jacobian
    #fit = curve_fit(model_homography!, jacobian_model_homography!, 𝐇, reshape(reinterpret(Float64,vec(pts)),(4*N,)) , 𝛉;  inplace = true, show_trace = false, maxIter = 5)
    sum(fit.resid.^2)
end

function model_homography(𝐇,𝛉)
    # 2 parameters per 2D point.
    N = Int(length(𝛉)/ 2)
    index₁ = SVector(1,2)
    index₂ = SVector(3,4)
    reprojections = Matrix{Float64}(undef,4,N)
    i = 1
    for n = 1:N
        # Extract 2D point and convert to homogeneous coordinates
        𝐦 = hom(SVector{2,Float64}(𝛉[i],𝛉[i+1]))
        reprojections[index₁,n] = hom⁻¹(𝐦)
        reprojections[index₂,n] = hom⁻¹(𝐇 * 𝐦)
        i = i + 2
    end
    reshape(reinterpret(Float64,vec(reprojections)),(4*N,))
end

function model_homography!(reprojections::Array{Float64,1},𝐇,𝛉)
    # 2 parameters per 2D point.
    N = Int(length(𝛉)/ 2)
    index₁ = SVector(1,2)
    index₂ = SVector(3,4)
    reprojections_view = reshape(reinterpret(Float64,reprojections),(4,N))
    i = 1
    for n = 1:N
        # Extract 2D point and convert to homogeneous coordinates
        𝐦 = hom(SVector{2,Float64}(𝛉[i],𝛉[i+1]))
        reprojections_view[index₁,n] = hom⁻¹(𝐦)
        reprojections_view[index₂,n] = hom⁻¹(𝐇 * 𝐦)
        i = i + 2
    end
    reprojections
end

function jacobian_model_homography(𝐇,𝛉)
    # 2 parameters per 2D point.
    N = Int(length(𝛉) / 2)
    index₁ = SVector(1,2)
    index₂ = SVector(3,4)
    𝐉 = zeros(4*N,2*N)
    # Create a view of the jacobian matrix 𝐉 and reshape it so that
    # it will be more convenient to index into the appropriate entries
    # whilst looping over all of the data points.
    𝐉v = reshape(reinterpret(Float64,𝐉), 4, N, 2*N)
    𝐀 = SMatrix{2,3,Float64,6}(1,0,0,1,0,0)
    𝐈₃ = SMatrix{3,3}(1.0I)
    i = 1
    for n = 1:N
        # Extract 3D point and convert to homogeneous coordinates.
        𝐦 = hom(SVector{2,Float64}(𝛉[i], 𝛉[i+1]))

        # Derivative of residual in first and second image w.r.t 2D point in the
        # first image.
        ∂𝐫₁_d𝐦 = 𝐀 * 𝐈₃
        ∂𝐫₂_d𝐦 = 𝐀 * ∂hom⁻¹(𝐇 * 𝐦) * 𝐇
    @.  𝐉v[index₁,n,i:i+1] = ∂𝐫₁_d𝐦[:,index₁]
    @.  𝐉v[index₂,n,i:i+1] = ∂𝐫₂_d𝐦[:,index₁]
        i = i + 2
    end
    𝐉
end

function jacobian_model_homography!(𝐉::Array{Float64,2}, 𝐇,𝛉)
    Base.display(𝐉)
    pause
    # 2 parameters per 2D point.
    N = Int(length(𝛉) / 2)
    index₁ = SVector(1,2)
    index₂ = SVector(3,4)
    # Create a view of the jacobian matrix 𝐉 and reshape it so that
    # it will be more convenient to index into the appropriate entries
    # whilst looping over all of the data points.
    𝐉v = reshape(reinterpret(Float64,𝐉), 4, N, 2*N)
    𝐀 = SMatrix{2,3,Float64,6}(1,0,0,1,0,0)
    𝐈₃ = SMatrix{3,3}(1.0I)
    i = 1
    for n = 1:N
        # Extract 3D point and convert to homogeneous coordinates.
        𝐦 = hom(SVector{2,Float64}(𝛉[i], 𝛉[i+1]))

        # Derivative of residual in first and second image w.r.t 2D point in the
        # first image.
        ∂𝐫₁_d𝐦 = 𝐀 * 𝐈₃
        ∂𝐫₂_d𝐦 = 𝐀 * ∂hom⁻¹(𝐇 * 𝐦) * 𝐇
    @.  𝐉v[index₁,n,i:i+1] = ∂𝐫₁_d𝐦[:,index₁]
    @.  𝐉v[index₂,n,i:i+1] = ∂𝐫₂_d𝐦[:,index₁]
        i = i + 2
    end
    𝐉
end
