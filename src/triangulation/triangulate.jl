function triangulate(method::DirectLinearTransform, 𝐅::AbstractArray, 𝒟::Tuple{AbstractArray, Vararg{AbstractArray}})
    ℳ, ℳʹ =  𝒟
    𝐏₁, 𝐏₂ = construct(ProjectionMatrix(),𝐅)
    N = length(ℳ)
    𝒴 = Array{Point3D}(undef,N)
    for n = 1:N
        𝐦 = ℳ[n]
        𝐦ʹ = ℳʹ[n]
        eq1 = 𝐦[1] * 𝐏₁[3,:] - 𝐏₁[1,:]
        eq2 = 𝐦[2] * 𝐏₁[3,:] - 𝐏₁[2,:]
        eq3 = 𝐦ʹ[1] * 𝐏₂[3,:] - 𝐏₂[1,:]
        eq4 = 𝐦ʹ[2] * 𝐏₂[3,:] - 𝐏₂[2,:]
        𝐀 = SMatrix{4,4}(transpose(hcat(eq1,eq2,eq3,eq4)))
        # 𝐀 = [ (𝐦[1] * 𝐏₁[3,:] - 𝐏₁[1,:])' ;
        #       (𝐦[2] * 𝐏₁[3,:] - 𝐏₁[2,:])' ;
        #       (𝐦ʹ[1] * 𝐏₂[3,:] - 𝐏₂[1,:])' ;
        #       (𝐦ʹ[2] * 𝐏₂[3,:] - 𝐏₂[2,:])' ]
        # 𝐀₁ = vec2antisym(𝐦)*𝐏₁
        # 𝐀₂ = vec2antisym(𝐦ʹ)*𝐏₂
        # @show typeof(𝐀₁)
        # @show size(𝐀₁)
        # @show size(𝐀₂)
        # 𝐀 = vcat(𝐀₁,𝐀₂)
        # @show size(𝐀)
        # λ, f = smallest_eigenpair(Symmetric(Array(𝐀'*𝐀)))
        # @show λ
        # 𝒴[n] = Point3DH(𝑛(f))
        U,S,V = svd(𝐀)
        𝒴[n] = hom⁻¹(V[:,4])
    end
    𝒴
end

function triangulate(method::DirectLinearTransform, 𝐏₁::AbstractArray, 𝐏₂::AbstractArray, 𝒟::Tuple{AbstractArray, Vararg{AbstractArray}})
    ℳ, ℳʹ =  𝒟
    N = length(ℳ)
    𝒴 = Array{Point3D}(undef,N)
    for n = 1:N
        𝐦 = ℳ[n]
        𝐦ʹ = ℳʹ[n]
        eq1 = 𝐦[1] * 𝐏₁[3,:] - 𝐏₁[1,:]
        eq2 = 𝐦[2] * 𝐏₁[3,:] - 𝐏₁[2,:]
        eq3 = 𝐦ʹ[1] * 𝐏₂[3,:] - 𝐏₂[1,:]
        eq4 = 𝐦ʹ[2] * 𝐏₂[3,:] - 𝐏₂[2,:]
        𝐀 = SMatrix{4,4}(transpose(hcat(eq1,eq2,eq3,eq4)))
        # 𝐀 = [  (𝐦[1] * 𝐏₁[3,:] - 𝐏₁[1,:])' ;
        #         (𝐦[2] * 𝐏₁[3,:] - 𝐏₁[2,:])' ;
        #         (𝐦ʹ[1] * 𝐏₂[3,:] - 𝐏₂[1,:])' ;
        #         (𝐦ʹ[2] * 𝐏₂[3,:] - 𝐏₂[2,:])' ]

        # 𝐀₁ = vec2antisym(𝐦)*𝐏₁
        # 𝐀₂ = vec2antisym(𝐦ʹ)*𝐏₂
        # @show typeof(𝐀₁)
        # @show size(𝐀₁)
        # @show size(𝐀₂)
        # 𝐀 = vcat(𝐀₁,𝐀₂)
        # @show size(𝐀)
        # λ, f = smallest_eigenpair(Symmetric(Array(𝐀'*𝐀)))
        # @show λ
        # 𝒴[n] = Point3DH(𝑛(f))
        U,S,V = svd(𝐀)
        𝒴[n] = hom⁻¹(V[:,4])
    end
    𝒴
end
