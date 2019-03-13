using MultipleViewGeometry, Test, LinearAlgebra
using MultipleViewGeometry.ModuleTypes
using StaticArrays, Calculus
# Tests for homography matrix estimation

# Normal to the planar surface.
𝐧 = [0.0, 0.0, -1.0]
# Distance of the plane from the origin.
d = 100
# Sample points on the planar surface.
𝒳 = [Point3D(x,y,d) for x = -100:5:100 for y = -100:5:100]
𝒳 = 𝒳[1:50:end]
# Intrinsic and extrinsic parameters of camera one.
𝐊₁ = Matrix{Float64}(I,3,3)
𝐑₁ = Matrix{Float64}(I,3,3)
𝐭₁ = [0.0, 0.0, 0.0]

# Intrinsic and extrinsic parameters of camera two.
𝐊₂ = Matrix{Float64}(I,3,3)
𝐑₂ = Matrix{Float64}(I,3,3)
𝐭₂ = [100.0, 2.0, -100.0]

# Camera projection matrices.
𝐏₁ = construct(ProjectionMatrix(),𝐊₁,𝐑₁,𝐭₁)
𝐏₂ = construct(ProjectionMatrix(),𝐊₂,𝐑₂,𝐭₂)

# Set of corresponding points.
ℳ = project(Pinhole(),𝐏₁,𝒳)
ℳʹ = project(Pinhole(),𝐏₂,𝒳)

# Estimate of the homography matrix and the true homography matrix.
𝐇 = estimate(HomographyMatrix(), DirectLinearTransform(), (ℳ, ℳʹ))
𝐇₀ = construct(HomographyMatrix(),𝐊₁,𝐑₁,𝐭₁,𝐊₂,𝐑₂,𝐭₂,𝐧,d)

𝐇₀ = 𝐇₀ / norm(𝐇₀)
𝐇₀ = 𝐇₀ / sign(𝐇₀[3,3])

𝐇 = 𝐇 / norm(𝐇)
𝐇 = 𝐇 / sign(𝐇[3,3])

for i = zip(ℳ,ℳʹ)
    m, mʹ =  i
    𝐦  = hom(m)
    𝐦ʹ = hom(mʹ)
    residual = vec2antisym(𝐦ʹ)*𝐇₀*𝐦
    @test isapprox(sum(residual), 0.0; atol = 1e-7)
end

for i = zip(ℳ,ℳʹ)
    m, mʹ =  i
    𝐦  = hom(m)
    𝐦ʹ = hom(mʹ)
    residual = vec2antisym(𝐦ʹ)*𝐇*𝐦
    @show residual
    @test isapprox(sum(residual), 0.0; atol = 1e-7)
end


#@test isapprox(sum(residual), 0.0; atol = 1e-7)
#dot(𝒳[10],𝐧) + d
