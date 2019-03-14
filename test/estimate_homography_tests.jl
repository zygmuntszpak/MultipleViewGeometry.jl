using MultipleViewGeometry, Test, LinearAlgebra
using MultipleViewGeometry.ModuleTypes
using StaticArrays, Calculus, GeometryTypes
using MultipleViewGeometry.ModuleSyntheticData

# Fix random seed.
Random.seed!(1234)
# Construct two camera matrices and parametrise two planar surfaces.
f = 50
image_width = 640 / 10
image_height = 480 / 10
𝐊₁ = @SMatrix [f 0 0 ;
               0 f 0 ;
               0 0 1 ]
𝐑₁ = SMatrix{3,3,Float64,9}(rotxyz(0, 25*(pi/180), 0))
𝐭₁ = [-30.0, 0.0, -5.0]

𝐊₂ = @SMatrix [f 0 0 ;
               0 f 0 ;
               0 0 1 ]

𝐑₂ = SMatrix{3,3,Float64,9}(rotxyz(0, -25*(pi/180), 0))
𝐭₂ = [30.0, 0.0, 5.0]

# Normals and distance from origin
𝐧₁ = [0.0, 0.0, 1.0]
d₁ = 55.0

# Normals and distance from origin
𝐧₁ = [0.0, 0.0, 1.0]
d₁ = 55.0

𝐧₂ = [0.5, -0.2, 2.0]
d₂ = 145.0

𝒳₁ = generate_planar_points(𝐧₁,d₁, 20, 50)
𝒳₂ = generate_planar_points(𝐧₂,d₂, 20, 50)


world_basis = (Vec(1.0, 0.0, 0.0), Vec(0.0, 1.0, 0.0), Vec(0.0, 0.0, 1.0))
camera_basis = (Point(0.0, 0.0, 0.0), Vec(-1.0, 0.0, 0.0), Vec(0.0, -1.0, 0.0), Vec(0.0, 0.0, 1.0))
picture_basis = (Point(0.0, 0.0), Vec(-1.0, 0.0), Vec(0.0, -1.0))

camera₁ = Pinhole(image_width, image_height, f, camera_basis..., picture_basis...)
camera₂ = Pinhole(image_width, image_height, f, camera_basis..., picture_basis...)
relocate!(camera₁, 𝐑₁, 𝐭₁)
relocate!(camera₂, 𝐑₂, 𝐭₂)

𝐑₁′, 𝐭₁′ = ascertain_pose(camera₁, world_basis... )
𝐊₁′ = obtain_intrinsics(camera₁, CartesianSystem())
𝐑₂′, 𝐭₂′ = ascertain_pose(camera₂, world_basis... )
𝐊₂′ = obtain_intrinsics(camera₂, CartesianSystem())

𝐏₁ = construct(ProjectionMatrix(),𝐊₁′,𝐑₁′,𝐭₁′)
𝐏₂ = construct(ProjectionMatrix(),𝐊₂′,𝐑₂′,𝐭₂′)

# Set of corresponding points.
ℳ₁ = project(camera₁,𝐏₁,𝒳₁)
ℳ₁ʹ= project(camera₂,𝐏₂,𝒳₁)
ℳ₂ = project(camera₁,𝐏₁,𝒳₂)
ℳ₂ʹ= project(camera₂,𝐏₂,𝒳₂)

𝐇₁ = estimate(HomographyMatrix(), DirectLinearTransform(), (ℳ₁, ℳ₁ʹ))
𝐇₂ = estimate(HomographyMatrix(), DirectLinearTransform(), (ℳ₂, ℳ₂ʹ))


for i = zip(ℳ₁,ℳ₁ʹ)
    m, mʹ =  i
    𝐦  = hom(m)
    𝐦ʹ = hom(mʹ)
    residual = vec2antisym(𝐦ʹ)*𝐇₁*𝐦
    @test isapprox(sum(residual), 0.0; atol = 1e-7)
end

for i = zip(ℳ₂,ℳ₂ʹ)
    m, mʹ =  i
    𝐦  = hom(m)
    𝐦ʹ = hom(mʹ)
    residual = vec2antisym(𝐦ʹ)*𝐇₂*𝐦
    @test isapprox(sum(residual), 0.0; atol = 1e-7)
end

# Tests for homography matrix estimation
#
# # Normal to the planar surface.
# 𝐧 = [0.0, 0.0, -1.0]
# # Distance of the plane from the origin.
# d = 100
# # Sample points on the planar surface.
# 𝒳 = [Point3D(x,y,d) for x = -100:5:100 for y = -100:5:100]
# 𝒳 = 𝒳[1:50:end]
# # Intrinsic and extrinsic parameters of camera one.
# 𝐊₁ = Matrix{Float64}(I,3,3)
# 𝐑₁ = Matrix{Float64}(I,3,3)
# 𝐭₁ = [0.0, 0.0, 0.0]
#
# # Intrinsic and extrinsic parameters of camera two.
# 𝐊₂ = Matrix{Float64}(I,3,3)
# 𝐑₂ = Matrix{Float64}(I,3,3)
# 𝐭₂ = [100.0, 2.0, -100.0]
#
# # Camera projection matrices.
# 𝐏₁ = construct(ProjectionMatrix(),𝐊₁,𝐑₁,𝐭₁)
# 𝐏₂ = construct(ProjectionMatrix(),𝐊₂,𝐑₂,𝐭₂)
#
# # Set of corresponding points.
# ℳ = project(Pinhole(),𝐏₁,𝒳)
# ℳʹ = project(Pinhole(),𝐏₂,𝒳)
#
# # Estimate of the homography matrix and the true homography matrix.
# 𝐇 = estimate(HomographyMatrix(), DirectLinearTransform(), (ℳ, ℳʹ))
# 𝐇₀ = construct(HomographyMatrix(),𝐊₁,𝐑₁,𝐭₁,𝐊₂,𝐑₂,𝐭₂,𝐧,d)
#
# 𝐇₀ = 𝐇₀ / norm(𝐇₀)
# 𝐇₀ = 𝐇₀ / sign(𝐇₀[3,3])
#
# 𝐇 = 𝐇 / norm(𝐇)
# 𝐇 = 𝐇 / sign(𝐇[3,3])
#
# for i = zip(ℳ,ℳʹ)
#     m, mʹ =  i
#     𝐦  = hom(m)
#     𝐦ʹ = hom(mʹ)
#     residual = vec2antisym(𝐦ʹ)*𝐇₀*𝐦
#     @test isapprox(sum(residual), 0.0; atol = 1e-7)
# end
#
# for i = zip(ℳ,ℳʹ)
#     m, mʹ =  i
#     𝐦  = hom(m)
#     𝐦ʹ = hom(mʹ)
#     residual = vec2antisym(𝐦ʹ)*𝐇*𝐦
#     @show residual
#     @test isapprox(sum(residual), 0.0; atol = 1e-7)
# end
#

#@test isapprox(sum(residual), 0.0; atol = 1e-7)
#dot(𝒳[10],𝐧) + d
