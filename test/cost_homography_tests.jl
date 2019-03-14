using MultipleViewGeometry, Test, Random
using MultipleViewGeometry.ModuleCostFunction
using MultipleViewGeometry.ModuleTypes
using LinearAlgebra
using StaticArrays
using GeometryTypes
using MultipleViewGeometry.ModuleSyntheticData

# Test for cost functions.

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
#
# 𝐇₁ = estimate(HomographyMatrix(), DirectLinearTransform(), (ℳ₁, ℳ₁ʹ))
# 𝐇₂ = estimate(HomographyMatrix(), DirectLinearTransform(), (ℳ₂, ℳ₂ʹ))

𝐇₁ = construct(HomographyMatrix(),𝐊₁′,𝐑₁′,-𝐭₁,𝐊₂′,𝐑₂′,-𝐭₂, 𝐧₁, d₁)
𝐇₂ = construct(HomographyMatrix(),𝐊₁′,𝐑₁′,-𝐭₁,𝐊₂′,𝐑₂′,-𝐭₂, 𝐧₂, d₂)
𝐡₁ = vec(𝐇₁)
𝐡₂ = vec(𝐇₂)

Λ₁ =  [SMatrix{3,3}(Matrix(Diagonal([1.0,1.0,0.0]))) for i = 1:length(ℳ₁)]
Λ₂ =  [SMatrix{3,3}(Matrix(Diagonal([1.0,1.0,0.0]))) for i = 1:length(ℳ₁)]
Jₐₘₗ =  cost(AML(),HomographyMatrix(), 𝐡₁, (Λ₁,Λ₂), (ℳ₁, ℳ₁ʹ))

@test isapprox(Jₐₘₗ, 0.0; atol = 1e-14)

# Verify that the vectorised fundamental matrix is in the null space of X
𝐗 = X(AML(),HomographyMatrix(), 𝐡₁, (Λ₁,Λ₂), (ℳ₁, ℳ₁ʹ))

# The true parameters should lie in the null space of the matrix X.
@test all(isapprox.(𝐗 * 𝐡₁, 0.0; atol = 1e-9))
