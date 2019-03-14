using MultipleViewGeometry, Test
using MultipleViewGeometry.ModuleCostFunction
using MultipleViewGeometry.ModuleTypes
using LinearAlgebra
using StaticArrays
using GeometryTypes

# Test for cost functions.

# Test cost function on Fundamental matrix estimation.

𝒳 = [Point3D(x,y,rand(50:100)) for x = -100:5:100 for y = -100:5:100]
𝒳 = 𝒳[1:50:end]


# Specify the coordinate systems of the world, the camera frame and the picture
# plane.
world_basis = (Vec(1.0, 0.0, 0.0), Vec(0.0, 1.0, 0.0), Vec(0.0, 0.0, 1.0))
camera_basis = (Point(0.0, 0.0, 0.0), Vec(-1.0, 0.0, 0.0), Vec(0.0, -1.0, 0.0), Vec(0.0, 0.0, 1.0))
picture_basis = (Point(0.0, 0.0), Vec(-1.0, 0.0), Vec(0.0, -1.0))

# The focal length for both cameras is one.
f = 1
image_width = 640 / 10
image_height = 480 / 10

camera₁ = Pinhole(image_width, image_height, f, camera_basis..., picture_basis...)
camera₂ = Pinhole(image_width, image_height, f, camera_basis..., picture_basis...)

# Rotate and translate camera one.
𝐑₁ = Matrix{Float64}(I,3,3)
𝐭₁ = [-50.0, -2.0, 0.0]
relocate!(camera₁, 𝐑₁, 𝐭₁)

# Rotate and translate camera two.
𝐑₂ = Matrix{Float64}(I,3,3)
𝐭₂ = [50.0, 2.0, 0.0]
relocate!(camera₂, 𝐑₂, 𝐭₂)


𝐑₁′, 𝐭₁′ = ascertain_pose(camera₁, world_basis... )
𝐊₁′ = obtain_intrinsics(camera₁, CartesianSystem())
𝐑₂′, 𝐭₂′ = ascertain_pose(camera₂, world_basis... )
𝐊₂′ = obtain_intrinsics(camera₂, CartesianSystem())

# Camera projection matrices.
𝐏₁ = construct(ProjectionMatrix(),𝐊₁′,𝐑₁′,𝐭₁′)
𝐏₂ = construct(ProjectionMatrix(),𝐊₂′,𝐑₂′,𝐭₂′)

# Set of corresponding points.
ℳ = project(camera₁,𝐏₁,𝒳)
ℳʹ = project(camera₂,𝐏₂,𝒳)

# Construct fundamental matrix from projection matrices.
𝐅 = construct(FundamentalMatrix(), 𝐏₁, 𝐏₂)

# Ensure the estimated and true matrix have the same scale and sign.
𝐅 = 𝐅 / norm(𝐅)
𝐅 = 𝐅 / sign(𝐅[3,1])
𝐟 = vec(𝐅)

Λ₁ =  [SMatrix{3,3}(Matrix(Diagonal([1.0,1.0,0.0]))) for i = 1:length(ℳ)]
Λ₂ =  [SMatrix{3,3}(Matrix(Diagonal([1.0,1.0,0.0]))) for i = 1:length(ℳ)]
Jₐₘₗ =  cost(AML(),FundamentalMatrix(), SVector{9}(𝐟), (Λ₁,Λ₂), (ℳ, ℳʹ))

@test isapprox(Jₐₘₗ, 0.0; atol = 1e-14)

# Verify that the vectorised fundamental matrix is in the null space of X
𝐗 = X(AML(),FundamentalMatrix(), vec(𝐅), (Λ₁,Λ₂), (ℳ, ℳʹ))

# The true parameters should lie in the null space of the matrix X.
@test all(isapprox.(𝐗 * 𝐟, 0.0; atol = 1e-9))

# Verify that the the vectorised fundamental matrix is in the null space of H.
# H represents the Hessian matrix of the AML cost function.
𝐇 = H(AML(),FundamentalMatrix(), vec(𝐅), (Λ₁,Λ₂), (ℳ, ℳʹ))

# The true parameters should lie in the null space of the matrix H.
@test all(isapprox.(𝐇 * 𝐟, 0.0; atol = 1e-9))
