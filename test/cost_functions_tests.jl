using MultipleViewGeometry, Test
using MultipleViewGeometry.ModuleCostFunction
using MultipleViewGeometry.ModuleTypes
using BenchmarkTools, LinearAlgebra
using StaticArrays

# Test for cost functions.

# Test cost function on Fundamental matrix estimation.

𝒳 = [Point3DH(x,y,z,1.0)
                        for x=-1:0.5:10 for y=-1:0.5:10 for z=2:-0.1:1]

# Intrinsic and extrinsic parameters of camera one.
𝐊₁ = SMatrix{3,3}(Matrix{Float64}(I, 3, 3))
𝐑₁ = SMatrix{3,3}(Matrix{Float64}(I, 3, 3))
𝐭₁ =  @SVector [0.0, 0.0, -10]

# Intrinsic and extrinsic parameters of camera two.
𝐊₂ = SMatrix{3,3}(Matrix{Float64}(I, 3, 3))
𝐑₂ = SMatrix{3,3}(Matrix{Float64}(I, 3, 3)) #SMatrix{3,3,Float64,9}(rotxyz(pi/10,pi/10,pi/10))
𝐭₂ = @SVector [10.0, 10.0, -10.0]

# Camera projection matrices.
𝐏₁ = construct(ProjectionMatrix(),𝐊₁,𝐑₁,𝐭₁)
𝐏₂ = construct(ProjectionMatrix(),𝐊₂,𝐑₂,𝐭₂)

# Set of corresponding points.
ℳ = project(Pinhole(),𝐏₁,𝒳)
ℳʹ = project(Pinhole(),𝐏₂,𝒳)

𝐅 = construct(FundamentalMatrix(),𝐊₁,𝐑₁,𝐭₁,𝐊₂,𝐑₂,𝐭₂)

# Ensure the estimated and true matrix have the same scale and sign.
𝐅 = 𝐅 / norm(𝐅)
𝐅 = 𝐅 / sign(𝐅[1,3])
𝐟 = reshape(𝐅,9,1)

Λ₁ =  [SMatrix{3,3}(Matrix(Diagonal([1.0,1.0,0.0]))) for i = 1:length(ℳ)]
Λ₂ =  [SMatrix{3,3}(Matrix(Diagonal([1.0,1.0,0.0]))) for i = 1:length(ℳ)]
Jₐₘₗ =  cost(AML(),FundamentalMatrix(), SVector{9}(𝐟), (Λ₁,Λ₂), (ℳ, ℳʹ))

@test isapprox(Jₐₘₗ, 0.0; atol = 1e-14)

# Verify that the vectorised fundamental matrix is in the null space of X
𝐗 = X(AML(),FundamentalMatrix(), reshape(𝐅,9,1), (Λ₁,Λ₂), (ℳ, ℳʹ))

# The true parameters should lie in the null space of the matrix X.
@test all(isapprox.(𝐗 * 𝐟, 0.0; atol = 1e-10))

# Verify that the the vectorised fundamental matrix is in the null space of H.
# H represents the Hessian matrix of the AML cost function.
𝐇 = H(AML(),FundamentalMatrix(), reshape(𝐅,9,1), (Λ₁,Λ₂), (ℳ, ℳʹ))

# The true parameters should lie in the null space of the matrix H.
@test all(isapprox.(𝐇 * 𝐟, 0.0; atol = 1e-10))
