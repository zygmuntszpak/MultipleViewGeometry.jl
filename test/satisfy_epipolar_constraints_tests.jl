using MultipleViewGeometry, Test
using MultipleViewGeometry.ModuleCostFunction
using MultipleViewGeometry.ModuleTypes
using MultipleViewGeometry.ModuleConstraints
using BenchmarkTools
using StaticArrays, Random, LinearAlgebra

# Fix random seed.
Random.seed!(1234)

# Test for cost functions.

# Test cost function on Fundamental matrix estimation.

𝒳 = [Point3DH(x,y,z,1.0)
                        for x=-1:0.5:10 for y=-1:0.5:10 for z=2:-0.1:1]

# Intrinsic and extrinsic parameters of camera one.
𝐊₁ = SMatrix{3,3}(1.0I)
𝐑₁ = SMatrix{3,3}(1.0I)
𝐭₁ =  @SVector [0.0, 0.0, -10]

# Intrinsic and extrinsic parameters of camera two.
𝐊₂ = SMatrix{3,3}(1.0I)
𝐑₂ = SMatrix{3,3}(1.0I) #SMatrix{3,3,Float64,9}(rotxyz(pi/10,pi/10,pi/10))
𝐭₂ = @SVector [10.0, 10.0, -10.0]

# Camera projection matrices.
𝐏₁ = construct(ProjectionMatrix(),𝐊₁,𝐑₁,𝐭₁)
𝐏₂ = construct(ProjectionMatrix(),𝐊₂,𝐑₂,𝐭₂)

# Set of corresponding points.
ℳ = project(Pinhole(),𝐏₁,𝒳)
ℳʹ = project(Pinhole(),𝐏₂,𝒳)

𝐅 = construct(FundamentalMatrix(),𝐊₁,𝐑₁,𝐭₁,𝐊₂,𝐑₂,𝐭₂)

# Verify that the algorithm returns the correct answer when the
# constraint is already satisfied.
𝒪 ,𝒪ʹ = satisfy(FundamentalMatrix(), EpipolarConstraint(), 𝐅, (ℳ, ℳʹ))

# Verify that the original corresponding points satisfy the epipolar constraint.
N = length(ℳ)
for n = 1:N
    𝐦 = ℳ[n]
    𝐦ʹ = ℳʹ[n]
    @test  isapprox(𝐦'*𝐅*𝐦ʹ, 0.0; atol = 1e-14)
end

# Verify that the 'corrected' points satisfy the epipolar constraint.
N = length(ℳ)
for n = 1:N
    𝐦 = 𝒪[n]
    𝐦ʹ = 𝒪ʹ[n]
    @test  isapprox(𝐦'*𝐅*𝐦ʹ, 0.0; atol = 1e-14)
end

# Perturb the original corresponding points slightly so that they no-longer
# satisfy the epipolar constraint.
N = length(ℳ)
σ = 1e-7
for n = 1:N
    ℳ[n] = ℳ[n] + SVector{3}(σ * vcat(rand(2,1),0))
    ℳʹ[n] = ℳʹ[n] + SVector{3}(σ * vcat(rand(2,1),0))
    𝐦 = ℳ[n]
    𝐦ʹ = ℳʹ[n]
    @test abs(𝐦'*𝐅*𝐦ʹ) > 1e-12
end


# Verify that the algorithm returns the correct answer when applied
# to sets of correspondences that do not satisfy the epipolar constraint.
𝒪 ,𝒪ʹ = satisfy(FundamentalMatrix(), EpipolarConstraint(), 𝐅, (ℳ, ℳʹ))

# Verify that the 'corrected' points satisfy the epipolar constraint.
N = length(ℳ)
for n = 1:N
    𝐦 = 𝒪[n]
    𝐦ʹ = 𝒪ʹ[n]
    @test  isapprox(𝐦'*𝐅*𝐦ʹ, 0.0; atol = 1e-14)
end

# σ = 1e-7
# SVector{3}(σ * vcat(rand(2,1),0))
#
# ℳ[1] - 𝒪[1]
