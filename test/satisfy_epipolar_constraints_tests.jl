using MultipleViewGeometry, Test
using MultipleViewGeometry.ModuleCostFunction
using MultipleViewGeometry.ModuleTypes
using MultipleViewGeometry.ModuleConstraints
using StaticArrays, Random, LinearAlgebra

# Fix random seed.
Random.seed!(1234)

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

𝐅 = construct(FundamentalMatrix(),𝐊₁′,𝐑₁′,-𝐭₁′,𝐊₂′,𝐑₂′,-𝐭₂′)

# Verify that the algorithm returns the correct answer when the
# constraint is already satisfied.
𝒪,𝒪ʹ = satisfy(FundamentalMatrix(), EpipolarConstraint(), 𝐅, (ℳ, ℳʹ))

# Verify that the original corresponding points satisfy the epipolar constraint.
N = length(ℳ)
for n = 1:N
    𝐦 = hom(ℳ[n])
    𝐦ʹ = hom(ℳʹ[n])
    @test  isapprox(𝐦'*𝐅*𝐦ʹ, 0.0; atol = 1e-12)
end

# Verify that the 'corrected' points satisfy the epipolar constraint.
N = length(ℳ)
for n = 1:N
    𝐦 = hom(𝒪[n])
    𝐦ʹ = hom(𝒪ʹ[n])
    @test  isapprox(𝐦'*𝐅*𝐦ʹ, 0.0; atol = 1e-12)
end

# Perturb the original corresponding points slightly so that they no-longer
# satisfy the epipolar constraint.
N = length(ℳ)
σ = 1e-7
for n = 1:N
    ℳ[n] = ℳ[n] + SVector{2}(σ * rand(2,1))
    ℳʹ[n] = ℳʹ[n] + SVector{2}(σ * rand(2,1))
    𝐦 = hom(ℳ[n])
    𝐦ʹ = hom(ℳʹ[n])
    @test abs(𝐦'*𝐅*𝐦ʹ) > 1e-12
end


# Verify that the algorithm returns the correct answer when applied
# to sets of correspondences that do not satisfy the epipolar constraint.
𝒪 ,𝒪ʹ = satisfy(FundamentalMatrix(), EpipolarConstraint(), 𝐅, (ℳ, ℳʹ))

# Verify that the 'corrected' points satisfy the epipolar constraint.
N = length(ℳ)
for n = 1:N
    𝐦 = hom(𝒪[n])
    𝐦ʹ = hom(𝒪ʹ[n])
    @test  isapprox(𝐦'*𝐅*𝐦ʹ, 0.0; atol = 1e-12)
end
