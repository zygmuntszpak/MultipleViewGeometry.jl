using MultipleViewGeometry, Test, Random
using MultipleViewGeometry.ModuleCostFunction
using MultipleViewGeometry.ModuleTypes
using MultipleViewGeometry.ModuleConstraints
using MultipleViewGeometry.ModuleConstruct
using LinearAlgebra
using StaticArrays
using GeometryTypes

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

# 3D points corresponding to the first and second planar surface
𝒴₁ = triangulate(DirectLinearTransform(),𝐏₁,𝐏₂,(ℳ₁,ℳ₁ʹ))
𝒴₂ = triangulate(DirectLinearTransform(),𝐏₁,𝐏₂,(ℳ₂,ℳ₂ʹ))

# Triangulating with the same projection matrices that were used to construct
# (ℳ,ℳʹ) should yield 3D points that lie on a plane.
N = length(𝒴₁)
for n = 1:N
    X = 𝒴₁[n]
    @test isapprox(dot(𝐧₁,X) - d₁, 0.0; atol = 1e-11)
end

N = length(𝒴₂)
for n = 1:N
    X = 𝒴₂[n]
    @test isapprox(dot(𝐧₂,X) - d₂, 0.0; atol = 1e-11)
end

# Triangulating with the same projection matrices that were used to construct
# (ℳ,ℳʹ) should yield the same 3D points as the original 𝒳.
N = length(𝒴₁)
for n = 1:N
    @test  isapprox(sum(abs.(𝒳₁[n]-𝒴₁[n])/3), 0.0; atol = 1e-12)
end

N = length(𝒴₂)
for n = 1:N
    @test  isapprox(sum(abs.(𝒳₂[n]-𝒴₂[n])/3), 0.0; atol = 1e-12)
end
