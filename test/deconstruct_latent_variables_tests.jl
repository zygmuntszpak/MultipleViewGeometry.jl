using MultipleViewGeometry,Test, LinearAlgebra
using MultipleViewGeometry.ModuleTypes
using StaticArrays, Calculus, GeometryTypes
using Random

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

𝐇₁ = construct(HomographyMatrix(),𝐊₁′,𝐑₁′,-𝐭₁′,𝐊₂′,𝐑₂′,-𝐭₂′,𝐧₁,d₁)
𝐇₂ = construct(HomographyMatrix(),𝐊₁′,𝐑₁′,-𝐭₁′,𝐊₂′,𝐑₂′,-𝐭₂′,𝐧₂,d₂)

𝛈 = deconstruct(LatentVariables(HomographyMatrices()), (𝐇₁,𝐇₂))
ℋ = compose(LatentVariables(HomographyMatrices()),𝛈)

𝐇₁′ = ℋ[1]
𝐇₂′ = ℋ[2]

@test 𝐇₁ ≈ 𝐇₁′
@test 𝐇₂ ≈ 𝐇₂′
