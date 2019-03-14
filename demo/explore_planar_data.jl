using Makie
using MultipleViewGeometry, Test, Random
using MultipleViewGeometry.ModuleCostFunction
using MultipleViewGeometry.ModuleTypes
using MultipleViewGeometry.ModuleConstraints
using MultipleViewGeometry.ModuleConstruct
using MultipleViewGeometry.ModuleDraw
using MultipleViewGeometry.ModuleMove
using MultipleViewGeometry.ModuleSyntheticData
using LinearAlgebra
using StaticArrays
using GeometryTypes


# Fix random seed.
Random.seed!(1234)
# Construct two camera matrices and parametrise two planar surfaces.
f = 50
image_width = 640
image_height = 480
𝐊₁ = @SMatrix [f 0 0 ;
               0 f 0 ;
               0 0 1 ]
𝐑₁ = SMatrix{3,3,Float64,9}(rotxyz(0, 1*(pi/180), 0))
𝐭₁ = [-300.0, 0.0, -50.0]

𝐊₂ = @SMatrix [f 0 0 ;
               0 f 0 ;
               0 0 1 ]

𝐑₂ = SMatrix{3,3,Float64,9}(rotxyz(0, -1*(pi/180), 0))
𝐭₂ = [300.0, 0.0, 5.0]

# # Normals and distance from origin
# 𝐧₁ = [0.0, 0.0, 1.0]
# d₁ = 55.0

z₁ = 100.0
z₂ = 200.0
x_range = -1000.0:1000.0
y_range = -1000.0:1000.0


N = 5000

𝒳₁ = generate_planar_points(-1000.0:1000.0, -1000.0:1000.0, z₁ , N)
𝒳₂ = generate_planar_points(-1500.0:1500.0, -1500.0:1500.0, z₂ , N)

scale = 20.0f0
x = Vec3f0(0); baselen = 0.2f0 * scale ; dirlen = 1f0 * scale
# create an array of differently colored boxes in the direction of the 3 axes
rectangles = [
    (HyperRectangle(Vec3f0(x), Vec3f0(dirlen, baselen, baselen)), RGBAf0(1,0,0,1)),
    (HyperRectangle(Vec3f0(x), Vec3f0(baselen, dirlen, baselen)), RGBAf0(0,1,0,1)),
    (HyperRectangle(Vec3f0(x), Vec3f0(baselen, baselen, dirlen)), RGBAf0(0,0,1,1))
]
meshes = map(GLNormalMesh, rectangles)
scene = mesh(merge(meshes))

𝐗₁ = reshape(reinterpret(Float64, 𝒳₁),3,N)
𝐗₂ = reshape(reinterpret(Float64, 𝒳₂),3,N)
scatter!(scene, 𝐗₁[1,:], 𝐗₁[2,:], 𝐗₁[3,:], markersize = 3, color = :red)
scatter!(scene, 𝐗₂[1,:], 𝐗₂[2,:], 𝐗₂[3,:], markersize = 3, color = :green)
# scatter!(scene, 𝐗₃[1,:], 𝐗₃[2,:], 𝐗₃[3,:], markersize = 3, color = :black)
# scatter!(scene, 𝐗₄[1,:], 𝐗₄[2,:], 𝐗₄[3,:], markersize = 3, color = :blue)


# Convention I: world_basis = (Vec(1.0, 0.0, 0.0), Vec(0.0, 1.0, 0.0), Vec(0.0, 0.0, 1.0))
world_basis = (Vec(1.0, 0.0, 0.0), Vec(0.0, 1.0, 0.0), Vec(0.0, 0.0, 1.0))
# Convention I: camera_basis = (Point(0.0, 0.0, 0.0), Vec(-1.0, 0.0, 0.0), Vec(0.0, -1.0, 0.0), Vec(0.0, 0.0, 1.0))
camera_basis = (Point(0.0, 0.0, 0.0), Vec(-1.0, 0.0, 0.0), Vec(0.0, -1.0, 0.0), Vec(0.0, 0.0, 1.0))
# Convention I: picture_basis = (Point(0.0, 0.0), Vec(1.0, 0.0), Vec(0.0, 1.0))
picture_basis = (Point(0.0, 0.0), Vec(-1.0, 0.0), Vec(0.0, -1.0))

camera₁ = Pinhole(image_width, image_height, f, camera_basis..., picture_basis...)
camera₂ = Pinhole(image_width, image_height, f, camera_basis..., picture_basis...)

relocate!(camera₁, 𝐑₁, 𝐭₁)
draw!(camera₁, scene)

relocate!(camera₂, 𝐑₂, 𝐭₂)
draw!(camera₂, scene)


𝐑₁′, 𝐭₁′ = ascertain_pose(camera₁, world_basis... )
𝐊₁′ = obtain_intrinsics(camera₁, CartesianSystem())

𝐑₂′, 𝐭₂′ = ascertain_pose(camera₂, world_basis... )
𝐊₂′ = obtain_intrinsics(camera₂, CartesianSystem())


𝐏₁ = construct(ProjectionMatrix(), 𝐊₁′, 𝐑₁′, 𝐭₁′)
𝐏₂ = construct(ProjectionMatrix(),𝐊₂′,𝐑₂′,𝐭₂′)

# Set of corresponding points.
ℳ₁ = project(camera₁,𝐏₁,𝒳₁)
ℳ₁ʹ= project(camera₂,𝐏₂,𝒳₁)

ℳ₂ = project(camera₁,𝐏₁,𝒳₂)
ℳ₂ʹ = project(camera₂,𝐏₂,𝒳₂)



# # Red
# ℳ₁ = project(camera₁, 𝐏₁, [SVector(𝐗₁...)])
# Base.display(ℳ₁)
#
# # Green
# ℳ₂ = project(camera₁, 𝐏₁, [SVector(𝐗₂...)])
# Base.display(ℳ₂)
#
# # Black
# ℳ₃ = project(camera₁, 𝐏₁, [SVector(𝐗₃...)])
# Base.display(ℳ₃)
#
# # Blue
# ℳ₄ = project(camera₁, 𝐏₁, [SVector(𝐗₄...)])
# Base.display(ℳ₄)

ℳ₁, ℳ₁ʹ = crop(HyperRectangle(Vec(0,0),Vec(200,200)), (ℳ₁, ℳ₁ʹ))
ℳ₂, ℳ₂ʹ = crop(HyperRectangle(Vec(300,300),Vec(200,200)), (ℳ₂, ℳ₂ʹ))

# Visualize the set of corresponding points
scene = Scene()
M₁ = reshape(reinterpret(Float64,ℳ₁),(2,length(ℳ₁)))
M₂ = reshape(reinterpret(Float64,ℳ₂),(2,length(ℳ₂)))
# M₃ = reshape(reinterpret(Float64,ℳ₃),(2,length(ℳ₃)))
# M₄ = reshape(reinterpret(Float64,ℳ₄),(2,length(ℳ₄)))

scatter!(scene,M₁[1,:], M₁[2,:], markersize = 5, color = :red, limits = FRect(0, 0, 640, 480.0))
scatter!(scene,M₂[1,:], M₂[2,:], markersize = 5, color = :green, limits = FRect(0, 0, 640, 480.0))
# scatter!(scene,M₃[1,:], M₃[2,:], markersize = 1, color = :black, limits = FRect(0, 0, 64, 48.0))
# scatter!(scene,M₄[1,:], M₄[2,:], markersize = 1, color = :blue, limits = FRect(0, 0, 64, 48.0))
