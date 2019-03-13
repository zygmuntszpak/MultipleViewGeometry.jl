using Makie
using MultipleViewGeometry, Test, Random
using MultipleViewGeometry.ModuleCostFunction
using MultipleViewGeometry.ModuleTypes
using MultipleViewGeometry.ModuleConstraints
using MultipleViewGeometry.ModuleConstruct
using MultipleViewGeometry.ModuleDraw
using MultipleViewGeometry.ModuleMove
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



# Generates a random point on the plane centered around a point on the plane
# that is closest to the origin.
function random_points_on_plane(𝐧::AbstractArray, d::Real, extent::Real, N::Int)
    # Generate vector 𝐰 on a plane through the origin with normal vector 𝐧.
    first(𝐧) == 0 ? 𝐰 = cross(𝐧,SVector(1.0,0.0,0.0)) : 𝐰 = cross(𝐧,SVector(0.0,0.0,1.0))
    points = Array{SVector{3,Float64},1}(undef,N)
    for n = 1:N
        # Rotate 𝐰 randomly around the axis 𝐧.
        θ = rand() * 2*pi
        𝐤 = 𝐧 / norm(𝐧)
        𝐯 = 𝐰*cos(θ) + cross(𝐤,𝐰)*sin(θ) + 𝐤*dot(𝐤,𝐰)*(1-cos(θ))
        # Scale the vector so that it lies in the interval [0, extent)
        𝐯 = (rand() * extent) * 𝐯
        # Translate the vector so that it lies on the plane parametrised by 𝐧 and d.
        𝐯 = 𝐯 + d*(𝐧/norm(𝐧)^2)
        points[n] = 𝐯
    end
    points
end

𝒳₁ = random_points_on_plane(𝐧₁,d₁, 20, 50)
𝒳₂ = random_points_on_plane(𝐧₂,d₂, 20, 50)
X₁ = reshape(reinterpret(Float64,𝒳₁),(3,length(𝒳₁)))
X₂ = reshape(reinterpret(Float64,𝒳₂),(3,length(𝒳₂)))

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
scatter!(scene, X₁[1,:],X₁[2,:], X₁[3,:], markersize = 3, color = :red)
scatter!(scene, X₂[1,:],X₂[2,:], X₂[3,:], markersize = 3, color = :blue)

world_basis = (Vec(1.0, 0.0, 0.0), Vec(0.0, 1.0, 0.0), Vec(0.0, 0.0, 1.0))
camera_basis = (Point(0.0, 0.0, 0.0), Vec(-1.0, 0.0, 0.0), Vec(0.0, -1.0, 0.0), Vec(0.0, 0.0, 1.0))
picture_basis = (Point(0.0, 0.0), Vec(-1.0, 0.0), Vec(0.0, -1.0))

camera₁ = Pinhole(image_width, image_height, f, camera_basis..., picture_basis...)
camera₂ = Pinhole(image_width, image_height, f, camera_basis..., picture_basis...)
relocate!(camera₁, 𝐑₁, 𝐭₁)
relocate!(camera₂, 𝐑₂, 𝐭₂)

draw!(camera₁, scene)
draw!(camera₂, scene)


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


M₁ = reshape(reinterpret(Float64,ℳ₁),(2,length(ℳ₁)))
M₂ = reshape(reinterpret(Float64,ℳ₂),(2,length(ℳ₂)))
M₁ʹ = reshape(reinterpret(Float64,ℳ₁ʹ),(2,length(ℳ₁ʹ)))
M₂ʹ = reshape(reinterpret(Float64,ℳ₂ʹ),(2,length(ℳ₂ʹ)))

# Visualize the set of corresponding points
scene = Scene()
scatter!(scene,M₁[1,:], M₁[2,:], markersize = 1, color = :red, limits = FRect(0, 0, 64, 48.0))
scatter!(scene,M₂[1,:], M₂[2,:], markersize = 1, color = :blue, limits = FRect(0, 0, 64, 48.0))

scene = Scene()
scatter!(scene,M₁ʹ[1,:], M₁ʹ[2,:], markersize = 1, color = :red, limits = FRect(0, 0, 64, 48.0))
scatter!(scene,M₂ʹ[1,:], M₂ʹ[2,:], markersize = 1, color = :blue, limits = FRect(0, 0, 64, 48.0))


# We will construct a pair of homography matrices and then construct a pair of
# projection matrices from the homographies.

𝐇₁ = construct(HomographyMatrix(),𝐊₁′,𝐑₁′,-𝐭₁′,𝐊₂′,𝐑₂′,-𝐭₂′,𝐧₁,d₁)
𝐇₂ = construct(HomographyMatrix(),𝐊₁′,𝐑₁′,-𝐭₁′,𝐊₂′,𝐑₂′,-𝐭₂′,𝐧₂,d₂)


context = ProjectionMatrices(HomographyMatrices(), Chojnacki(), TwoViews())
𝐎₁, 𝐎₂ = construct(context, (𝐇₁,𝐇₂))
𝐅₁ = construct(FundamentalMatrix(), 𝐏₁, 𝐏₂)
𝐅₂ = construct(FundamentalMatrix(), 𝐎₁, 𝐎₂)

𝐅₁ = 𝐅₁ / 𝐅₁[3,3]
𝐅₂ = 𝐅₂ / 𝐅₂[3,3]

# Use homographies to determine corresponding coordinates in the second image.
𝒪₁ = similar(ℳ₁ʹ)
𝒪₂ = similar(ℳ₂ʹ)
for n = 1:length(ℳ₁)
    𝐦 = hom(ℳ₁[n])
    𝐦ʹ = 𝐇₁*𝐦
    𝒪₁[n] = hom⁻¹(𝐦ʹ)
end
for n = 1:length(ℳ₂)
    𝐦 = hom(ℳ₂[n])
    𝐦ʹ = 𝐇₂*𝐦
    𝒪₂[n] = hom⁻¹(𝐦ʹ)
end

#TODO Test that 𝒪₁ .== ℳ₁ʹ and 𝒪₂ .== M₂ʹ
