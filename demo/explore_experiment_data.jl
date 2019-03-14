using Makie
using MultipleViewGeometry, Test, Random
using MultipleViewGeometry.ModuleCostFunction
using MultipleViewGeometry.ModuleTypes
using MultipleViewGeometry.ModuleConstraints
using MultipleViewGeometry.ModuleConstruct
using MultipleViewGeometry.ModuleDraw
using MultipleViewGeometry.ModuleMove
using MultipleViewGeometry.ModuleSyntheticData
using MultipleViewGeometry.ModuleNoise
using MultipleViewGeometry.ModuleExperiment
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

# Normals and (negative) distance from origin
𝐧₁ = [0.0, 0.0, 1.0]
d₁ = -100.0

𝐧₂ = [0.0, 0.0, 1.0]
d₂ = -200.0

z₁ = -d₁
z₂ = -d₂
x_range = -1000.0:1000.0
y_range = -1000.0:1000.0

N = 5000

𝒳₁ = generate_planar_points(-1000.0:1000.0, -1000.0:1000.0, z₁ , N)
𝒳₂ = generate_planar_points(-1500.0:1500.0, -1500.0:1500.0, z₂ , N)

𝐇₁ = construct(HomographyMatrix(),𝐊₁′,𝐑₁′,𝐭₁′,𝐊₂′,𝐑₂′,𝐭₂′, 𝐧₁, d₁)
𝐇₂ = construct(HomographyMatrix(),𝐊₁′,𝐑₁′,𝐭₁′,𝐊₂′,𝐑₂′,𝐭₂′, 𝐧₂, d₂)

𝐏₁ = construct(ProjectionMatrix(), 𝐊₁′, 𝐑₁′, 𝐭₁′)
𝐏₂ = construct(ProjectionMatrix(),𝐊₂′,𝐑₂′,𝐭₂′)

# Set of corresponding points.
ℳ₁ = project(camera₁,𝐏₁,𝒳₁)
ℳ₁ʹ= project(camera₂,𝐏₂,𝒳₁)

ℳ₂ = project(camera₁,𝐏₁,𝒳₂)
ℳ₂ʹ = project(camera₂,𝐏₂,𝒳₂)


# Discard corresponding points which fall outside specific rectangular regions in the first image.
𝒪₁, 𝒪₁ʹ = crop(HyperRectangle(Vec(0,0),Vec(200,200)), (ℳ₁, ℳ₁ʹ))
𝒪₂, 𝒪₂ʹ = crop(HyperRectangle(Vec(300,300),Vec(200,200)), (ℳ₂, ℳ₂ʹ))

trials = Vector{ExperimentTrial}(undef,10)
for t = 1:10
    𝒫₁, 𝒫₁ʹ = perturb(GaussianNoise(), 1, (𝒪₁, 𝒪₁ʹ))
    𝒫₂, 𝒫₂ʹ = perturb(GaussianNoise(), 1, (𝒪₂, 𝒪₂ʹ))
    trial = ExperimentTrial( (𝒪₁, 𝒪₁ʹ),  (𝒫₁, 𝒫₁ʹ), (ℳ₁, ℳ₁ʹ),  @SMatrix zeros(3,3))
    trials[t] = trial
end

condition = ExperimentCondition("σ = 1", trials)
participant₁ = Participant(DirectLinearTransform(), Dict(condition.description => condition))
participant₂ = Participant(BundleAdjustment(DirectLinearTransform(), 5, 1e-10), Dict(condition.description => condition))

experiment = Experiment("Increasing Noise Level", HomographyMatrix(), [participant₁, participant₂], ReprojectionError())

#experiment = construct_experiment(2,10,1:3)

experiment = construct_experiment(PlanarScene(2, 10, tuple(100:200, 100:200)),10,1:3)

z = conduct_experiment(experiment, ReprojectionError())


estimate(HomographyMatrix(), DirectLinearTransform(), p)

tabulate(z)

p = ([rand(2,4), rand(2,4)], [rand(2,4), rand(2,4)])


z.participant_results

r₁ = assess(ReprojectionError(), HomographyMatrix(), 𝐇₁, (𝒪₁, 𝒪₁ʹ))

@time assess(ReprojectionError(), HomographyMatrix(), 𝐇₁, (𝒪₁, 𝒪₁ʹ))

@code_warntype assess(ReprojectionError(), HomographyMatrix(), 𝐇₁, (𝒪₁, 𝒪₁ʹ))

@enter assess(ReprojectionError(), HomographyMatrix(), 𝐇₁, (𝒪₁, 𝒪₁ʹ))

J1 = r₁.jacobian

J2 = r₁.jacobian

zeros()

# hom⁻¹(𝐇₁*hom(ℳ₁[1]))
#
# ℳ₁ʹ[1]
#
# 𝐧₁'*𝒳₁[1] + d₁
#
# hom⁻¹(𝐇₂*hom(ℳ₂[1]))
#
# ℳ₂ʹ[1]
#
# 𝐧₂'*𝒳₂[1] + d₂
