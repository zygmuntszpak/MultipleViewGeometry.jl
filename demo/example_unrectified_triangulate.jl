using MultipleViewGeometry, Base.Test
using MultipleViewGeometry.ModuleCostFunction
using MultipleViewGeometry.ModuleTypes
using MultipleViewGeometry.ModuleConstraints
using MultipleViewGeometry.ModuleConstruct
using MultipleViewGeometry.ModuleDraw
using MultipleViewGeometry.ModuleTriangulation
using BenchmarkTools, Compat
using StaticArrays
using  MAT, Plots

# Load MATLAB matrices that represent a pair of images and that contain
# a set of manually matched corresponding points.
file = matopen("./data/teapot.mat")
X = read(file,"pts3D")
X = X[:,1:10:end]
close(file)

# Fix random seed.
srand(1234)
plotlyjs()

#𝒳 = [Point3DH(x,y,z,1.0) for x=-1:0.5:10 for y=-1:0.5:10 for z=2:-0.1:1]

𝒳 = [Point3DH(X[1,i],X[2,i],X[3,i],1.0) for i = 1:size(X,2)]

#X = reinterpret(Float64,map(SVector{4,Float64},𝒳),(4,length(𝒳)))
#Z = reinterpret(SVector{4,Float64},(3000,4))
# reinterpret(SVector{4,Float64}, X, (size(X,2),))#
# reinterpret(SVector{4,Float64}, Z, (size(Z,2),))#
# reinterpret(SVector{4,Float64}, X, (3000,1))

# Intrinsic and extrinsic parameters of camera one.
𝐊₁ = eye(3)
𝐑₁ = eye(3)
𝐭₁ = [-250.0, 0.0, 2500.0]
𝐜₁ = [0.0, 0.0, 0.0]
𝐄₁ = [𝐑₁ 𝐭₁]

# Intrinsic and extrinsic parameters of camera two.
𝐊₂ = eye(3)
𝐑₂ =  eye(3)
𝐭₂ = [250.0,   0.0, 2500.0]
𝐜₂ = [0.0, 0.0, 0.0]
𝐄₂ = [𝐑₂ 𝐭₂]

# Camera projection matrices.
𝐏₁ = construct(ProjectionMatrix(),𝐊₁,𝐑₁,𝐭₁)
𝐏₂ = construct(ProjectionMatrix(),𝐊₂,𝐑₂,𝐭₂)

# Set of corresponding points.
ℳ = project(Pinhole(),𝐏₁,𝒳)
ℳʹ = project(Pinhole(),𝐏₂,𝒳)

# Convert the array of Point3DH into a 4 x N matrix to simplify
# the plotting of the data points.
#X = reinterpret(Float64,map(SVector{4,Float64},𝒳),(4,length(𝒳)))

# Visualise the data points
p1 = Plots.plot(X[1,:],X[2,:],X[3,:],seriestype = :scatter, ms=1,grid = false, box = :none, legend = false)
draw!(WorldCoordinateSystem3D(), 450, p1)
draw!(Camera3D(), 𝐊₁, 𝐑₁, 𝐭₁, 250, p1)
draw!(Camera3D(), 𝐊₂, 𝐑₂, 𝐭₂, 250, p1)


# Plot the projections of the point cloud in the image pair.
p2 =  Plots.plot();
for n = 1:length(ℳ)
    m = ℳ[n]
    Plots.plot!([m[1]],[m[2]], grid = false, box = :none, legend = false,
                    seriestype = :scatter, ms = 2, markercolor=:Black,
                    aspect_ratio = :equal)
end

p3 =  Plots.plot();
for n = 1:length(ℳʹ)
    mʹ = ℳʹ[n]
    Plots.plot!([mʹ[1]],[mʹ[2]], grid = false, box = :none, legend = false,
                    seriestype = :scatter, ms = 2, markercolor=:Black,
                    aspect_ratio = :equal)
end

# Visualise the 3D point cloud, as well as the projected images.
l = @layout [ a; [b c] ]
p4 = Plots.plot(p1,p2, p3, layout = l)


# 𝒴 = triangulate(DirectLinearTransform(),𝐏₁,𝐏₂,(ℳ,ℳʹ))
#
# Y = reinterpret(Float64,map(SVector{4,Float64},𝒴),(4,length(𝒴)))
# # Visualise the data points
# p5 = Plots.plot(Y[1,:],Y[2,:],Y[3,:],seriestype = :scatter, ms=1,grid = false, box = :none, legend = false)
# draw!(WorldCoordinateSystem3D(), 450, p1)
# draw!(Camera3D(), 𝐊₁, 𝐑₁, 𝐭₁, 250, p1)
# draw!(Camera3D(), 𝐊₂, 𝐑₂, 𝐭₂, 250, p1)

# Estimate of the fundamental matrix and the true fundamental matrix.
𝐅 = estimate(FundamentalMatrix(), DirectLinearTransform(), (ℳ, ℳʹ))
𝐄 = construct(EssentialMatrix(), 𝐅,  𝐊₁, 𝐊₂)

𝐏₁, 𝐏₂ = construct(ProjectionMatrix(), 𝐄, (ℳ, ℳʹ))

𝒴 = triangulate(DirectLinearTransform(),𝐏₁,𝐏₂,(ℳ,ℳʹ))
Y = reinterpret(Float64,map(SVector{4,Float64},𝒴),(4,length(𝒴)))
# Visualise the data points
p5 = Plots.plot(Y[1,:],Y[2,:],Y[3,:],seriestype = :scatter, ms=1,grid = false, box = :none, legend = false)
draw!(WorldCoordinateSystem3D(), 450/1000, p1)
# draw!(Camera3D(), 𝐊₁, 𝐏₁[1:3,1:3], 𝐏₁[:,4], 250/1000, p1)
# draw!(Camera3D(), 𝐊₂, 𝐏₂[1:3,1:3], 𝐏₂[:,4], 250/1000, p1)
# draw!(Camera3D(), 𝐊₁, 𝐑₁, 𝐭₁, 250, p1)
# draw!(Camera3D(), 𝐊₂, 𝐑₂, 𝐭₂, 250, p1)
