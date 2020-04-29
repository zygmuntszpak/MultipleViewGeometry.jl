using MultipleViewGeometry
using StaticArrays
using Parameters
using GeometryBasics
import Makie


function construct_calibration_world()
    coordinate_system  = CartesianSystem(Point(0.0, 0.0, 0.0), Vec(1.0, 0.0, 0.0), Vec(0.0, 1.0, 0.0), Vec(0.0, 0.0, 1.0))
    N = 10
    limits = range(-500; stop = 500, length = N)
    points = vec([Point3(x, y, 0.0) for x in limits, y in limits])
    planes = [EuclideanPlane3D(CartesianSystem(Point(0.0, 0.0, 0.0), Vec(1.0, 0.0, 0.0), Vec(0.0, 1.0, 0.0), Vec(0.0, 0.0, 1.0)))]
    groups =  [IntervalAllotment(1:N^2)]
    PlanarWorld(coordinate_system = coordinate_system, points = points, planes = planes, groups = groups)
end

world = construct_calibration_world()

pinhole₁ = Pinhole(intrinsics = IntrinsicParameters(width = 640, height = 480, focal_length = 100))
analogue_image₁ = AnalogueImage(coordinate_system = OpticalSystem())
camera₁ = Camera(image_type = analogue_image₁, model = pinhole₁)
𝐑₁ = SMatrix{3,3,Float64,9}(rotxyz(180*(pi/180), 0*(pi/180), 0*(pi/180)))
𝐭₁ = [-500.0,0.0, 3000.0]
camera₁ = relocate(camera₁, 𝐑₁, 𝐭₁)

pinhole₂ = Pinhole(intrinsics = IntrinsicParameters(width = 640, height = 480, focal_length = 100))
analogue_image₂ = AnalogueImage(coordinate_system = OpticalSystem())
camera₂ = Camera(image_type = analogue_image₂, model = pinhole₂)
𝐑₂ = SMatrix{3,3,Float64,9}(rotxyz(180*(pi/180), 10*(pi/180), 10*(pi/180)))
𝐭₂ = [500.0,0.0, 3000.0]
camera₂ = relocate(camera₂, 𝐑₂, 𝐭₂)

pinhole₃ = Pinhole(intrinsics = IntrinsicParameters(width = 640, height = 480, focal_length = 100))
analogue_image₃ = AnalogueImage(coordinate_system = OpticalSystem())
camera₃ = Camera(image_type = analogue_image₃, model = pinhole₃)
𝐑₃ = SMatrix{3,3,Float64,9}(rotxyz(180*(pi/180), -10*(pi/180), 10*(pi/180)))
𝐭₃ = [0.0,0.0, 2500.0]
camera₃ = relocate(camera₃, 𝐑₃, 𝐭₃)

pinhole₄ = Pinhole(intrinsics = IntrinsicParameters(width = 640, height = 480, focal_length = 100))
analogue_image₄ = AnalogueImage(coordinate_system = OpticalSystem())
camera₄ = Camera(image_type = analogue_image₄, model = pinhole₄)
𝐑₄ = SMatrix{3,3,Float64,9}(rotxyz(180*(pi/180), 0*(pi/180), -10*(pi/180)))
𝐭₄ = [0.0,500.0, 3500.0]
camera₄ = relocate(camera₄, 𝐑₄, 𝐭₄)


visualize =  VisualizeWorld(; visual_properties = MakieVisualProperties(scale = 150, markersize = 25))
cameras = [camera₁, camera₂, camera₃, camera₄]
visualize(world, cameras)
@unpack scene = visualize
display(scene)


calibrate = CalibrateCamera()
calibrate(world, cameras)


@unpack points = world
ℳ′ = [Point(p[1], p[2]) for p in points]

aquire = AquireImage()
# Determine projections of the 3D points in each camera view.
𝓜 = [aquire(world, camera) for camera in cameras]

ℋ = [fit_homography(ℳ′,ℳ,  DirectLinearTransform()) for ℳ in 𝓜]
𝐇₁ = matrix(ℋ[1])
𝐇₁ = 𝐇₁ / norm(𝐇₁)

𝐇₂ = matrix(ℋ[2])
𝐇₂ = 𝐇₂ / norm(𝐇₂)

# world_plane = collect(reshape(reinterpret(Float64,ℳ′),(2,100)))
# image_1 = collect(reshape(reinterpret(Float64,𝓜[1]),(2,100)))
# image_2 = collect(reshape(reinterpret(Float64,𝓜[2]),(2,100)))
# image_3 = collect(reshape(reinterpret(Float64,𝓜[3]),(2,100)))
# image_4 = collect(reshape(reinterpret(Float64,𝓜[4]),(2,100)))
# file = matopen("matfile.mat", "w")
# write(file, "world_plane", world_plane)
# write(file, "image_1", image_1)
# write(file, "image_2", image_2)
# write(file, "image_3", image_3)
# write(file, "image_4", image_4)
# close(file)



# aquire = AquireImage()
# pts₁ = aquire(world, camera₁)
# pts₂ = aquire(world, camera₂)

# # Estimate a homography matrix between corresponding points
# homography = fit_homography(pts₁, pts₂, DirectLinearTransform())
# 𝐇 = matrix(homography)
#
# s₁ = [vec2antisym(hom(pts₂[i])) * 𝐇  * hom(pts₁[i]) for i = 1:length(pts₁)]
#
# for 𝐱 in s₁
#     @test all(isapprox.(𝐱, [0.0, 0.0, 0.0] ; atol = 1e-12))
# end

# Estimate a homography matrix between corresponding points
# fundamental_matrix = fit_fundamental_matrix(pts₁, pts₂, DirectLinearTransform())
# 𝐅 = matrix(fundamental_matrix)
