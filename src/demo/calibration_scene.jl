using MultipleViewGeometry
using StaticArrays
using Parameters
using GeometryBasics
using Setfield
using LinearAlgebra
#import Makie


function construct_calibration_world()
    coordinate_system  = CartesianSystem(Point(0.0, 0.0, 0.0), Vec(1.0, 0.0, 0.0), Vec(0.0, 1.0, 0.0), Vec(0.0, 0.0, 1.0))
    N = 3
    limits = range(-500; stop = 500, length = N)
    points = vec([Point3(x, y, 0.0) for x in limits, y in limits])
    planes = [EuclideanPlane3D(CartesianSystem(Point(0.0, 0.0, 0.0), Vec(1.0, 0.0, 0.0), Vec(0.0, 1.0, 0.0), Vec(0.0, 0.0, 1.0)))]
    groups =  [IntervalAllotment(1:N^2)]
    PlanarWorld(coordinate_system = coordinate_system, points = points, planes = planes, groups = groups)
end

world = construct_calibration_world()
world_coordinate_system = world.coordinate_system

pinhole₁ = Pinhole(intrinsics = IntrinsicParameters(width = 640, height = 480, focal_length = 100),
                   extrinsics = ExtrinsicParameters(world_coordinate_system))
analogue_image₁ = AnalogueImage(coordinate_system = OpticalSystem())
camera₁ = Camera(image_type = analogue_image₁, model = pinhole₁)
𝐑₁ = SMatrix{3,3,Float64,9}(rotxyz(0*(pi/180), 0*(pi/180), 0*(pi/180)))
𝐭₁ = [-500.0,0.0, -3000.0]
camera₁ = relocate(camera₁, 𝐑₁, 𝐭₁)

pinhole₂ = Pinhole(intrinsics = IntrinsicParameters(width = 640, height = 480, focal_length = 100),
                   extrinsics = ExtrinsicParameters(world_coordinate_system))
analogue_image₂ = AnalogueImage(coordinate_system = OpticalSystem())
camera₂ = Camera(image_type = analogue_image₂, model = pinhole₂)
𝐑₂ = SMatrix{3,3,Float64,9}(rotxyz(0*(pi/180), 10*(pi/180), 10*(pi/180)))
𝐭₂ = [500.0, 0.0, -3000.0]
camera₂ = relocate(camera₂, 𝐑₂, 𝐭₂)

pinhole₃ = Pinhole(intrinsics = IntrinsicParameters(width = 640, height = 480, focal_length = 100),
                   extrinsics = ExtrinsicParameters(world_coordinate_system))
analogue_image₃ = AnalogueImage(coordinate_system = OpticalSystem())
camera₃ = Camera(image_type = analogue_image₃, model = pinhole₃)
𝐑₃ = SMatrix{3,3,Float64,9}(rotxyz(0*(pi/180), -10*(pi/180), 10*(pi/180)))
𝐭₃ = [0.0, 0.0, -2500.0]
camera₃ = relocate(camera₃, 𝐑₃, 𝐭₃)

pinhole₄ = Pinhole(intrinsics = IntrinsicParameters(width = 640, height = 480, focal_length = 100),
                   extrinsics = ExtrinsicParameters(world_coordinate_system))
analogue_image₄ = AnalogueImage(coordinate_system = OpticalSystem())
𝐑₄ = determine_rotation(0*(pi/180), 0*(pi/180), -10*(pi/180), world_coordinate_system)
𝐭₄ = [0.0, 500.0, -3500.0]
camera₄ = relocate(camera₄, 𝐑₄, 𝐭₄)

visualize =  VisualizeWorld(; visual_properties = MakieVisualProperties(scale = 150, markersize = 25))
cameras = [camera₁, camera₂, camera₃, camera₄]
visualize(world, cameras)
@unpack scene = visualize
display(scene)

# Determine projections of the 3D points in each camera view.
aquire = AquireImage()
𝓜 = [aquire(world, camera) for camera in cameras]
Λ = SMatrix{2,2}(I) * 0.005
# Add slight noise.
𝓞 = [apply_noise(𝓜[i], Λ) for i = 1:length(𝓜)]

calibrate = CalibrateCamera()
#calibrated_cameras = calibrate(world, cameras)
calibrated_cameras = calibrate(world, 𝓞)

visualize(world, calibrated_cameras)
@unpack scene = visualize
display(scene)
