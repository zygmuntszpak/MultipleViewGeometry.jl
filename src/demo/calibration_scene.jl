using MultipleViewGeometry
using StaticArrays
using Parameters
using GeometryBasics
import Makie


function construct_calibration_world()
    coordinate_system  = CartesianSystem(Point(0.0, 0.0, 0.0), Vec(1.0, 0.0, 0.0), Vec(0.0, 1.0, 0.0), Vec(0.0, 0.0, 1.0))
    N = 10
    limits = range(-1000; stop = 1000, length = N)
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
𝐑₂ = SMatrix{3,3,Float64,9}(rotxyz(180*(pi/180), 0*(pi/180), 0*(pi/180)))
𝐭₂ = [500.0,0.0, 3000.0]
camera₂ = relocate(camera₂, 𝐑₂, 𝐭₂)

visualize =  VisualizeWorld(; visual_properties = MakieVisualProperties(scale = 150, markersize = 25))
visualize(world, [camera₁, camera₂])
@unpack scene = visualize
display(scene)


aquire = AquireImage()
pts₁ = aquire(world, camera₁)
pts₂ = aquire(world, camera₂)

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
fundamental_matrix = fit_fundamental_matrix(pts₁, pts₂, DirectLinearTransform())
𝐅 = matrix(fundamental_matrix)
