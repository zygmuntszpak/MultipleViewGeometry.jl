
function construct_world(coordinate_system::AbstractCoordinateSystem)
    N = 3
    limits = range(-500; stop = 500, length = N)
    points = vec([Point3(x, y, 0.0) for x in limits, y in limits])
    planes = [EuclideanPlane3D(CartesianSystem(Point(0.0, 0.0, 0.0), Vec(-1.0, 0.0, 0.0), Vec(0.0, -1.0, 0.0), Vec(0.0, 0.0, 1.0)))]
    groups =  [IntervalAllotment(1:N^2)]
    PlanarWorld(coordinate_system = coordinate_system, points = points, planes = planes, groups = groups)
end

@testset "Canonical World Coordinate System & Camera " begin

    coordinate_system  = CartesianSystem(Point(0.0, 0.0, 0.0), Vec(1.0, 0.0, 0.0), Vec(0.0, 1.0, 0.0), Vec(0.0, 0.0, 1.0))
    world = construct_world(coordinate_system)
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

    𝐑, 𝐭 = ascertain_pose(camera₂, world_coordinate_system)
    pinholeₓ = Pinhole(intrinsics = IntrinsicParameters(width = 640, height = 480, focal_length = 100),
                       extrinsics = ExtrinsicParameters(world_coordinate_system))
    analogue_imageₓ = AnalogueImage(coordinate_system = OpticalSystem())
    cameraₓ = Camera(image_type = analogue_imageₓ, model = pinholeₓ)
    cameraₓ = relocate(cameraₓ, 𝐑, 𝐭)

    @test norm(extrinsics(camera₂).coordinate_system.𝐞₁ - extrinsics(cameraₓ).coordinate_system.𝐞₁) ≈ 0.0
    @test norm(extrinsics(camera₂).coordinate_system.𝐞₂ - extrinsics(cameraₓ).coordinate_system.𝐞₂) ≈ 0.0
    @test norm(extrinsics(camera₂).coordinate_system.𝐞₃ - extrinsics(cameraₓ).coordinate_system.𝐞₃) ≈ 0.0

    𝐄₂ = matrix(extrinsics(camera₂), world_coordinate_system)
    𝐄ₓ = matrix(extrinsics(cameraₓ), world_coordinate_system)
    @test norm(𝐄₂ - 𝐄ₓ) ≈ 0.0
end
