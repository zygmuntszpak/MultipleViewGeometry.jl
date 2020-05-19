# TODO Place scene generation code into function to avoid dubplication
function generate_multi_planar_world()
    x₁ = 0.0
    x₁′ = 0.0
    y₁ = -1000.0
    y₁′ = 2000.0
    z₁ = -1000.0
    z₁′ = 1000.0
    points₁ = [Point3(rand(x₁:x₁′), rand(y₁:y₁′), rand(z₁:z₁′)) for n = 1:250]
    #EuclideanPlane3D(CartesianSystem(Point(0.0, 0.0, 0.0), Vec(1.0, 0.0, 0.0), Vec(0.0, 1.0, 0.0), Vec(0.0, 0.0, 1.0)))
    plane₁ = [Plane(Vec3(1.0, 0.0, 0.0), 0.0)]
    p₁ = [x₁, y₁, z₁]
    q₁ = [x₁, y₁′, z₁]
    r₁ = [x₁, y₁′, z₁′]
    s₁ = [x₁, y₁, z₁′]
    segment₁ = [p₁ => q₁ , q₁ => r₁ , r₁ => s₁ , s₁ => p₁]
    plane_segment₁ = PlaneSegment(first(plane₁), segment₁)

    x₂ = 0.0
    x₂′ = 3000.0
    y₂ = 2000.0
    y₂′ = 2000.0
    z₂ = -1000.0
    z₂′ = 1000.0
    points₂ = [Point3(rand(x₂:x₂′), rand(y₂:y₂′), rand(z₂:z₂′)) for n = 1:250]
    plane₂ = [Plane(Vec3(0.0, 1.0, 0.0), 2000.0)]
    p₂ = [x₂, y₂, z₂]
    q₂ = [x₂, y₂, z₂′]
    r₂ = [x₂′, y₂, z₂′]
    s₂ = [x₂′, y₂, z₂]
    segment₂ = [p₂ => q₂ , q₂ => r₂ , r₂ => s₂ , s₂ => p₂]
    plane_segment₂ = PlaneSegment(first(plane₂), segment₂)

    planes = vcat(plane_segment₁, plane_segment₂)
    points = vcat(points₁, points₂)
    groups = [IntervalAllotment(1:250), IntervalAllotment(251:500)]

    world = PrimitiveWorld(points = points, planes = planes, groups = groups)

    pinhole₁ = Pinhole(intrinsics = IntrinsicParameters(width = 640, height = 480, focal_length = 100.0))
    analogue_image₁ = AnalogueImage(coordinate_system = PlanarCartesianSystem())
    camera₁ = Camera(image_type = analogue_image₁, model = pinhole₁)
    𝐑₁ = SMatrix{3,3,Float64,9}(rotxyz(90*(pi/180), -130*(pi/180), 0*(pi/180)))
    𝐭₁ = [3000.0,0.0, 0.0]
    camera₁ = relocate(camera₁, 𝐑₁, 𝐭₁)

    pinhole₂ = Pinhole(intrinsics = IntrinsicParameters(width = 640, height = 480, focal_length = 100.0))
    analogue_image₂ = AnalogueImage(coordinate_system = PlanarCartesianSystem())
    camera₂ = Camera(image_type = analogue_image₂, model = pinhole₂)
    𝐑₂ = SMatrix{3,3,Float64,9}(rotxyz(90*(pi/180), -110*(pi/180), 0*(pi/180)))
    𝐭₂ = [4000.0,0.0, 0.0]
    camera₂ = relocate(camera₂, 𝐑₂, 𝐭₂)


    return world, camera₁, camera₂
end

function create_rectification_world()
    points = [Point3(100,100,500), Point3(100,-100,500) , Point3(-100,-100,500) , Point3(-100,100,500)]
    planes = [EuclideanPlane3D(CartesianSystem(Point(0.0, 0.0, 500.0), Vec(1.0, 0.0, 0.0), Vec(0.0, 1.0, 0.0), Vec(0.0, 0.0, 1.0)))]
    groups = [IntervalAllotment(1:4)]
    world = PlanarWorld(points = points, planes = planes, groups = groups)
    return world
end
