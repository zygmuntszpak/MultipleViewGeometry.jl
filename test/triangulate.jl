
function create_triangulation_test_world()
    x₁ = 0.0
    x₁′ = 0.0
    y₁ = -1000.0
    y₁′ = 2000.0
    z₁ = -1000.0
    z₁′ = 1000.0
    points₁ = [Point3(rand(x₁:x₁′), rand(y₁:y₁′), rand(z₁:z₁′)) for n = 1:250]
    plane₁ = [Plane(Vec3(1.0, 0.0, 0.0), 0)]
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
    plane₂ = [Plane(Vec3(0.0, 1.0, 0.0), 2000)]
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
end
world = create_triangulation_test_world()
@unpack points, planes = world

@testset "Two-View Triangulation" begin
    plane = first(planes)
    𝐧 = normal(plane)
    d = distance(plane)

    𝛑 = push(𝐧, -d) # 𝛑 =[n -d]
    for 𝐱 in points[1:250]
        @test isapprox(dot(𝛑, hom(𝐱)), 0.0; atol = 1e-14)
    end

    pinhole₁ = Pinhole(intrinsics = IntrinsicParameters(width = 640, height = 480, focal_length = 100))
    analogue_image₁ = AnalogueImage(coordinate_system = OpticalSystem())
    camera₁ = MultipleViewGeometry.Camera(image_type = analogue_image₁, model = pinhole₁)
    𝐑₁ = SMatrix{3,3,Float64,9}(rotxyz(40*(pi/180), -130*(pi/180), 10*(pi/180)))
    𝐭₁ = [3000.0,0.0, 0.0]
    camera₁ = relocate(camera₁, 𝐑₁, 𝐭₁)

    pinhole₂ = Pinhole(intrinsics = IntrinsicParameters(width = 640, height = 480, focal_length = 100))
    analogue_image₂ = AnalogueImage(coordinate_system = OpticalSystem())
    camera₂ = MultipleViewGeometry.Camera(image_type = analogue_image₂, model = pinhole₂)
    𝐑₂ = SMatrix{3,3,Float64,9}(rotxyz(90*(pi/180), -110*(pi/180), 0*(pi/180)))
    𝐭₂ = [4000.0,0.0, 0.0]
    camera₂ = relocate(camera₂, 𝐑₂, 𝐭₂)

    aquire = AquireImage()
    ℳ = aquire(world, camera₁)
    ℳ′ = aquire(world, camera₂)

    𝒞 = Observations(tuple(ℳ,ℳ′))

    triangulate =  Triangulate(DirectLinearTriangulation())
    estimated_points = triangulate(camera₁, camera₂, 𝒞)

    # Verify that the triangulated points are close to the true points.
    for couple in zip(estimated_points, points)
        @test isapprox(norm(first(couple)-last(couple)), 0.0; atol = 1e-7)
    end
end
