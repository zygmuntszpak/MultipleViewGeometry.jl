function construct_test_world()
    # Generate points on two planar surfaces
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

    planes = vcat(plane₁, plane₂)
    points = vcat(points₁, points₂)
    groups = [IntervalAllotment(1:250), IntervalAllotment(251:500)]
    world = PrimitiveWorld(points = points, planes = planes, groups = groups)
    return world
end

@testset "World Coordinate System Transformation" begin

    world = construct_test_world()

    pinhole₁ = Pinhole(intrinsics = IntrinsicParameters(width = 640, height = 480, focal_length = 100))
    analogue_image₁ = AnalogueImage(coordinate_system = OpticalSystem())
    camera₁ = MultipleViewGeometry.Camera(image_type = analogue_image₁, model = pinhole₁)
    𝐑₁ = SMatrix{3,3,Float64,9}(rotxyz(90*(pi/180), -110*(pi/180), 0*(pi/180)))
    𝐭₁ = [3000.0,0.0, 0.0]
    camera₁ = relocate(camera₁, 𝐑₁, 𝐭₁)

    pinhole₂ = Pinhole(intrinsics = IntrinsicParameters(width = 640, height = 480, focal_length = 100))
    analogue_image₂ = AnalogueImage(coordinate_system = OpticalSystem())
    camera₂ = MultipleViewGeometry.Camera(image_type = analogue_image₂, model = pinhole₂)
    𝐑₂ = SMatrix{3,3,Float64,9}(rotxyz(90*(pi/180), -110*(pi/180), 0*(pi/180)))
    𝐭₂ = [4000.0,0.0, 0.0]
    camera₂ = relocate(camera₂, 𝐑₂, 𝐭₂)


    aquire = AquireImage()

    # Project 3D points onto the cameras.
    ℳ = aquire(world, camera₁)
    ℳ′ = aquire(world, camera₂)

    cameraₐ  = deepcopy(camera₁)
    cameraᵦ  = deepcopy(camera₂)
    world₂ = deepcopy(world)

    default_world_system = CartesianSystem(Point(0.0, 0.0, 0.0), Vec(1.0, 0.0, 0.0), Vec(0.0, 1.0, 0.0), Vec(0.0, 0.0, 1.0))
    alternative_world_system = extrinsics(camera₁).coordinate_system

    transformation_context = WorldSystemTransformation(CoordinateTransformation(source = default_world_system, target = alternative_world_system))
    cameraₐ = transformation_context(cameraₐ)
    cameraᵦ = transformation_context(cameraᵦ)
    world₂ = transformation_context(world₂)

    # Project transformed 3D points onto the cameras.
    ℳₐ = aquire(world₂, cameraₐ)
    ℳᵦ′ = aquire(world₂, cameraᵦ)


    # Verify that the coordinates of the image points are the same irrespective
    # of the choice of the world coordinate system.
    for couple in zip(ℳ, ℳₐ)
        @test isapprox(norm(first(couple)-last(couple)), 0.0; atol = 1e-10)
    end
    for couple in zip(ℳ′, ℳᵦ′)
        @test isapprox(norm(first(couple)-last(couple)), 0.0; atol = 1e-10)
    end

    # Verify that the original 3D points lie on their corresponding planes.
    points3D = world.points
    planes3D = world.planes
    groups = world.groups
    for (i, plane3D) in enumerate(planes3D)
        @unpack interval = groups[i]
        subset = points3D[interval]
        for pt in subset
            @test on_plane(pt, plane3D; tol = 1e-10)
        end
    end

    # Verify that the transformed 3D points lie on their corresponding transformed planes.
    points3Dᵦ = world₂.points
    planes3Dᵦ = world₂.planes
    groups = world₂.groups
    for (i, plane3Dᵦ) in enumerate(planes3Dᵦ)
        @unpack interval = groups[i]
        subset = points3Dᵦ[interval]
        for pt in subset
            @test on_plane(pt, plane3Dᵦ; tol = 1e-10)
        end
    end

    # Verify that the fundamental matrices are the same irrespective of how we choose
    # the world coordinate system.
    𝐅₁ = matrix(FundamentalMatrix(camera₁, camera₂))
    𝐅₁ = 𝐅₁ / norm(𝐅₁)
    𝐅₂ = matrix(FundamentalMatrix(cameraₐ, cameraᵦ))
    𝐅₂ = 𝐅₂ / norm(𝐅₂)
    @test norm(𝐅₁ - 𝐅₂) < 1e-15


    # Verify that the homography matrices are the same irrespective of how we choose
    # the world coordinate system.
    ℋ = matrices(HomographyMatrices(camera₁, camera₂, world.planes))
    𝐇₁ = ℋ[1] / norm(ℋ[1])
    𝐇₂ = ℋ[2] / norm(ℋ[2])

    ℋ₂ = matrices(HomographyMatrices(cameraₐ, cameraᵦ, world₂.planes))
    𝐇ₐ = ℋ₂[1] / norm(ℋ₂[1])
    𝐇ᵦ = ℋ₂[2] / norm(ℋ₂[2])

    @test min(norm(𝐇₁ - 𝐇ₐ), norm(𝐇₁ + 𝐇ₐ)) < 1e-15
    @test min(norm(𝐇₂ - 𝐇ᵦ), norm(𝐇₂ + 𝐇ᵦ)) < 1e-15

    # Verify that the points are triangulated correctly for both choices of
    # world coordinate systems.
    𝒞₁ = Observations(tuple(ℳ,ℳ′))
    triangulate =  Triangulate(DirectLinearTriangulation())
    estimated_points₁ = triangulate(camera₁, camera₂, 𝒞₁)
    # Verify that the triangulated points are close to the true points.
    for couple in zip(estimated_points₁, world.points)
        @test isapprox(norm(first(couple)-last(couple)), 0.0; atol = 1e-7)
    end

    𝒞₂ = Observations(tuple(ℳₐ,ℳᵦ′))
    triangulate =  Triangulate(DirectLinearTriangulation())
    estimated_points₂ = triangulate(cameraₐ, cameraᵦ, 𝒞₂)
    # Verify that the triangulated points are close to the true points.
    for couple in zip(estimated_points₂, world₂.points)
        @test isapprox(norm(first(couple)-last(couple)), 0.0; atol = 1e-7)
    end
end
