
@testset "Calibrated Rectificaton (Fusiello et al.)" begin

    world = create_rectification_world()
    @unpack points, planes = world
    plane = first(planes)
    𝐧 = normal(plane)
    d = distance(plane)

    # Verify that the default points lie on the plane
    𝛑 = push(𝐧, -d) # 𝛑 =[n -d]
    for 𝐱 in points
        @test isapprox(dot(𝛑, hom(𝐱)), 0.0; atol = 1e-14)
    end

    # Position two cameras in the scene
    camera₁ = Camera(image_type = AnalogueImage(coordinate_system = OpticalSystem()),
                     model = Pinhole(intrinsics = IntrinsicParameters(height = 400, width = 400)))
    𝐑₁ = SMatrix{3,3,Float64,9}(rotxyz(10*(pi/180), 0*(pi/180), 0*(pi/180)))
    ₁ = [0.0, 200.0, -300.0]
    camera₁ = relocate(camera₁, 𝐑₁, 𝐭₁)

    camera₂ = Camera(image_type = AnalogueImage(coordinate_system = OpticalSystem()),
                     model = Pinhole(intrinsics = IntrinsicParameters(height = 400, width = 400)))
    𝐑₂ = SMatrix{3,3,Float64,9}(rotxyz(-10*(pi/180), 5*(pi/180), 5*(pi/180)))
    𝐭₂ = [200.0, 0.0, -300.0]
    camera₂ = relocate(camera₂, 𝐑₂, 𝐭₂)

    # Project 3D points onto the cameras.
    aquire = AquireImage()
    ℳ = aquire(world, camera₁)
    ℳ′ = aquire(world, camera₂)

    # Rectify the camera pair so that all epipolar lines are horizontal.
    rectify = Rectify(camera₁, camera₂, FusielloCalibratedRectification())

    # Project the 3D points onto the rectified cameras.
    𝒩₁ = aquire(world, rectify.camera₁)
    𝒩₂ = aquire(world, rectify.camera₂)

    # Obtain the homography matrices that map the image plane of the unrectified
    # cameras to the image plane of corresponding rectified cameras.
    @unpack 𝐇₁, 𝐇₂ = rectify

    # Map points from the unrectified image to the rectified image.
    𝒪₁ = resolve(ℳ, rectify.𝐇₁)
    𝒪₂ = resolve(ℳ′, rectify.𝐇₂)

    # Projection the 3D points onto the rectified images, or mapping the
    # points from the unrectified image to the rectified image using the
    # special homography matrix should yield equivalent results.
    for i in eachindex(𝒩)
        @test isapprox(norm(𝒩₁[i] - 𝒪₁[i]), 0.0; atol = 1e-12)
        @test isapprox(norm(𝒩₂[i] - 𝒪₂[i]), 0.0; atol = 1e-12)
    end

    # The y-coordinates
    for i in eachindex(𝒩)
        @test isapprox(norm(last(𝒩₁[i]) - last(𝒩₂[i])), 0.0; atol = 1e-12)
    end

end
