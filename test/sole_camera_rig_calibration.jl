function construct_calibration_world()
    coordinate_system  = CartesianSystem(Point(0.0, 0.0, 0.0), Vec(1.0, 0.0, 0.0), Vec(0.0, 1.0, 0.0), Vec(0.0, 0.0, 1.0))
    N = 3
    limits = range(-500; stop = 500, length = N)
    points = vec([Point3(x, y, 0.0) for x in limits, y in limits])
    planes = [EuclideanPlane3D(CartesianSystem(Point(0.0, 0.0, 0.0), Vec(1.0, 0.0, 0.0), Vec(0.0, 1.0, 0.0), Vec(0.0, 0.0, 1.0)))]
    groups =  [IntervalAllotment(1:N^2)]
    PlanarWorld(coordinate_system = coordinate_system, points = points, planes = planes, groups = groups)
end


@testset "Sole Camera Calibration" begin
    world = construct_calibration_world()

    pinhole₁ = Pinhole(intrinsics = IntrinsicParameters(width = 640, height = 480, focal_length = 100))
    analogue_image₁ = AnalogueImage(coordinate_system = OpticalSystem())
    camera₁ = Camera(image_type = analogue_image₁, model = pinhole₁)
    𝐑₁ = SMatrix{3,3,Float64,9}(rotxyz(0*(pi/180), 0*(pi/180), 0*(pi/180)))
    𝐭₁ = [-500.0,0.0, -3000.0]
    camera₁ = relocate(camera₁, 𝐑₁, 𝐭₁)

    pinhole₂ = Pinhole(intrinsics = IntrinsicParameters(width = 640, height = 480, focal_length = 100))
    analogue_image₂ = AnalogueImage(coordinate_system = OpticalSystem())
    camera₂ = Camera(image_type = analogue_image₂, model = pinhole₂)
    𝐑₂ = SMatrix{3,3,Float64,9}(rotxyz(0*(pi/180), 10*(pi/180), 10*(pi/180)))
    𝐭₂ = [500.0,0.0, -3000.0]
    camera₂ = relocate(camera₂, 𝐑₂, 𝐭₂)

    pinhole₃ = Pinhole(intrinsics = IntrinsicParameters(width = 640, height = 480, focal_length = 100))
    analogue_image₃ = AnalogueImage(coordinate_system = OpticalSystem())
    camera₃ = Camera(image_type = analogue_image₃, model = pinhole₃)
    𝐑₃ = SMatrix{3,3,Float64,9}(rotxyz(0*(pi/180), -10*(pi/180), 10*(pi/180)))
    𝐭₃ = [0.0,0.0, -2500.0]
    camera₃ = relocate(camera₃, 𝐑₃, 𝐭₃)

    pinhole₄ = Pinhole(intrinsics = IntrinsicParameters(width = 640, height = 480, focal_length = 100))
    analogue_image₄ = AnalogueImage(coordinate_system = OpticalSystem())
    camera₄ = Camera(image_type = analogue_image₄, model = pinhole₄)
    𝐑₄ = SMatrix{3,3,Float64,9}(rotxyz(0*(pi/180), 0*(pi/180), -10*(pi/180)))
    𝐭₄ = [0.0,500.0, -3500.0]
    camera₄ = relocate(camera₄, 𝐑₄, 𝐭₄)

    reference_cameras = [camera₁, camera₂, camera₃, camera₄]

    calibrate = CalibrateCamera()
    # Determine projections of the 3D points in each camera view.
    aquire = AquireImage()
    𝓜 = [aquire(world, camera) for camera in reference_cameras]
    # Test on noiseless data
    calibrated_cameras = calibrate(world, 𝓜)
    for i = 1:4
        ref_camᵢ = reference_cameras[i]
        calib_camᵢ = calibrated_cameras[i]

        ref_intrinsics = intrinsics(ref_camᵢ)
        calib_intrinsics = intrinsics(calib_camᵢ)
        # Verify intrinsics are close to the ground truth.
        @test isapprox(ref_intrinsics.focal_length, calib_intrinsics.focal_length; atol = 1e-7)
        @test isapprox(ref_intrinsics.skewedness, calib_intrinsics.skewedness; atol = 1e-7)
        @test isapprox(ref_intrinsics.scale_x, calib_intrinsics.scale_x; atol = 1e-7)
        @test isapprox(ref_intrinsics.scale_y, calib_intrinsics.scale_y; atol = 1e-7)
        @test isapprox(ref_intrinsics.principal_point, calib_intrinsics.principal_point; atol = 1e-7)

        ref_extrinsics = extrinsics(ref_camᵢ)
        calib_extrinsics = extrinsics(calib_camᵢ)

        ref_origin = MultipleViewGeometry.origin(ref_extrinsics)
        calib_origin = MultipleViewGeometry.origin(calib_extrinsics)
        @test isapprox(ref_origin, calib_origin; atol = 1e-7)

        ref_basis_vectors = MultipleViewGeometry.basis_vectors(ref_extrinsics)
        calib_basis_vectors  = MultipleViewGeometry.basis_vectors(calib_extrinsics)
        @test isapprox(ref_basis_vectors[1], calib_basis_vectors[1]; atol = 1e-7)
        @test isapprox(ref_basis_vectors[2], calib_basis_vectors[2]; atol = 1e-7)
        @test isapprox(ref_basis_vectors[3], calib_basis_vectors[3]; atol = 1e-7)

        @test isapprox(matrix(ref_extrinsics), matrix(calib_extrinsics), atol = 1e-9)
    end
    # Test on noisy data
    Λ = SMatrix{2,2}(I) * 0.005
    Random.seed!(1234)
    𝓞 = [apply_noise(𝓜[i], Λ) for i = 1:length(𝓜)]
    calibrated_cameras = calibrate(world, 𝓞)
    for i = 1:4
        ref_camᵢ = reference_cameras[i]
        calib_camᵢ = calibrated_cameras[i]

        ref_intrinsics = intrinsics(ref_camᵢ)
        calib_intrinsics = intrinsics(calib_camᵢ)
        # Verify intrinsics are close to the ground truth.
        @test isapprox(ref_intrinsics.focal_length, calib_intrinsics.focal_length; atol = 5)
        @test isapprox(ref_intrinsics.skewedness, calib_intrinsics.skewedness; atol = 0.003)
        @test isapprox(ref_intrinsics.scale_x, calib_intrinsics.scale_x; atol = 1e-7)
        @test isapprox(ref_intrinsics.scale_y, calib_intrinsics.scale_y; atol = 0.002)
        @test isapprox(ref_intrinsics.principal_point, calib_intrinsics.principal_point; atol = 0.6)

        ref_extrinsics = extrinsics(ref_camᵢ)
        calib_extrinsics = extrinsics(calib_camᵢ)

        ref_origin = MultipleViewGeometry.origin(ref_extrinsics)
        calib_origin = MultipleViewGeometry.origin(calib_extrinsics)
        @test isapprox(ref_origin, calib_origin; atol = 120)

        ref_basis_vectors = MultipleViewGeometry.basis_vectors(ref_extrinsics)
        calib_basis_vectors  = MultipleViewGeometry.basis_vectors(calib_extrinsics)
        @show isapprox(matrix(ref_extrinsics), matrix(calib_extrinsics); atol = 120)
    end
end

@testset "Sole Camera Calibration Analytic Jacobian" begin
    world = construct_calibration_world()

    pinhole₁ = Pinhole(intrinsics = IntrinsicParameters(width = 640, height = 480, focal_length = 100))
    analogue_image₁ = AnalogueImage(coordinate_system = OpticalSystem())
    camera₁ = Camera(image_type = analogue_image₁, model = pinhole₁)
    𝐑₁ = SMatrix{3,3,Float64,9}(rotxyz(0*(pi/180), 0*(pi/180), 0*(pi/180)))
    𝐭₁ = [-500.0,0.0, -3000.0]
    camera₁ = relocate(camera₁, 𝐑₁, 𝐭₁)

    pinhole₂ = Pinhole(intrinsics = IntrinsicParameters(width = 640, height = 480, focal_length = 100))
    analogue_image₂ = AnalogueImage(coordinate_system = OpticalSystem())
    camera₂ = Camera(image_type = analogue_image₂, model = pinhole₂)
    𝐑₂ = SMatrix{3,3,Float64,9}(rotxyz(0*(pi/180), 10*(pi/180), 10*(pi/180)))
    𝐭₂ = [500.0,0.0, -3000.0]
    camera₂ = relocate(camera₂, 𝐑₂, 𝐭₂)

    pinhole₃ = Pinhole(intrinsics = IntrinsicParameters(width = 640, height = 480, focal_length = 100))
    analogue_image₃ = AnalogueImage(coordinate_system = OpticalSystem())
    camera₃ = Camera(image_type = analogue_image₃, model = pinhole₃)
    𝐑₃ = SMatrix{3,3,Float64,9}(rotxyz(0*(pi/180), -10*(pi/180), 10*(pi/180)))
    𝐭₃ = [0.0,0.0, -2500.0]
    camera₃ = relocate(camera₃, 𝐑₃, 𝐭₃)

    pinhole₄ = Pinhole(intrinsics = IntrinsicParameters(width = 640, height = 480, focal_length = 100))
    analogue_image₄ = AnalogueImage(coordinate_system = OpticalSystem())
    camera₄ = Camera(image_type = analogue_image₄, model = pinhole₄)
    𝐑₄ = SMatrix{3,3,Float64,9}(rotxyz(0*(pi/180), 0*(pi/180), -10*(pi/180)))
    𝐭₄ = [0.0,500.0, -3500.0]
    camera₄ = relocate(camera₄, 𝐑₄, 𝐭₄)

    reference_cameras = [camera₁, camera₂, camera₃, camera₄]

    calibrate = CalibrateCamera()
    # Determine projections of the 3D points in each camera view.
    aquire = AquireImage()
    𝓜 = [aquire(world, camera) for camera in reference_cameras]

    # Test on noiseless data
    calibrated_cameras = calibrate(world, 𝓜)

    # Compare the analytic jacobian for the parameter refinement
    # against the numerical jacobian.
    ℰ =  [matrix(extrinsics(camera)) for camera in calibrated_cameras]
    𝐀 = matrix(intrinsics(calibrated_cameras[1]))
    𝐤 = coefficients(distortion(calibrated_cameras[1]))

    task = MultipleViewGeometry.CameraCalibrationTask()
    # The total number of views.
    M = length(ℰ)
    # Camera intrinsics, lens distortion and extrinsics.
    𝛈 = MultipleViewGeometry.compose_parameter_vector(𝐀, 𝐤, ℰ)
    # The total number of image points across all views.
    N = sum(map(x-> length(x), 𝓜))
    # The projections of the 3D points onto each image, and the actual 3D points.
    @unpack points = world
    observations = Observations(tuple(tuple(𝓜, points)))

    # Initialise the residual vector so that it need not be recreated for
    # each iteration of the LevenbergMarquardt optimization step.
    residuals = zeros(Float64, 2*N*M)
    objective = MultipleViewGeometry.SumOfSquares(task, MultipleViewGeometry.VectorValuedObjective(task, residuals))

    # TODO Instantiate the Jacobian matrix so that it need not be recreated
    # for each step of the optimization loop.
    jacobian_matrix = MultipleViewGeometry.JacobianMatrix(task, objective, observations, zeros(Float64, 2*N, length(𝛈)))
    𝐉 = jacobian_matrix(𝛈)

    @unpack vector_valued_objective = objective
    g = x-> vector_valued_objective(x, observations)
    𝐉₂ = FiniteDiff.finite_difference_jacobian(g, 𝛈)

    @test isapprox(𝐉, 𝐉₂,;atol = 1e-4)
end
