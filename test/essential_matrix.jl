@testset "Essential Matrix Instantiation" begin
    world, camera₁, camera₂ = generate_multi_planar_world()


    𝐄 = matrix(EssentialMatrix(camera₁, camera₂))

    # Result 9.17 of R. Hartley and A. Zisserman, “Two-View Geometry,” Multiple View Geometry in Computer Vision
    # A 3 by 3 matrix is an essential matrix if and only if two of its singular values
    # are equal, and the third is zero.
    U, S, V = svd(𝐄)

    @test isapprox(abs(S[1] - S[2]), 0.0; atol = 1e-14)
    @test isapprox(S[3], 0.0; atol = 1e-10)
end
