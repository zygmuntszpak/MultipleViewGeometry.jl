
@testset "Fundamental Matrix Instantiation" begin
    world, camera₁, camera₂ = generate_multi_planar_world()

    @unpack points, planes = world
    plane = last(planes)
    𝐧 = normal(plane)
    d = distance(plane)

    aquire = AquireImage()
    ℳ = aquire(world, camera₁)
    ℳ′ = aquire(world, camera₂)

    𝐅₁ = matrix(FundamentalMatrix(camera₁, camera₂))
    𝐅₁ = 𝐅₁ / norm(𝐅₁)
    𝐅₁ = 𝐅₁ / sign(𝐅₁[3,3])

    𝐅₂ = matrix(FundamentalMatrix(Projection(camera₁), Projection(camera₂)))
    𝐅₂ = 𝐅₂ / norm(𝐅₂)
    𝐅₂ = 𝐅₂ / sign(𝐅₂[3,3])

    # Verify that the epipolar constraint is satsfied.
    r = [hom(ℳ′[i])' * 𝐅₁  * hom(ℳ[i]) for i = 1:length(ℳ)]
    @test all(isapprox.(r, 0.0; atol = 1e-12))

    r = [hom(ℳ′[i])' * 𝐅₂  * hom(ℳ[i]) for i = 1:length(ℳ)]
    @test all(isapprox.(r, 0.0; atol = 1e-12))
end

@testset "Fundamental Matrix Estimation" begin
    world, camera₁, camera₂ = generate_multi_planar_world()

    aquire = AquireImage()
    ℳ = aquire(world, camera₁)
    ℳ′ = aquire(world, camera₂)

    # Estimate a fundamental matrix between corresponding points
    fundamental_matrix = fit_fundamental_matrix(ℳ , ℳ′, DirectLinearTransform())
    𝐅 = matrix(fundamental_matrix)

    # Verify that the epipolar constraint is satsfied.
    r = [hom(ℳ′[i])' * 𝐅  * hom(ℳ[i]) for i = 1:length(ℳ)]
    @test all(isapprox.(r, 0.0; atol = 1e-12))

end
