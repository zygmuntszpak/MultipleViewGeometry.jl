using MultipleViewGeometry, Base.Test


# Tests for a set of two-dimensional Cartesian points represented by homogeneous
# coordinates.
ℳ = map(HomogeneousPoint,
        [(-10.0, -10.0, 1.0),
         (-10.0,  10.0, 1.0),
         ( 10.0, -10.0, 1.0),
         ( 10.0,  10.0, 1.0)])

ℳʹ, 𝐓 = transform(HomogeneousCoordinates(),CanonicalToHartley(),ℳ)
@test ℳʹ == map(HomogeneousPoint,
                                        [(-1.0,-1.0, 1.0),
                                         (-1.0, 1.0, 1.0),
                                         (1.0, -1.0, 1.0),
                                         (1.0,  1.0, 1.0)])

Λ = [eye(2) for i = 1:length(ℳ)]
Λʹ = transform(CovarianceMatrices(), CanonicalToHartley(), Λ , tuple(𝐓))
for 𝚲 ∈ Λ
    @test 𝚲 == eye(2)
end

for 𝚲 ∈ Λʹ
    @test 𝚲 ≈ eye(2)/100
end

Λ = [eye(4) for i = 1:length(ℳ)]
Λʹ = transform(CovarianceMatrices(), CanonicalToHartley(), Λ , tuple(𝐓,𝐓))
for 𝚲 ∈ Λ
    @test 𝚲 == eye(4)
end

for 𝚲 ∈ Λʹ
    @test 𝚲 ≈ eye(4)/100
end

# Tests for a set of three-dimensional Cartesian points represented by homogeneous
# coordinates.
ℳ = map(HomogeneousPoint,
           [(-10.0, -10.0, -10.0, 1.0),
            (-10.0, -10.0,  10.0, 1.0),
            (-10.0,  10.0, -10.0, 1.0),
            (-10.0,  10.0,  10.0, 1.0),
            ( 10.0, -10.0, -10.0, 1.0),
            ( 10.0, -10.0,  10.0, 1.0),
            ( 10.0,  10.0, -10.0, 1.0),
            ( 10.0,  10.0,  10.0, 1.0)])

ℳʹ, 𝐓 = transform(HomogeneousCoordinates(),CanonicalToHartley(),ℳ)
@test ℳʹ == map(HomogeneousPoint,
                                         [(-1.0,-1.0, -1.0, 1.0),
                                          (-1.0,-1.0,  1.0, 1.0),
                                          (-1.0, 1.0, -1.0, 1.0),
                                          (-1.0, 1.0,  1.0, 1.0),
                                          (1.0, -1.0, -1.0, 1.0),
                                          (1.0, -1.0,  1.0, 1.0),
                                          (1.0,  1.0, -1.0, 1.0),
                                          (1.0,  1.0,  1.0, 1.0)])
@test 𝐓 == [0.1 0.0 0.0 -0.0;
            0.0 0.1 0.0 -0.0;
            0.0 0.0 0.1 -0.0;
            0.0 0.0 0.0 1.0]

Λ = [eye(3) for i = 1:length(ℳ)]
Λʹ = transform(CovarianceMatrices(), CanonicalToHartley(), Λ , tuple(𝐓))
for 𝚲 ∈ Λ
    @test 𝚲 == eye(3)
end
for 𝚲 ∈ Λʹ
    @test 𝚲 ≈ eye(3)/100
end
