@testset "Planar World Instantiation" begin
    world = PrimitiveWorld()
    @inferred PrimitiveWorld()
    @unpack points, planes = world
    plane = first(planes)
    
    𝐧 = normal(plane)
    d = distance(plane)

    # Verify that the default points lie on the plane
    𝛑 = push(𝐧, -d) # 𝛑 =[n -d]
    for 𝐱 in points
        @test isapprox(dot(𝛑, hom(𝐱)), 0.0; atol = 1e-14)
    end
end
