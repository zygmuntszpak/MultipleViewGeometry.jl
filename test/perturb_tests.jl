using MultipleViewGeometry, Base.Test
using MultipleViewGeometry.ModuleCostFunction
using MultipleViewGeometry.ModuleTypes
using MultipleViewGeometry.ModuleConstraints
using MultipleViewGeometry.ModuleConstruct
using MultipleViewGeometry.ModuleNoise
using BenchmarkTools, Compat
using StaticArrays

# Fix random seed.
srand(1234)

𝒳 = [Point3DH(x,y,z,1.0)
                        for x=-1:1:10 for y=-1:1:10 for z=-1:1:10]
𝒟 = perturb(GaussianNoise(), 1.0, tuple(𝒳) )
𝒳ʹ = 𝒟[1]

N = length(𝒳ʹ)
for n = 1:N
    @test !isapprox(sum(abs.(𝒳[1][1:3]-𝒳ʹ[1][1:3])/4), 0.0; atol = 1e-12)
    # No noise should have been added to the last coordinate.
    @test isapprox(sum(abs.(𝒳[1][4]-𝒳ʹ[1][4])), 0.0; atol = 1e-12)
end
