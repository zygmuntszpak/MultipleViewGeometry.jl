using MAT
using MultipleViewGeometry, Test
using MultipleViewGeometry.ModuleCostFunction
using MultipleViewGeometry.ModuleTypes
using BenchmarkTools
using StaticArrays, LinearAlgebra, SparseArrays

file = matopen("debug/mPts1.mat")
mpts1 = read(file,"mPts1")
close(file)

file = matopen("debug/mPts2.mat")
mpts2 = read(file,"mPts2")
close(file)

file = matopen("debug/covs.mat")
covs = read(file,"covs")
close(file)

file = matopen("debug/Ct.mat")
Ct = read(file,"Ct")
close(file)

npts, dim  = size(mpts1)

s²  = 1e-7
s = sqrt(s² )
Λ₁ =  [SMatrix{3,3}(s² * Matrix(Diagonal([1.0,1.0,0.0]))) for i = 1:npts]
Λ₂ =  [SMatrix{3,3}(s² *Matrix(Diagonal([1.0,1.0,0.0]))) for i = 1:npts]

𝒪 = [Point2DH(vcat(mpts1[i,:],[1])) for i =1:npts]
𝒪ʹ= [Point2DH(vcat(mpts2[i,:],[1])) for i =1:npts]

𝐅₀ = estimate(FundamentalMatrix(),DirectLinearTransform(), (𝒪, 𝒪ʹ))
𝐅 = estimate(FundamentalMatrix(),
                           FundamentalNumericalScheme(reshape(𝐅₀,9,1), 5, 1e-10),
                                                             (Λ₁,Λ₂), (𝒪, 𝒪ʹ))
𝐟 = reshape(𝐅,9,1)
𝐟 = 𝐟 / norm(𝐟)

# # Validate the covariance matrix of an estimate based on the AML cost function.
C1 = covariance_matrix(AML(), HessianApproximation(), FundamentalMatrix(), reshape(𝐅,9,1), (Λ₁,Λ₂), (𝒪 , 𝒪ʹ))
C2 = covariance_matrix(AML(), CanonicalApproximation(), FundamentalMatrix(),  𝐟, (Λ₁,Λ₂), (𝒪 , 𝒪ʹ))

@test norm(C1-Ct) / norm(Ct) * 100 < 0.6
@test norm(C2-Ct) / norm(Ct) * 100 < 0.6
