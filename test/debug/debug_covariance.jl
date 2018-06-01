using MAT
using MultipleViewGeometry, Base.Test
using MultipleViewGeometry.ModuleCostFunction
using MultipleViewGeometry.ModuleTypes
using BenchmarkTools, Compat
using StaticArrays

#"C:\Users\Spock\Desktop\tempdata"

file = matopen("test/debug/F.mat")
F = read(file,"F")
close(file)

file = matopen("test/debug/mPts1List.mat")
mpts1list = read(file,"mPts1List")
close(file)

file = matopen("test/debug/mPts2List.mat")
mpts2list = read(file,"mPts2List")
close(file)

file = matopen("test/debug/covs.mat")
covs = read(file,"covs")
close(file)

file = matopen("test/debug/Ct.mat")
Ct = read(file,"Ct")
close(file)

npts, dim  = size(mpts1list[1])



s²  = 1e-7
s = sqrt(s² )
Λ₁ =  [SMatrix{3,3}(s² *diagm([1.0,1.0,0.0])) for i = 1:npts]
Λ₂ =  [SMatrix{3,3}(s² *diagm([1.0,1.0,0.0])) for i = 1:npts]

ntrial = 10000
ℱ = zeros(9,ntrial)
for itrial = 1:ntrial
  mpts1 = mpts1list[itrial]
  mpts2 = mpts2list[itrial]
  𝒪 = [Point2DH(vcat(mpts1[i,:],[1])) for i =1:npts]
  𝒪ʹ= [Point2DH(vcat(mpts2[i,:],[1])) for i =1:npts]

  𝐅₀ = estimate(FundamentalMatrix(),DirectLinearTransform(), (𝒪, 𝒪ʹ))
  𝐅 = estimate(FundamentalMatrix(),
                           FundamentalNumericalScheme(reshape(𝐅₀,9,1), 5, 1e-10),
                                                             (Λ₁,Λ₂), (𝒪, 𝒪ʹ))
  𝐟 = reshape(𝐅,9,1)
  𝐟 = 𝐟 / norm(𝐟)
  ℱ[:,itrial] = 𝐟
end

norm(F[:,1])
norm(ℱ[:,1])

F[:,1]
ℱ[:,1]


index = 2
mpts1 = mpts1list[index]
mpts2 = mpts2list[index]
𝒪 = [Point2DH(vcat(mpts1[i,:],[1])) for i =1:npts]
𝒪ʹ= [Point2DH(vcat(mpts2[i,:],[1])) for i =1:npts]
T1 = hartley_transformation(𝒪)
T2 = hartley_transformation(𝒪ʹ)

𝐅₀  = estimate(FundamentalMatrix(), DirectLinearTransform(), (𝒪, 𝒪ʹ))
𝐅 = estimate(FundamentalMatrix(),
                        FundamentalNumericalScheme(reshape(𝐅₀,9,1), 5, 1e-10),
                                                          (Λ₁,Λ₂), (𝒪, 𝒪ʹ))
𝐟 = reshape(𝐅,9,1)
𝐟 = 𝐟 / norm(𝐟)

𝐅 = 𝐅 / norm(𝐅)


𝛍 = mean(ℱ,2)
𝛍 = 𝛍 / norm(𝛍)
d = length(𝛍)
𝐏 = eye(d) - norm(𝛍)^-2 * (𝛍*𝛍')
Cₜ = zeros((d,d))
for itrial = 1:ntrial
    𝐟 =  ℱ[:,itrial]
    Cₜ = Cₜ + 𝐏*(𝐟-𝛍) *(𝐟-𝛍)'*𝐏'
end
Cₜ = Cₜ / ntrial

#
# 𝛍 = mean(F,2)
# 𝛍 = 𝛍 / norm(𝛍)
# d = length(𝛍)
# 𝐏 = eye(d) - norm(𝛍)^-2 * (𝛍*𝛍')
# Cₜ = zeros((d,d))
# for itrial = 1:ntrial
#     𝐟 =  F[:,itrial]
#     Cₜ = Cₜ + 𝐏*(𝐟-𝛍) *(𝐟-𝛍)'*𝐏'
# end
# Cₜ = Cₜ / ntrial

# Q  = ℱ;
# for itrial = 1:ntrial
#     q =  Q[:,itrial]
#     Q[:,itrial] = Q[:,itrial] / sign(q[end-1])
# end
#
# 𝛍 = mean(Q,2)
# 𝛍 = 𝛍 / norm(𝛍)
# d = length(𝛍)
# 𝐏 = eye(d) - norm(𝛍)^-2 * (𝛍*𝛍')
# Cₜ = zeros((d,d))
# for itrial = 1:ntrial
#     𝐟 =  Q[:,itrial]
#     Cₜ = Cₜ + 𝐏*(𝐟-𝛍) *(𝐟-𝛍)'*𝐏'
# end
# Cₜ = Cₜ / ntrial

#
# A = Cₜ
#
#
# # Validate the covariance matrix of an estimate based on the AML cost function.
C = covariance_matrix(AML(),FundamentalMatrix(), reshape(𝐅,9,1), (Λ₁,Λ₂), (𝒪 , 𝒪ʹ))
C1 = covariance_matrix(AML(),FundamentalMatrix(),  𝐟, (Λ₁,Λ₂), (𝒪 , 𝒪ʹ))
C2 = covariance_matrix_debug(AML(),FundamentalMatrix(), reshape(𝐅,9,1), (Λ₁,Λ₂), (𝒪 , 𝒪ʹ))


C3 = covariance_matrix_debug(AML(),FundamentalMatrix(), 𝐟, (Λ₁,Λ₂), (𝒪 , 𝒪ʹ))

Cn = covariance_matrix_normalised(AML(),FundamentalMatrix(), 𝐟, (Λ₁,Λ₂), (𝒪 , 𝒪ʹ))

# (diag(C3)./diag(Ct) .-1)*100
#
# norm(diag(C3)-diag(Ct)) / norm(diag(Ct)) * 100
#
norm(C3-Ct) / norm(Ct) * 100
norm(C2-Ct) / norm(Ct) * 100
norm(Cn-Ct) / norm(Ct) * 100
norm(C1-Ct) / norm(Ct) * 100


# Cx = covariance_matrix(AML(),FundamentalMatrix(),  𝐟 , (Λ₁,Λ₂), (𝒪 , 𝒪ʹ))
#
#
# P = (1/norm(𝐟)) * (eye(9) - ((𝐟*𝐟') / norm(𝐟)^2) )
#
# P2 = eye(9) - norm(𝐟)^-2 * (𝐟*𝐟')
#
# cost(AML(),FundamentalMatrix(), 𝐟, (Λ₁,Λ₂), (𝒪 , 𝒪ʹ))


# A = Float64.([3 0 2; 2 0 -2; 0 1 1])
# size(A)
# B = similar(A)
# rows, cols = size(A)
# for r = 1:rows
#     for c =1:cols
#         B[r,c] = det(A)
#     end
# end

Q = reshape(𝐟,(3,3))
