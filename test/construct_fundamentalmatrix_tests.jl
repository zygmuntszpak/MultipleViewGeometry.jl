using MultipleViewGeometry, Base.Test

𝐊₁ = eye(3)
𝐑₁ = eye(3)
𝐭₁ = [1.0, 1.0, 1.0]
𝐊₂ = eye(3)
𝐑₂ = eye(3)
𝐭₂ = [2.0, 2.0, 2.0]

@test construct(FundamentalMatrix(),𝐊₁,𝐑₁,𝐭₁,𝐊₂,𝐑₂,𝐭₂) == vec2antisym([-1,-1,-1])
