using MultipleViewGeometry, Base.Test

𝐊 = eye(3)
𝐑 = eye(3)
𝐭 = [1.0, 1.0, 1.0]

@test construct(ProjectionMatrix(),𝐊,𝐑,𝐭) == [eye(3) -ones(3)]
