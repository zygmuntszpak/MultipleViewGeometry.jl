using MultipleViewGeometry, Base.Test

# Vectors are scaled so that the last component is unity.
v = [2,2,2]
@test 𝑛(v) == [1.0, 1.0, 1.0]
# Vectors which represent points at infinity are unchanged.
v = [2,2,0]
@test 𝑛(v) == [2,2,0]
