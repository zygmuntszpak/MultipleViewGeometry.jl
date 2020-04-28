const  √ = sqrt
const  ∑ = sum
const  ⊗ = kron

function hom⁻¹(v::StaticVector)
    if isapprox(v[end], 0.0; atol = 1e-14)
        pop(v)
    else
        pop(v / v[end])
    end
end

function hom(v::StaticVector)
    push(v,1)
end

function ∂hom⁻¹(𝐧::StaticVector)
    k = length(𝐧)
    𝐞ₖ = push(zeros(SVector{k-1}),1.0)
    𝐈 = SMatrix{3,3}(1.0I)
    1/𝐧[k]*𝐈 - 1/𝐧[k]^2 * 𝐧 * 𝐞ₖ'
end

"""
    𝑛(v::Vector{T}) where T<:Number
Scales a length-``n`` vector ``v``  such that the last component
of the vector is one, provided that the last component is not zero. If the last
component is zero then the vector is left unchanged.
# Details
Suppose the length-``n`` vector ``v`` represents the homogeneous coordinates  of
a point in a projective space. The corresponding Cartesian coordinates  usually
just the first ``n-1`` numbers of homogeneous coordinates divided by  the last
component. So if the last component is one, then the first  ``n-1`` homogeneous
coordinates can be interpreted as Cartesian.  The exceptional case is when the
last component of the homogenenous coordinates is zero. These homogeneous
coordinates are associated with so-called *points at infinity* and have no
Cartesian counterparts.
# Example
```julia
h = [4, 4 , 2]
c = 𝑛(h)
3-element Array{Float64,1}:
 2.0
 2.0
 1.0
```
"""
function 𝑛(v::AbstractArray)
    if v[end] != 0 && v[end] != 1
        v .= v ./ v[end]
    else
        v
    end
end

function 𝑛(v::StaticVector)
    if v[end] != 0 && v[end] != 1
        v / v[end]
    else
        v
    end
end

function ∂𝑛(𝐧::AbstractArray)
    k = length(𝐧)
    𝐞ₖ = fill(0.0,(k,1))
    𝐞ₖ[k] = 1
    1/𝐧[k]*Matrix{Float64}(I, k, k) - 1/𝐧[k]^2 * 𝐧 * 𝐞ₖ'
end


function vec2antisym(v::AbstractArray)
    if length(v) != 3
         throw(ArgumentError("The operation is only defined for a length-3 vector."))
    end
    𝐒  = @SMatrix [   0  -v[3]    v[2] ;
                    v[3]    0    -v[1] ;
                   -v[2]  v[1]      0]

end

function smallest_eigenpair(A::AbstractArray; absolute::Bool = false)
    F = eigen(A)
    index = absolute ? argmin(abs.(F.values)) : argmin(F.values)
    (F.values[index], F.vectors[:,index])
end

function smallest_eigenpair(A::AbstractArray, B::AbstractArray; absolute::Bool = false)
    F = eigen(A, B)
    index = absolute ? argmin(abs.(F.values)) : argmin(F.values)
    (F.values[index], F.vectors[:,index])
end

function smallest_eigenpairs(A::AbstractArray; pairs::Integer = 2, absolute::Bool = false)
    F = eigen(A)
    indices = absolute ? sortperm(abs.(F.values)) : sortperm(F.values)
    (F.values[indices[Base.OneTo(pairs)]], F.vectors[:, indices[Base.OneTo(pairs)]])
end


function minor(A, i, j)
    m, n = size(A)
    B = similar(A, m-1, n-1)
    for j′=1:j-1, i′=1:i-1; B[i′,j′] = A[i′,j′]; end
    for j′=1:j-1, i′=i+1:m; B[i′-1,j′]= A[i′,j′]; end
    for j′=j+1:n, i′=1:i-1; B[i′,j′-1] = A[i′,j′]; end
    for j′=j+1:n, i′=i+1:m; B[i′-1,j′-1] = A[i′,j′]; end
    return B
end

function evaluate_difference(𝛉₁::AbstractVector, 𝛉₂::AbstractVector)
    val₁ = norm( 𝛉₁ / norm(𝛉₁) - 𝛉₂/norm(𝛉₂))
    val₂ = norm( 𝛉₁ / norm(𝛉₁) + 𝛉₂/norm(𝛉₂))
    return min(val₁, val₂)
end

function sigmoid(x::Number, slope::Number = 0.5)
    # NB Slope must be a positive number
    return 1 / (1 + exp(-slope*x))
end
