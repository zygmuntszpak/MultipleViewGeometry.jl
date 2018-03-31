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
function 𝑛(v::Vector{T}) where T<:Real
    if v[end] != 0 && v[end] != 1
        v = v ./ v[end]
    else
        v
    end
end

function smallest_eigenpair(A::AbstractArray{T,2}) where T <:Real
    F = eigfact(A)::Base.LinAlg.Eigen{Float64,Float64,Array{Float64,2},Array{Float64,1}};
    (F[:values][1], F[:vectors][:,1])
end

function vec2antisym(v::AbstractVector{T})  where T<:Real
    if length(v) != 3
         throw(ArgumentError("The operation is only defined for a length-3 vector."))
    end
    𝐒::Matrix{T} = [  0  -v[3]    v[2] ;
                    v[3]    0    -v[1] ;
                   -v[2]  v[1]      0]
end
