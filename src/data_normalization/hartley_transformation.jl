"""
    hartley_transformation(pts::AbstractArray{T}) where T<:HomogeneousPoint

Returns a matrix which can be used to map a set of ``d``-dimensional Cartesian
points  which are represented by ``\\text{length-}(d+1)`` homogeneous coordinates into a
data-dependent coordinate system. In the data-dependent coordinate system the
origin is the center of mass (centroid) of the points  and the root-mean-square
distance of the points to the origin is equal to ``\\sqrt{d}``.

# Details

A point in ``\\mathbb{R}^d`` with
Cartesian coordinates  ``\\left(m_1, m_2, \\ldots, m_d \\right)`` can also be
expressed in homogeneous coordinates with the vector ``\\mathbf{m} =
\\left[m_1, m_2, \\ldots, m_d , 1 \\right]^\\top``.

Suppose one has a set ``\\left \\{ \\mathbf{m}_n \\right \\}_{n = 1}^{N} `` of
Cartesian points which are represented by homogeneous coordinates.
Let
```math
\\overline{\\mathbf{m}} = \\frac{1}{N} \\sum_{n = 1}^{N} \\mathbf{m}_n
\\quad  \\text{and} \\quad
\\sigma = \\left( \\frac{1}{d \\times n} \\sum_{n = 1}^{N}  \\left \\| \\mathbf{m}_n -
\\overline{\\mathbf{m}} \\right \\|^{2} \\right)^{1/2}
```
represent the centroid of the points and the root-mean-square distance of the
points to the centroid, respectively.

This function returns the matrix
```math
\\mathbf{T} =
\\begin{bmatrix}
\\sigma^{-1} & 0              &   0     & \\ldots        & -\\sigma^{-1} \\overline{m}_1 \\\\
           0 & \\sigma^{-1}   &   0     & \\ldots        & -\\sigma^{-1} \\overline{m}_2 \\\\
           0 & 0              & \\ddots &  0             &  \\vdots                      \\\\
           \\vdots & \\vdots  &   0     &  \\sigma^{-1}  & -\\sigma^{-1} \\overline{m}_d \\\\
           0 & 0              &   0     &           0    &                             1
\\end{bmatrix}
```
such that a transformed point ``\\tilde{\\mathbf{m}}_n = \\mathbf{T} \\mathbf{m}_n``
has a root-mean-square distance to the origin of a new coordinate system
equal to ``\\sqrt{d}``.


"""
function hartley_transformation(pts::AbstractArray{T}) where T<:HomogeneousPoint
    if isempty(pts)
        throw(ArgumentError("Array cannot be empty."))
    end
    # Convert list of homogeneous coordinates into an npts x ndim matrix.
    npts = length(pts);
    ndim = length(pts[1].coords);
    array = reinterpret(Float64,pts,(npts,ndim))
    𝐌 = transpose(reshape(array,(ndim,npts)))
    _hartley_transformation(view(𝐌,:,1:ndim-1))

end

function _hartley_transformation(𝐌::AbstractArray{T,2})::Matrix{T} where T<:Number
    if isempty(𝐌)
        throw(ArgumentError("Array cannot be empty."))
    end
    𝐜 = mean(𝐌,1)
    ndim = length(𝐜)
    # Compute root mean square distance of each point to the centroid.
    σ = √((1/length(𝐌)) * ∑((𝐌 .- 𝐜).^2))
    σ⁻¹ = 1./σ
    𝐓::Matrix{T} = [σ⁻¹*eye(ndim) -σ⁻¹*transpose(𝐜);
                    zeros(1,ndim)                1]

end

"""
    hartley_normalization(pts::AbstractArray{T}) where T<:HomogeneousPoint

Maps a set of ``d``-dimensional Cartesian points  which are represented by
``\\text{length-}(d+1)`` homogeneous coordinates into a data-dependent coordinate system.
In the data-dependent coordinate system the origin is the center of mass
(centroid) of the points  and the root-mean-square distance of the points to the
origin is equal to ``\\sqrt{d}``.

# Details

A point in ``\\mathbb{R}^d`` with
Cartesian coordinates  ``\\left(m_1, m_2, \\ldots, m_d \\right)`` can also be
expressed in homogeneous coordinates with the vector ``\\mathbf{m} =
\\left[m_1, m_2, \\ldots, m_d , 1 \\right]^\\top``.

Suppose one has a set ``\\left \\{ \\mathbf{m}_n \\right \\}_{n = 1}^{N} `` of
Cartesian points which are represented by homogeneous coordinates.
Let
```math
\\overline{\\mathbf{m}} = \\frac{1}{N} \\sum_{n = 1}^{N} \\mathbf{m}_n
\\quad  \\text{and} \\quad
\\sigma = \\left( \\frac{1}{d \\times n} \\sum_{n = 1}^{N}  \\left \\| \\mathbf{m}_n -
\\overline{\\mathbf{m}} \\right \\|^{2} \\right)^{1/2}
```
represent the centroid of the points and the root-mean-square distance of the
points to the centroid, respectively.

This function returns a new set of points ``\\left \\{ \\tilde{\\mathbf{m}}_n \\right \\}_{n = 1}^{N} ``
where ``\\tilde{\\mathbf{m}}_n = \\mathbf{T} \\mathbf{m}_n`` for each ``n``, and
```math
\\mathbf{T} =
\\begin{bmatrix}
\\sigma^{-1} & 0              &   0     & \\ldots        & -\\sigma^{-1} \\overline{m}_1 \\\\
           0 & \\sigma^{-1}   &   0     & \\ldots        & -\\sigma^{-1} \\overline{m}_2 \\\\
           0 & 0              & \\ddots &  0             &  \\vdots                      \\\\
           \\vdots & \\vdots  &   0     &  \\sigma^{-1}  & -\\sigma^{-1} \\overline{m}_d \\\\
           0 & 0              &   0     &           0    &                             1
\\end{bmatrix}.
```
These new points have the property that their root-mean-square distance to
the origin of the coordinate system is equal to ``\\sqrt{d}``.


"""
function hartley_normalization(pts::AbstractArray{T}) where T<:HomogeneousPoint
    𝐓  = hartley_transformation(pts)
    hartley_normalization!(copy(pts))
end

function hartley_normalization!(pts::AbstractArray{T}) where T<:HomogeneousPoint
    𝐓 = hartley_transformation(pts)
    map!(pts , pts) do p
        𝐦 = collect(p.coords)
        𝐦 = 𝑛(𝐓 * 𝐦)
        HomogeneousPoint(tuple(𝐦...))
    end
    (pts,𝐓)
end
