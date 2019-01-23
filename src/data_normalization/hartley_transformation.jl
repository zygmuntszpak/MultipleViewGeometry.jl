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
function hartley_transformation(ℳ::Vector{T})::SMatrix where T <:AbstractArray
    if isempty(ℳ)
        throw(ArgumentError("Array cannot be empty."))
    end
    npts = length(ℳ)
    ndim = length(ℳ[1])
    𝐜 = centroid(ℳ)
    σ = root_mean_square(ℳ, 𝐜)
    σ⁻¹ = 1 / σ
    𝐓 = SMatrix{ndim+1,ndim+1,Float64, (ndim+1)^2}([σ⁻¹*Matrix{Float64}(I,ndim,ndim) -σ⁻¹*𝐜 ; zeros(1,ndim) 1.0])
end

function centroid(positions::Vector{T}) where T <: AbstractArray
    x = zeros(T)
    for pos ∈ positions
        x = x + pos
    end
    return x / length(positions)
end

function root_mean_square(ℳ::Vector{T}, 𝐜::T ) where  T <: AbstractArray
    total = 0.0
    npts = length(ℳ)
    ndim = length(ℳ[1])
    for 𝐦 ∈ ℳ
         total  = total + ∑((𝐦-𝐜).^2)
    end
    σ = √( (1/(npts*ndim)) * total)
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
function hartley_normalization(ℳ::Vector{<:AbstractArray})
    𝒪, 𝐓 = hartley_normalization!(copy(ℳ))
end

function hartley_normalization!(ℳ::Vector{<:AbstractArray})
    𝐓 = hartley_transformation(ℳ)
    map!(ℳ , ℳ) do 𝐦
         hom⁻¹(𝐓 * hom(𝐦))
    end
     ℳ, 𝐓
end
