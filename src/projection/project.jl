function project(e::Pinhole, 𝐏::AbstractArray{T1,2}, 𝒳::AbstractArray{T2}) where {T1<:Real,T2<:HomogeneousPoint}

    if size(𝐏) != (3,4)
        throw(ArgumentError("Expect 3 x 4 projection matrix."))
    end
    ℳ = map(𝒳) do X
        𝐗 = collect(X.coords)
        𝐦 = 𝑛(𝐏 * 𝐗)
        HomogeneousPoint(tuple(𝐦...))
    end
end
