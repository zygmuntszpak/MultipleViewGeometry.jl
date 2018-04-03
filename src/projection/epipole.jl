function epipole(𝐅::Matrix{T}) where T<:Real
𝐔, 𝐒, 𝐕 = svd(𝐅)
𝐞 = 𝑛(𝐕[:,end])
HomogeneousPoint(tuple(𝐕[:,end]...))
end
