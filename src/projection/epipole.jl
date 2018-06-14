function epipole(𝐅::AbstractArray)
    𝐔, 𝐒, 𝐕 = svd(𝐅)
    𝐞 = 𝑛(MVector{3}(𝐕[:,end]))
    Point2DH(𝐞)
end
