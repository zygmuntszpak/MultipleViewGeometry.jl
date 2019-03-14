# Assume homogeneous coordinates
function perturb(noise::GaussianNoise, σ::Real, 𝒟::Tuple{AbstractArray, Vararg{AbstractArray}})
    𝓔 = deepcopy(𝒟)
    S = length(𝓔)
    for s = 1:S
        ℳ = 𝓔[s]
        N = length(ℳ)
        for n = 1:N
            𝐦 = ℳ[n]
            D = length(𝐦)
            ℳ[n] = 𝐦 + σ*SVector(randn((D,1))...)
        end
    end
    𝓔
end
