function compose(lv::LatentVariables, 𝛈::AbstractArray)
        N = div(length(𝛈) - 12,  4)
        𝐚  = @view 𝛈[1:9]
        𝐀 = reshape(𝐚, (3,3))
        𝐛 = @view 𝛈[10:12]
        𝐰 = @view 𝛈[end-(N-1):end]
        r = range(13, step = 3, length = N+1)
        𝐯 = reshape(view(𝛈,first(r):last(r)-1), (3,N))
        ℋ = Array{Array{Float64,2},1}(undef,(N,))
        for n = 1:N
            ℋ[n] = 𝐰[n]*𝐀 + 𝐛*𝐯[:,n]'
        end
        ℋ
end
