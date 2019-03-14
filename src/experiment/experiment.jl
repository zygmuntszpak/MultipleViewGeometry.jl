# TODO split into smaller functions
function conduct_experiment(xp::Experiment, cost_function::CostFunction)
    task = xp.task
    P = length(xp.participants)
    xp_result = ExperimentResult(xp.description, xp.task, Vector{ParticipantResult}(undef, P), xp.cost_function)
    for (p, participant) in pairs(xp.participants)
        algorithm = participant.algorithm
        @show typeof(algorithm)

        results = ParticipantResult(participant, Dict{String, ConditionResult}(), 0.0)
        for (str, condition) in participant.conditions
            @show condition.description
            N = length(condition.trials)
            trial_results = Vector{TrialResult}(undef,N)
            for (t, trial) in pairs(condition.trials)
                    𝚹 = estimate(task, algorithm, trial.perturbed_training_data)
                    #residual = assess(cost_function, task, 𝚹, trial.pure_training_data)
                    data = trial.perturbed_training_data
                    #@show length(𝚹)
                    residual = assess(cost_function, task, 𝚹, data)
                    Base.display(residual)
                    trial_results[t] = TrialResult(cost_function, length(first(data)), residual)
            end
            results.conditions[str] = ConditionResult(condition.description, trial_results)
        end
        xp_result.results[p] = results
    end
    xp_result
end

function initialize_camera_pair()
    # Fix random seed.
    Random.seed!(1234)

    # Construct two camera matrices
    f = 50
    image_width = 640
    image_height = 480
    𝐊₁ = @SMatrix [f 0 0 ;
                   0 f 0 ;
                   0 0 1 ]
    𝐑₁ = SMatrix{3,3,Float64,9}(rotxyz(0, 1*(pi/180), 0))
    𝐭₁ = [-300.0, 0.0, -50.0]

    𝐊₂ = @SMatrix [f 0 0 ;
                   0 f 0 ;
                   0 0 1 ]

    𝐑₂ = SMatrix{3,3,Float64,9}(rotxyz(0, -1*(pi/180), 0))
    𝐭₂ = [300.0, 0.0, 5.0]


    world_basis = (Vec(1.0, 0.0, 0.0), Vec(0.0, 1.0, 0.0), Vec(0.0, 0.0, 1.0))
    camera_basis = (Point(0.0, 0.0, 0.0), Vec(-1.0, 0.0, 0.0), Vec(0.0, -1.0, 0.0), Vec(0.0, 0.0, 1.0))
    picture_basis = (Point(0.0, 0.0), Vec(-1.0, 0.0), Vec(0.0, -1.0))

    camera₁ = Pinhole(image_width, image_height, f, camera_basis..., picture_basis...)
    camera₂ = Pinhole(image_width, image_height, f, camera_basis..., picture_basis...)

    relocate!(camera₁, 𝐑₁, 𝐭₁)
    relocate!(camera₂, 𝐑₂, 𝐭₂)

    𝐑₁′, 𝐭₁′ = ascertain_pose(camera₁, world_basis... )
    𝐊₁′ = obtain_intrinsics(camera₁, CartesianSystem())

    𝐑₂′, 𝐭₂′ = ascertain_pose(camera₂, world_basis... )
    𝐊₂′ = obtain_intrinsics(camera₂, CartesianSystem())

    return camera₁, 𝐊₁′, 𝐑₁′, 𝐭₁′ , camera₂, 𝐊₂′, 𝐑₂′,  𝐭₂′
end

function construct_experiment(planar_scene::PlanarScene, trial_count::Int, noise::AbstractRange)

    camera₁, 𝐊₁′, 𝐑₁′,  𝐭₁′, camera₂, 𝐊₂′, 𝐑₂′, 𝐭₂′ = initialize_camera_pair()
    𝐏₁ = construct(ProjectionMatrix(), 𝐊₁′, 𝐑₁′, 𝐭₁′)
    𝐏₂ = construct(ProjectionMatrix(), 𝐊₂′, 𝐑₂′, 𝐭₂′)
    K = planar_scene.plane_count
    𝓜 = Vector{Array{Point{2,Float64},1}}(undef,K)
    𝓜ʹ = Vector{Array{Point{2,Float64},1}}(undef,K)
    𝓞 = Vector{Array{Point{2,Float64},1}}(undef,K)
    𝓞ʹ = Vector{Array{Point{2,Float64},1}}(undef,K)
    𝓗 =  Vector{SArray{Tuple{3,3},Float64,2,9}}(undef,K)
    for k = 1:K
        N = 5000
        # Normals and (negative) distance from origin
        𝐧 = [0.0, 0.0, 1.0]
        d = -100.0 * k

        z = -d
        x_range = -1000.0:1000.0
        y_range = -1000.0:1000.0

        𝒳 = generate_planar_points(-1000.0:1000.0, -1000.0:1000.0, z , N)
        𝐇 = construct(HomographyMatrix(),𝐊₁′,𝐑₁′,𝐭₁′,𝐊₂′,𝐑₂′,𝐭₂′, 𝐧, d)
        𝓗[k] = 𝐇

        # Set of corresponding points.
        ℳ = project(camera₁,𝐏₁,𝒳)
        ℳʹ= project(camera₂,𝐏₂,𝒳)

        # Discard corresponding points which fall outside specific rectangular regions in the first image.
        𝒪, 𝒪ʹ = crop(HyperRectangle(Vec(0,0),Vec(200,200)), (ℳ, ℳʹ))

        𝓜[k] = ℳ
        𝓜ʹ[k] = ℳʹ

        𝓞[k] = 𝒪
        𝓞ʹ[k] = 𝒪ʹ
    end

    𝓟 = Vector{Array{Point{2,Float64},1}}(undef,K)
    𝓟ʹ = Vector{Array{Point{2,Float64},1}}(undef,K)
    conditions = Dict{String, ExperimentCondition}()
    for σ in noise
        trials = Vector{ExperimentTrial}(undef,trial_count)
        for t = 1:trial_count
            for k = 1:K
                𝒫, 𝒫ʹ = perturb(GaussianNoise(), σ, (𝓞[k], 𝓞ʹ[k]))
                𝓟[k] = 𝒫
                𝓟ʹ[k] =  𝒫ʹ
            end
            trial = ExperimentTrial( (𝓞, 𝓞ʹ),  (𝓟, 𝓟ʹ), (𝓜, 𝓜ʹ),  𝓗)
            trials[t] = trial
        end
        condition = ExperimentCondition("σ = $σ", trials)
        conditions["σ = $σ"] =  condition
    end

    participant₁ = Participant(DirectLinearTransform(), conditions)
    participant₂ = Participant(BundleAdjustment(DirectLinearTransform(), 5, 1e-10), conditions)
    experiment = Experiment("Increasing Noise Level", HomographyMatrix(), [participant₁, participant₂], ReprojectionError())
end
