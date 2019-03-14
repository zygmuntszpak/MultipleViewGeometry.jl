Base.convert(::Type{String}, x::Type{BundleAdjustment{T}}) where T <: EstimationAlgorithm = "BA" * " (" * convert(String,T)  * ")"

Base.convert(::Type{String}, x::Type{DirectLinearTransform}) = "DLT"

Base.convert(::Type{String}, x::Type{FundamentalNumericalScheme}) = "FNS"


get_outcomes(x::AbstractConditionResult) = x.outcomes
get_description(x::AbstractConditionResult) = x.description
get_results(x::AbstractExperimentResult) = x.results
get_conditions(x::AbstractParticipantResult) = x.conditions
get_participant(x::AbstractParticipantResult) = x.participant
get_algorithm(x::AbstractParticipant) = x.algorithm
get_pure_testing_data(x::AbstractExperimentTrial) = x.pure_testing_data
get_pure_training_data(x::AbstractExperimentTrial) = x.pure_training_data
get_correct_parameters(x::AbstractExperimentTrial) = x.correct_parameters

# TODO Revise this with getter methods instead of indexing fields directly.
function tabulate(data::ExperimentResult)
    results = get_results(data)
    P = length(results)
    # Set primary key by enumerating participants.
    t = table(Base.OneTo(P); names = [:id], pkey=:id)
    # Add column with a short description of each algorithm.
    algorithms = Vector{String}(undef,P)
    for p = 1:P
        participant = get_participant(results[p])
        algorithms[p] = convert(String, typeof(get_algorithm(participant)))
    end
    t = pushcol(t, :algorithm => algorithms)
    # Add column for the results of each condition.

    for c in keys(get_conditions(first(results)))
        scores = Vector{Float64}(undef,P)
        for p = 1:P
            conditions = get_conditions(results[p])
            outcomes = get_outcomes(conditions[c])
            T = length(outcomes)
            trial = first(outcomes)
            total = zeros(size(trial.val))
            N = trial.total_points
            for t = 1:T
                total .= total .+ outcomes[t].val
            end
            # (average) root-mean-square error
            scores[p] = mean(sqrt.(total ./ (N*T*4)))
        end
        t = pushcol(t, Symbol(c) => scores)
    end
    t
end
