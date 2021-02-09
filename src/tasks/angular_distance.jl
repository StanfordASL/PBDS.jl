struct AngularDistance{M,R1,S} <: TaskMap{M,R1,S}
    reference_angle::SVector{1,S} # M = S1
end
AngularDistance{M,R1}(c::SVector{n,S}) where {M,R1,S,n} = AngularDistance{M,R1,S}(c)

struct AngularDistanceT{M,R1,S} <: TaskMapT{M,R1,S}
    base_map::fulltype(AngularDistance{M,R1,S})
end
AngularDistanceT{M,R1}(c::SVector{n,S}) where {M,R1,S,n} =
    AngularDistanceT{M,R1,S}(AngularDistance{M,R1,S}(c))

# S1
domain_coord_rep(::AngularDistance{S1,R1}) = ChartRep()

function task_map_chart_emb(::ChartRep, ::EmbRep, xm, task_map::AngularDistance{S1,R1},
        ::Chart{<:AngleChart,S1})
    Δθ = xm - task_map.reference_angle
    wrapToPi.(Δθ)
end