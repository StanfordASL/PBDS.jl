struct AngularPositionAroundPoint{M,N,S,m} <: TaskMap{M,N,S}
    position_center::SVector{m,S}
end
AngularPositionAroundPoint{M,N}(c::SVector{m,S}) where {M,N,S,m} =
    AngularPositionAroundPoint{M,N,S,dim(M)}(c)

struct AngularPositionAroundPointT{M,N,S,m} <: TaskMapT{M,N,S}
    base_map::AngularPositionAroundPoint{M,N,S,m}
end
AngularPositionAroundPointT{M,N}(c::SVector{m,S}) where {M,N,S,m} =
    AngularPositionAroundPointT{M,N,S,m}(AngularPositionAroundPoint{M,N}(c))

# S1
codomain_coord_rep(::AngularPositionAroundPoint{R2,S1,S}) where S = ChartRep()

function task_map_emb_chart(::EmbRep, ::ChartRep, xme,
        task_map::AngularPositionAroundPoint{R2,S1,S}, ::Chart{Angleπ,S1}) where S
    a = xme - task_map.position_center
    â = a ./ norm(a)
    SA[atan(â[2], â[1])]
end

function task_map_emb_chart(::EmbRep, ::ChartRep, xme,
        task_map::AngularPositionAroundPoint{R2,S1,S}, C2π::Chart{Angle2π,S1}) where S
    Cπ = Chart{Angleπ,S1}()
    xnπ = task_map_emb_chart(xme, task_map, Cπ)
    chart_transition(xnπ, Cπ, C2π)
end

# S2
function task_map_emb_chart(::EmbRep, ::ChartRep, xme,
        task_map::AngularPositionAroundPoint{R3,S2,S}, CN::Chart{<:SterProj,S2}) where S
    a = xme - task_map.position_center
    â = a ./ norm(a)
    emb_to_chart(â, CN)
end

# S2 embedded
function task_map_emb(::EmbRep, ::EmbRep, xme,
        task_map::AngularPositionAroundPoint{R3,S2,S}) where S
    a = xme - task_map.position_center
    â = a ./ norm(a)
end