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
codomain_coord_rep(::AngularPositionAroundPoint{R2,S1}) = ChartRep()

function task_map_emb_chart(::EmbRep, ::ChartRep, xme, task_map::AngularPositionAroundPoint{R2,S1}, 
        ::Chart{Angleπ,S1})
    a = xme - task_map.position_center
    d = norm(a)
    if d == 0
        return SA[0.]
    else
        â = a ./ norm(a)
        return SA[atan(â[2], â[1])]
    end
end

function task_map_emb_chart(::EmbRep, ::ChartRep, xme, task_map::AngularPositionAroundPoint{R2,S1}, 
        C2π::Chart{Angle2π,S1})
    Cπ = Chart{Angleπ,S1}()
    xnπ = task_map_emb_chart(xme, task_map, Cπ)
    chart_transition(xnπ, Cπ, C2π)
end

# S2
codomain_coord_rep(::AngularPositionAroundPoint{R3,S2}) = EmbRep()

function task_map_emb(::EmbRep, ::EmbRep, xme, task_map::AngularPositionAroundPoint{R3,S2})
    Δx = xme - task_map.position_center
    d = norm(Δx)
    if d == 0
        return @SVector zeros(3)
    else
        return Δx/norm(Δx)
    end
end