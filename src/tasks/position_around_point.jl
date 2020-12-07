struct PositionAroundPoint{M,N,S,m} <: TaskMap{M,N,S}
    position_center::SVector{m,S}
end

PositionAroundPoint{M,N}(c::SVector{m,S}) where {M,N,S,m} = PositionAroundPoint{M,N,S,dim(M)}(c)

struct PositionAroundPointT{M,N,S,m} <: TaskMapT{M,N,S}
    base_map::PositionAroundPoint{M,N,S,m}
end
PositionAroundPointT{M,N}(c::SVector{m,S}) where {M,N,S,m} =
    PositionAroundPointT{M,N,S,m}(PositionAroundPoint{M,N}(c))

task_map_emb(::EmbRep, ::EmbRep, xm, task_map::PositionAroundPoint{ℝ{m},ℝ{m},S}) where {m,S} =
    xm - task_map.position_center