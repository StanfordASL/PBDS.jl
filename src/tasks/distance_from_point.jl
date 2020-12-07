struct DistanceFromPoint{M,R1,S,m} <: TaskMap{M,R1,S}
    position_center::SVector{m,S}
end
DistanceFromPoint{M,R1}(c::SVector{m,S}) where {M,R1,S,m} = DistanceFromPoint{M,R1,S,embdim(M)}(c)

struct DistanceFromPointT{M,R1,S,m} <: TaskMapT{M,R1,S}
    base_map::DistanceFromPoint{M,R1,S,m}
end
DistanceFromPointT{M,R1}(c::SVector{m,S}) where {M,R1,S,m} =
    DistanceFromPointT{M,R1,S,m}(DistanceFromPoint{M,R1}(c))

task_map_emb(::EmbRep, ::EmbRep, xme, task_map::DistanceFromPoint{ℝ{m},R1}) where m =
    SA[norm(xme - task_map.position_center)]

# Updateable version
struct DistanceFromPointDynamic{M,R1,S,m} <: TaskMap{M,R1,S}
    position_center::Vector{SVector{m,S}}
end
DistanceFromPointDynamic{M,R1}(c::Vector{SVector{m,S}}) where {M,R1,S,m} =
    DistanceFromPointDynamic{M,R1,S,embdim(M)}(c)

struct DistanceFromPointDynamicT{M,R1,S,m} <: TaskMapT{M,R1,S}
    base_map::DistanceFromPointDynamic{M,R1,S,m}
end
DistanceFromPointDynamicT{M,R1}(c::Vector{SVector{m,S}}) where {M,R1,S,m} =
    DistanceFromPointDynamicT{M,R1,S,m}(DistanceFromPointDynamic{M,R1}(c))

task_map_emb(::EmbRep, ::EmbRep, xme, task_map::DistanceFromPointDynamic{ℝ{m},R1}) where m =
    SA[norm(xme - task_map.position_center[1])]