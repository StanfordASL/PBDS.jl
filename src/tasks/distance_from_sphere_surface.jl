struct DistanceFromSphereSurface{M,R1,S,m} <: TaskMap{M,R1,S}
    center::SVector{m,S}
    radius::S
end
DistanceFromSphereSurface{M,R1}(c::SVector{m,S}, r::S) where {M,R1,S,m} =
    DistanceFromSphereSurface{M,R1,S,embdim(M)}(c, r)

struct DistanceFromSphereSurfaceT{M,R1,S,m} <: TaskMapT{M,R1,S}
    base_map::DistanceFromSphereSurface{M,R1,S,m}
end
DistanceFromSphereSurfaceT{M,R1}(c::SVector{m,S}, r::S) where {M,R1,S,m} =
    DistanceFromSphereSurfaceT{M,R1,S,m}(DistanceFromSphereSurface{M,R1}(c, r))
    
task_map_emb(::EmbRep, ::EmbRep, xme, task_map::DistanceFromSphereSurface{ℝ{m},R1}) where m =
    SA[abs(norm(xme - task_map.center) - task_map.radius)]