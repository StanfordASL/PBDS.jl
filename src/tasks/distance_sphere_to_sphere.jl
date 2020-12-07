struct DistanceSphereToSphere{M,R1,S,m} <: TaskMap{M,R1,S}
    center::SVector{m,S}
    radius1::S
    radius2::S
end
DistanceSphereToSphere{M,R1}(c::SVector{m,S}, r1::S, r2::S) where {M,R1,S,m} =
    DistanceSphereToSphere{M,R1,S,dim(M)}(c, r1, r2)

struct DistanceSphereToSphereT{M,R1,S,m} <: TaskMapT{M,R1,S}
    base_map::DistanceSphereToSphere{M,R1,S,m}
end
DistanceSphereToSphereT{M,R1}(c::SVector{m,S}, r1::S, r2::S) where {M,R1,S,m} =
    DistanceSphereToSphereT{M,R1,S,m}(DistanceSphereToSphere{M,R1}(c, r1, r2))
    
task_map_emb(::EmbRep, ::EmbRep, xme, task_map::DistanceSphereToSphere{ℝ{m},R1}) where m =
    SA[norm(xme - task_map.center) - task_map.radius1 - task_map.radius2]