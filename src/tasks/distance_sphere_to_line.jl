struct DistanceSphereToLine{M,R1,S,m} <: TaskMap{M,R1,S}
    radius::S
    point1::SVector{m,S}
    point2::SVector{m,S}
end
DistanceSphereToLine{M,R1}(r::S, p1, p2) where {M,R1,S} =
    DistanceSphereToLine{M,R1,S,dim(M)}(r, p1, p2)

struct DistanceSphereToLineT{M,R1,S,m} <: TaskMapT{M,R1,S}
    base_map::DistanceSphereToLine{M,R1,S,m}
end

DistanceSphereToLineT{M,R1}(r::S, p1, p2) where {M,R1,S} =
    DistanceSphereToLineT{M,R1,S,m}(DistanceSphereToLine{M,R1}(r, p1, p2))

function task_map_emb(::EmbRep, ::EmbRep, xme, task_map::DistanceSphereToLine{ℝ{m},R1,S}) where {m,S}
    point1, point2 = task_map.point1, task_map.point2
    p1p2 = point2 - point1
    p1x = xme - point1
    p2x = xme - point2
    if dot(p1x, p1p2) > 0.
        if dot(p2x, p1p2) < 0.
            normal_vec = p1x - p1p2*dot(p1x, p1p2)/norm(p1p2)^2
            return SA[norm(normal_vec) - task_map.radius]
        else
            return SA[norm(p2x) - task_map.radius]
        end
    else
        return SA[norm(p1x) - task_map.radius]
    end
end