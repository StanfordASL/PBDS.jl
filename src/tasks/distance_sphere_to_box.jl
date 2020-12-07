struct DistanceSphereToBox{M,R1,S,m} <: TaskMap{M,R1,S}
    radius::S
    box::Rect{m,S}
    rotation::RotMatrix{m,S}
end
DistanceSphereToBox{M,R1}(r, b::Rect{m,S}, R) where {M,R1,S,m} =
    DistanceSphereToBox{M,R1,S,dim(M)}(r, b, R)

struct DistanceSphereToBoxT{M,R1,S,m} <: TaskMapT{M,R1,S}
    base_map::DistanceSphereToBox{M,R1,S,m}
end
DistanceSphereToBoxT{M,R1}(r, b::Rect{m,S}, R) where {M,R1,S,m} =
    DistanceSphereToBoxT{M,R1,S,m}(DistanceSphereToBox{M,R1}(r, b, R))

function min_dist_dim(rect::Rect{N,T1}, p::Vec{N,T2}, dim::Int) where {N,T1,T2}
    T = promote_type(T1,T2)
    return max(zero(T), max(minimum(rect)[dim] - p[dim], p[dim] - maximum(rect)[dim]))
end

function min_euclideansq(rect::Rect{N,T1}, p::Union{Rect{N,T2},Vec{N,T2}}) where {N,T1,T2}
    T = promote_type(T1,T2)
    minimum_dist = T(0.0)
    for dim in 1:length(p)
        d = min_dist_dim(rect, p, dim)
        minimum_dist += d*d
    end
    return minimum_dist
end

function min_euclidean(rect::Rect{N,T1}, p::Union{Rect{N,T2},Vec{N,T2}}) where {N,T1,T2}
    return sqrt(min_euclideansq(rect, p))
end

function task_map_emb(::EmbRep, ::EmbRep, xme, task_map::DistanceSphereToBox{ℝ{m},R1}) where m
    box = task_map.box
    center = @. box.origin + box.widths/2
    xme_rectframe  = task_map.rotation' * (xme - center) + center
    SA[min_euclidean(box, Vec(xme_rectframe)) - task_map.radius]
end

# Updateable version
struct DistanceSphereToBoxDynamic{M,R1,S,m} <: TaskMap{M,R1,S}
    radius::S
    box::Vector{Rect{m,S}}
    rotation::Vector{RotMatrix{m,S}}
end
DistanceSphereToBoxDynamic{M,R1}(r, b::Vector{Rect{m,S}}, R) where {M,R1,S,m} =
    DistanceSphereToBoxDynamic{M,R1,S,dim(M)}(r, b, R)

struct DistanceSphereToBoxDynamicT{M,R1,S,m} <: TaskMapT{M,R1,S}
    base_map::DistanceSphereToBoxDynamic{M,R1,S,m}
end
DistanceSphereToBoxDynamicT{M,R1}(r, b::Vector{Rect{m,S}}, R) where {M,R1,S,m} =
    DistanceSphereToBoxDynamicT{M,R1,S,m}(DistanceSphereToBoxDynamic{M,R1}(r, b, R))

function task_map_emb(::EmbRep, ::EmbRep, xme,
        task_map::DistanceSphereToBoxDynamic{ℝ{m},R1}) where m
    box = task_map.box[1]
    center = @. box.origin + box.widths/2
    xme_rectframe  = task_map.rotation[1]' * (xme - center) + center
    SA[min_euclidean(box, Vec(xme_rectframe)) - task_map.radius]
end