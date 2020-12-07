include(joinpath("caches", "link_transform.jl"))
include(joinpath("caches", "frame_point.jl"))

struct ControllerCache{n,S}
    link_caches::Vector{LinkFrameCache{n,S}}
    point_caches::Vector{FramePointCache{n,S}}
end
ControllerCache(link_caches::Vector{<:LinkFrameCache{n,S}}, point_caches) where {n,S} =
    ControllerCache{n,S}(link_caches, point_caches)

function setdirty!(cache::ControllerCache)
    map(setdirty!, cache.link_caches)
    map(setdirty!, cache.point_caches)
    nothing
end