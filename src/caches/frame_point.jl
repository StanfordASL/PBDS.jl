if !(@isdefined FramePointCache)
    @computed struct FramePointCache{n,S}
        link_cache::fulltype(LinkFrameCache{n,S})
        point::Point3D{SVector{3,S}}
        world_position_jacobian::CacheElement{Vector{SMatrix{3,12,S,36}}}
        world_position_jacobian_dot::CacheElement{Vector{SMatrix{3,12,S,36}}}
    end
end

function FramePointCache(link_cache::LinkFrameCache{n,S}, position_on_link) where {n,S}
    point = Point3D(default_frame(link_cache.link), position_on_link)
    world_position_jacobian = CacheElement([@SMatrix zeros(3,12)])
    world_position_jacobian_dot = CacheElement([@SMatrix zeros(3,12)])
    FramePointCache{n,S}(link_cache, point, world_position_jacobian, world_position_jacobian_dot)
end

function setdirty!(point_cache::FramePointCache)
    CustomCollections.setdirty!(point_cache.world_position_jacobian)
    CustomCollections.setdirty!(point_cache.world_position_jacobian_dot)
    nothing
end

# Frame transform and position in frame to world frame position
function point_world_position(Rt, point_frame_position)
    R = SMatrix{3,3,eltype(Rt)}(Rt[static(1):static(9)])
    t = Rt[static(10):static(12)]
    R*point_frame_position + t
end

function point_world_position_jacobian(Rt, point_frame_position)
    f = (Rt) -> point_world_position(Rt, point_frame_position)
    ForwardDiff.jacobian(f, Rt)
end

function point_world_position_jacobian_dot(Rt, Rt_dot, point_frame_position)
    m, n = 12, 3
    ∇xf = SMatrix{n*m,m,eltype(Rt)}(ForwardDiff.jacobian(Rt -> reshape(
        point_world_position_jacobian(Rt, point_frame_position), m*n), Rt))
    Jfdot = reshape(∇xf*Rt_dot, Size(n,m))
end

@inline function update_point_world_position_jacobian!(world_position_jacobian::CacheElement, Rt,
        point_frame_position)
    isdirty(world_position_jacobian) && _update_point_world_position_jacobian!(world_position_jacobian, Rt, point_frame_position)
    nothing
end

function _update_point_world_position_jacobian!(world_position_jacobian::CacheElement, Rt,
        point_frame_position)
    f = (Rt) -> point_world_position(Rt, point_frame_position)
    world_position_jacobian.data[1] = ForwardDiff.jacobian(f, Rt)
    world_position_jacobian.dirty = false
    nothing
end

@inline function update_point_world_position_jacobian_dot!(
    world_position_jacobian_dot::CacheElement, Rt, Rt_dot, point_frame_position)
    isdirty(world_position_jacobian_dot) && _update_point_world_position_jacobian_dot!(world_position_jacobian_dot, Rt, Rt_dot, point_frame_position)
    nothing
end

function _update_point_world_position_jacobian_dot!(world_position_jacobian_dot::CacheElement,
        Rt, Rt_dot, point_frame_position)
    n, m = size(world_position_jacobian_dot.data[1])
    f = (Rt) -> point_world_position(Rt, point_frame_position)
    Jf = (Rt) -> ForwardDiff.jacobian(f, Rt)
    ∇xf = ForwardDiff.jacobian(Rt -> reshape(Jf(Rt), m*n), Rt)
    Jfdot = SMatrix{n,m,eltype(Rt)}(reshape(∇xf*Rt, n, m))::SArray{Tuple{n,m},eltype(Rt),2,m*n}

    world_position_jacobian_dot.data[1] = Jfdot
    world_position_jacobian_dot.dirty = false
    nothing
end