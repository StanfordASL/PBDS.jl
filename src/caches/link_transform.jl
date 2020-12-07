# n is number of robot joints
if !(@isdefined LinkFrameCache)
    @computed struct LinkFrameCache{n,S}
        link::RigidBody{S}
        path::RigidBodyDynamics.Graphs.TreePath
        transform_jacobian::CacheElement{Vector{SMatrix{12,n,S,n*12}}}
        transform_jacobian_dot::CacheElement{Vector{SMatrix{12,n,S,n*12}}}
    end
end

function LinkFrameCache(link, state::MechanismState{X,M,C,JC}) where {X,M,C,JC}
    n = num_positions(state)
    link_path = path(state.mechanism, root_body(state.mechanism), link)
    transform_jacobian = CacheElement([@SMatrix zeros(12,n)])
    transform_jacobian_dot = CacheElement([@SMatrix zeros(12,n)])
    LinkFrameCache{n,X}(link, link_path, transform_jacobian, transform_jacobian_dot)
end

function setdirty!(link_cache::LinkFrameCache)
    CustomCollections.setdirty!(link_cache.transform_jacobian)
    CustomCollections.setdirty!(link_cache.transform_jacobian_dot)
    nothing
end

function link_frame_transform(xm, state_cache, link)
    state = state_cache[eltype(xm)]
    set_configuration!(state, xm)
    T = transform_to_root(state, link)
    R, t = SVector{9,eltype(xm)}(reshape(rotation(T),9)), SVector{3,eltype(xm)}(translation(T))
    [R; t]
end

@inline function update_link_transform_jacobian!(transform_jacobian::CacheElement, xm,
        state_cache, link)
    isdirty(transform_jacobian) && _update_link_transform_jacobian!(transform_jacobian, xm,
        state_cache, link)
    nothing
end

function _update_link_transform_jacobian!(transform_jacobian::CacheElement, xm, state_cache, link)
    f = (xm) -> link_frame_transform(xm, state_cache, link)
    transform_jacobian.data[1] = ForwardDiff.jacobian(f, xm)
    transform_jacobian.dirty = false
    nothing
end

@inline function update_link_transform_jacobian_dot!(transform_jacobian_dot::CacheElement,
        xm, vm, state_cache, link)
    isdirty(transform_jacobian_dot) && _update_link_transform_jacobian_dot!(transform_jacobian_dot,
        xm, vm, state_cache, link)
    nothing
end

function _update_link_transform_jacobian_dot!(transform_jacobian_dot::CacheElement,
        xm, vm, state_cache, link)
    n, m = size(transform_jacobian_dot.data[1])
    f = (xm) -> link_frame_transform(xm, state_cache, link)
    Jf = (xm) -> ForwardDiff.jacobian(f, xm)
    ∇xf = ForwardDiff.jacobian(xm -> reshape(Jf(xm), m*n), xm)
    Jfdot = SMatrix{n,m,eltype(xm)}(reshape(∇xf*vm, n, m))::SArray{Tuple{n,m},eltype(xm),2,m*n}

    transform_jacobian_dot.data[1] = Jfdot
    transform_jacobian_dot.dirty = false
    nothing
end