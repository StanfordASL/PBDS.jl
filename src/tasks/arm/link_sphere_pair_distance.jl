if !(@isdefined LinkSpherePairDistance)
    @computed struct LinkSpherePairDistance{Rm,R1,S,JC} <: TaskMap{Rm,R1,S}
        state::MechanismState{S,S,S,JC}
        state_cache::StateCache{S,JC}
        point1_cache::fulltype(FramePointCache{dim(Rm),S})
        point2_cache::fulltype(FramePointCache{dim(Rm),S})
        radius1::S
        radius2::S
    end
end
struct LinkSpherePairDistanceT{Rm,R1,S,JC} <: TaskMapT{Rm,R1,S}
    base_map::LinkSpherePairDistance{Rm,R1,S,JC}
end
function LinkSpherePairDistance{ℝ{m},R1}(state::MechanismState{S,S,S,JC}, state_cache,
        point1_cache::FramePointCache{m,S}, point2_cache::FramePointCache{m,S},
        radius1, radius2) where {m,R1,S,JC}
    LinkSpherePairDistance{ℝ{m},R1,S,JC}(state, state_cache, point1_cache, point2_cache,
        radius1, radius2)
end

function task_map_emb(::EmbRep, ::EmbRep, xme, task_map::LinkSpherePairDistance{ℝ{m},R1}) where m
    rframe = root_frame(task_map.state.mechanism)
    x1 = transform(task_map.state, task_map.point1_cache.point, rframe).v
    x2 = transform(task_map.state, task_map.point2_cache.point, rframe).v
    SA[norm(x2 - x1) - task_map.radius1 - task_map.radius2]
end

function task_jacobian_emb(xme, task_map::LinkSpherePairDistance{ℝ{m},R1}, arg...) where m
    state_cache = task_map.state_cache
    point1_cache, point2_cache = task_map.point1_cache, task_map.point2_cache
    link1_cache, link2_cache = point1_cache.link_cache, point2_cache.link_cache
    link1, link2 = link1_cache.link, link2_cache.link
    point1_frame_position, point2_frame_position = point1_cache.point.v, point2_cache.point.v
    radius1, radius2 = task_map.radius1, task_map.radius2

    Jl1c, Jl2c = link1_cache.transform_jacobian, link2_cache.transform_jacobian
    update_link_transform_jacobian!(Jl1c, xme, state_cache, link1)
    update_link_transform_jacobian!(Jl2c, xme, state_cache, link2)
    Jl1, Jl2 = Jl1c.data[1], Jl2c.data[1]

    Jp1c, Jp2c = point1_cache.world_position_jacobian, point2_cache.world_position_jacobian
    Rt1 = link_frame_transform(xme, state_cache, link1)
    Rt2 = link_frame_transform(xme, state_cache, link2)
    update_point_world_position_jacobian!(Jp1c, Rt1, point1_frame_position)
    update_point_world_position_jacobian!(Jp2c, Rt2, point2_frame_position)
    Jp1, Jp2 = Jp1c.data[1], Jp2c.data[1]

    Jlp1, Jlp2 = Jp1*Jl1, Jp2*Jl2
    Jlp = [Jlp1; Jlp2]

    point1_world_position = point_world_position(Rt1, point1_frame_position)
    point2_world_position = point_world_position(Rt2, point2_frame_position)
    world_positions = [point1_world_position; point2_world_position]
    Jn = stacked_position_sphere_distance_jacobian(world_positions, radius1, radius2)

    Jf = Jn*Jlp
end
task_jacobian_chart(xm, task_map::LinkSpherePairDistance{ℝ{m},R1},
    CM::Chart{1,ℝ{m}}, CN::Chart{1,R1}) where m = task_jacobian_emb(xm, task_map)

function task_jacobian_emb_dot(xme, vme, task_map::LinkSpherePairDistance{ℝ{m},R1}, arg...) where m
    state_cache = task_map.state_cache
    point1_cache, point2_cache = task_map.point1_cache, task_map.point2_cache
    link1_cache, link2_cache = point1_cache.link_cache, point2_cache.link_cache
    link1, link2 = link1_cache.link, link2_cache.link
    point1_frame_position, point2_frame_position = point1_cache.point.v, point2_cache.point.v
    radius1, radius2 = task_map.radius1, task_map.radius2

    Jl1c, Jl2c = link1_cache.transform_jacobian, link2_cache.transform_jacobian
    Jl1_dotc, Jl2_dotc = link1_cache.transform_jacobian_dot, link2_cache.transform_jacobian_dot
    update_link_transform_jacobian!(Jl1c, xme, state_cache, link1)
    update_link_transform_jacobian!(Jl2c, xme, state_cache, link2)
    update_link_transform_jacobian_dot!(Jl1_dotc, xme, vme, state_cache, link1)
    update_link_transform_jacobian_dot!(Jl2_dotc, xme, vme, state_cache, link2)
    Jl1, Jl2 = Jl1c.data[1], Jl2c.data[1]
    Jl1_dot, Jl2_dot = Jl1_dotc.data[1], Jl2_dotc.data[1]

    Jp1c, Jp2c = point1_cache.world_position_jacobian, point2_cache.world_position_jacobian
    Jp1_dotc = point1_cache.world_position_jacobian_dot
    Jp2_dotc = point2_cache.world_position_jacobian_dot
    Rt1 = link_frame_transform(xme, state_cache, link1)
    Rt2 = link_frame_transform(xme, state_cache, link2)
    Rt1_dot, Rt2_dot = Jl1*vme, Jl2*vme
    update_point_world_position_jacobian!(Jp1c, Rt1, point1_frame_position)
    update_point_world_position_jacobian!(Jp2c, Rt2, point2_frame_position)
    update_point_world_position_jacobian_dot!(Jp1_dotc, Rt1, Rt1_dot, point1_frame_position)
    update_point_world_position_jacobian_dot!(Jp2_dotc, Rt2, Rt2_dot, point2_frame_position)
    Jp1, Jp2 = Jp1c.data[1], Jp2c.data[1]
    Jp1_dot, Jp2_dot = Jp1_dotc.data[1], Jp2_dotc.data[1]

    Jlp1, Jlp2 = Jp1*Jl1, Jp2*Jl2
    Jlp1_dot = Jp1_dot*Jl1 + Jp1*Jl1_dot
    Jlp2_dot = Jp2_dot*Jl2 + Jp2*Jl2_dot

    Jlp = [Jlp1; Jlp2]
    Jlp_dot = [Jlp1_dot; Jlp2_dot]

    point1_world_position = point_world_position(Rt1, point1_frame_position)
    point2_world_position = point_world_position(Rt2, point2_frame_position)
    point1_world_velocity, point2_world_velocity = Jlp1*vme, Jlp2*vme
    world_positions = [point1_world_position; point2_world_position]
    world_velocities = [point1_world_velocity; point2_world_velocity]
    Jn = stacked_position_sphere_distance_jacobian(world_positions, radius1, radius2)
    Jn_dot = stacked_position_sphere_distance_jacobian_dot(world_positions, world_velocities,
        radius1, radius2)

    Jf_dot = Jn_dot*Jlp + Jn*Jlp_dot
end
task_jacobian_chart_dot(xm, task_map::LinkSpherePairDistance{ℝ{m},R1},
    CM::Chart{1,ℝ{m}}, CN::Chart{1,R1}) where m = task_jacobian_emb_dot(xm, task_map)

# Distances between points, with positions input as a single stacked vector 
function stacked_position_distance(positions)
    SA[norm(positions[static(1):static(3)] - positions[static(4):static(6)])]
end

function stacked_position_distance_jacobian(positions)
    f = (positions) -> stacked_position_distance(positions)
    ForwardDiff.jacobian(f, positions)
end

function stacked_position_distance_jacobian_dot(xs, xs_dot)
    m, n = 6, 1
    ∇xf = ForwardDiff.jacobian(xs -> stacked_position_distance_jacobian(xs),
        xs)::SArray{Tuple{m,m},eltype(xs),2,m*m}
    Jfdot = reshape(∇xf*xs_dot, Size(n,m))
end

# Distances between spheres, with positions input as a single stacked vector 
function stacked_position_sphere_distance(positions, radius1, radius2)
    SA[norm(positions[static(1):static(3)] - positions[static(4):static(6)]) - radius1 - radius2]
end

function stacked_position_sphere_distance_jacobian(positions, radius1, radius2)
    f = (positions) -> stacked_position_sphere_distance(positions, radius1, radius2)
    ForwardDiff.jacobian(f, positions)
end

function stacked_position_sphere_distance_jacobian_dot(xs, xs_dot, radius1, radius2)
    m, n = 6, 1
    ∇xf = ForwardDiff.jacobian(xs -> stacked_position_sphere_distance_jacobian(xs, radius1,
        radius2), xs)::SArray{Tuple{m,m},eltype(xs),2,m*m}
    Jfdot = reshape(∇xf*xs_dot, Size(n,m))
end


# For verification
# A bit slower than custom function after precomputation and caching
function task_jacobian_emb2(xme, task_map::LinkSpherePairDistance{ℝ{m},R1}, arg...) where m
    S = eltype(xme)
    state = task_map.state
    rframe = root_frame(state.mechanism)
    point1_cache, point2_cache = task_map.point1_cache, task_map.point2_cache
    link1_cache, link2_cache = point1_cache.link_cache, point2_cache.link_cache
    point1, point2 = point1_cache.point, point2_cache.point
    point1_world = transform(state, point1, rframe)
    point2_world = transform(state, point2, rframe)
    x1, x2 = point1_world.v, point2_world.v

    Jf1 = SMatrix{3,m,S}(point_jacobian(state, link1_cache.path, point1_world).J)
    Jf2 = SMatrix{3,m,S}(point_jacobian(state, link2_cache.path, point2_world).J)

    SMatrix{1,m,S}((Jf2' - Jf1')*(x2 - x1)/norm(x2-x1))
end

function task_map_emb_for_jdot(xme, task_map::LinkSpherePairDistance{ℝ{m},R1}) where m
    state = task_map.state_cache[eltype(xme)]
    rframe = root_frame(state.mechanism)
    set_configuration!(state, xme)
    x1 = transform(state, task_map.point1_cache.point, rframe).v
    x2 = transform(state, task_map.point2_cache.point, rframe).v
    SA[norm(x2 - x1) - task_map.radius1 - task_map.radius2]
end 
task_jacobian_emb_for_jdot(xme, task_map::LinkSpherePairDistance{ℝ{m},R1}) where m =
    ForwardDiff.jacobian(xme -> task_map_emb_for_jdot(xme, task_map), xme)

function task_jacobian_emb_dot2(xme, vme, task_map::LinkSpherePairDistance{ℝ{m},R1}, arg...) where m
    n, S = embdim(codomain_manifold(task_map)), eltype(xme)
    ∇xf = ForwardDiff.jacobian(xme -> SVector{m*n,S}(
        reshape(task_jacobian_emb_for_jdot(xme, task_map), m*n)),
        xme)::SArray{Tuple{m*n,m},S,2,m*m*n}
    Jfdot_emb = reshape(∇xf*vme, Size(n,m))
end

task_jacobian_chart_dot(xm, vm, task_map::LinkSpherePairDistance{ℝ{m},R1},
    CM::Chart{1,ℝ{m}}, CN::Chart{1,R1}) where m = task_jacobian_emb_dot(xm, vm, task_map)