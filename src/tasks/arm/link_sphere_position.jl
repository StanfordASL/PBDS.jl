# Map from SE(3) embedded in R12
if !(@isdefined LinkSpherePosition)
    @computed struct LinkSpherePosition{R12,R3,S,n} <: TaskMap{R12,R3,S}
        point_cache::fulltype(FramePointCache{n,S})
    end
end
struct LinkSpherePositionT{R12,R3,S,n} <: TaskMapT{R12,R3,S}
    base_map::LinkSpherePosition{R12,R3,S,n}
end
LinkSpherePosition{R12,R3}(point_cache::FramePointCache{n,S}) where {n,S} =
    LinkSpherePosition{R12,R3,S,n}(point_cache)

task_map_emb(::EmbRep, ::EmbRep, xme, task_map::LinkSpherePosition{R12,R3,S}) where {n,S} =
    point_world_position(xme, task_map.point_cache.point.v)

function task_jacobian_emb(xme, task_map::LinkSpherePosition{R12,R3,S}, arg...) where {n,S}
    Jfc = task_map.point_cache.world_position_jacobian
    point_frame_position = task_map.point_cache.point.v
    update_point_world_position_jacobian!(Jfc, xme, point_frame_position)
    Jf = Jfc.data[1]
end

task_jacobian_chart(xm, task_map::LinkSpherePosition{R12,R3,S},
    CM::Chart{1,R12}, CN::Chart{1,R3}) where S = task_jacobian_emb(xm, task_map)

function task_jacobian_emb_dot(xme, vme, task_map::LinkSpherePosition{R12,R3,S},
        arg...) where {n,S}
    Jfc_dot = task_map.point_cache.world_position_jacobian_dot
    point_frame_position = task_map.point_cache.point.v
    update_point_world_position_jacobian_dot!(Jfc_dot, xme, vme, point_frame_position)
    Jf_dot = Jfc_dot.data[1]
end
task_jacobian_chart_dot(xm, vm, task_map::LinkSpherePosition{R12,R3,S},
    CM::Chart{1,R12}, CN::Chart{1,R3}) where S = task_jacobian_emb_dot(xm, vm, task_map)