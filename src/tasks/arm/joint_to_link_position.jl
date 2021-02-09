struct JointToLinkPosition{Rn,R3,S,JC} <: TaskMap{Rn,R3,S}
    state_cache::StateCache{S,JC}
    link::RigidBody{S}
end
struct JointToLinkPositionT{Rn,R3,S,JC} <: TaskMapT{Rn,R3,S}
    base_map::JointToLinkPosition{Rn,R3,S,JC}
end
JointToLinkPosition{ℝ{n},R3}(s::StateCache{S,JC}, l::RigidBody{S}) where {n,R3,S,JC} = 
    JointToLinkPosition{ℝ{n},R3,S,JC}(s,l)

function task_map_emb(::EmbRep, ::EmbRep, xme, task_map::JointToLinkPosition{ℝ{n},R3}) where n
    state = task_map.state_cache[eltype(xme)]
    set_configuration!(state, xme)
    translation(transform_to_root(state, task_map.link))
end