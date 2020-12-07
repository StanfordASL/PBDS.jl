# Map onto SE(3) embedded in R12
if !(@isdefined LinkFrameTransform)
    @computed struct LinkFrameTransform{Rn,R12,S,JC} <: TaskMap{Rn,R12,S}
        state_cache::StateCache{S,JC}
        link_cache::fulltype(LinkFrameCache{dim(Rn),S})
    end
end
struct LinkFrameTransformT{Rn,R12,S,JC} <: TaskMapT{Rn,R12,S}
    base_map::LinkFrameTransform{Rn,R12,S,JC}
end
LinkFrameTransform{ℝ{n},R12}(state_cache::StateCache{S,JC},
        link_cache::LinkFrameCache{n,S}) where {n,R12,S,JC} =
    LinkFrameTransform{ℝ{n},R12,S,JC}(state_cache, link_cache)

task_map_emb(::EmbRep, ::EmbRep, xme, task_map::LinkFrameTransform{ℝ{n},R12,S,JC}) where {n,S,JC} =
    link_frame_transform(xme, task_map.state_cache, task_map.link_cache.link)

function task_jacobian_emb(xme, task_map::LinkFrameTransform{ℝ{n},R12,S,JC},
        arg...) where {n,S,JC}
    state_cache, link_cache = task_map.state_cache, task_map.link_cache
    Jfc = link_cache.transform_jacobian
    update_link_transform_jacobian!(Jfc, xme, state_cache, link_cache.link)
    Jf = Jfc.data[1]
end

task_jacobian_chart(xm, task_map::LinkFrameTransform{ℝ{n},R12,S,JC},
    CM::Chart{1,ℝ{n}}, CN::Chart{1,R12}) where {n,S,JC} = task_jacobian_emb(xm, task_map)

function task_jacobian_emb_dot(xme, vme, task_map::LinkFrameTransform{ℝ{n},R12,S,JC},
        arg...) where {n,S,JC}
    state_cache, link_cache = task_map.state_cache, task_map.link_cache
    Jfc_dot = link_cache.transform_jacobian_dot
    update_link_transform_jacobian_dot!(Jfc_dot, xme, vme, state_cache, link_cache.link)
    Jf_dot = Jfc_dot.data[1]
end

task_jacobian_chart_dot(xm, vm, task_map::LinkFrameTransform{ℝ{n},R12,S,JC},
    CM::Chart{1,ℝ{n}}, CN::Chart{1,R12}) where {n,S,JC} = task_jacobian_emb_dot(xm, vm, task_map)