function apply_task_metric_chart(pn, vn, wn, task::Task{<:TaskMap{M,N,S}},
        C::Chart{I,N}) where {M,N,S,I}
    G = metric_chart(pn, task, C)
    vn'*G*wn
end

task_riemannian_norm_chart(pn, vn, task::Task{<:TaskMap{M,N,S}}, C::Chart{I,N}) where {M,N,S,I} =
    sqrt(apply_task_metric(pn, vn, vn, task, C))

function metric_chart_transition(p1, G1, C1::Chart{I,M}, C2::Chart{J,M}) where {M,I,J}
    p2 = chart_transition(p1, C1, C2)
    ∂p1_∂p2 = chart_transition_jacobian(p2, C2, C1)
    G2 = ∂p1_∂p2'*G1*∂p1_∂p2
end

function metric_from_home_chart(pni, task::Task{<:TaskMapT}, Ci)
    Ch = home_task_chart(task, Ci)
    pnh = chart_transition(pni, T(Ci), T(Ch))
    Gh = metric_chart(pnh, task, Ch)
    Gi = metric_chart_transition(pnh, Gh, T(Ch), T(Ci))
end

function metric_from_home_chart(pni, task::Task{<:BaseTaskMap}, Ci)
    Ch = home_task_chart(task, Ci)
    pnh = chart_transition(pni, Ci, Ch)
    Gh = metric_chart(pnh, task, Ch)
    Gi = metric_chart_transition(pnh, Gh, Ch, Ci)
end

function metric_from_emb(pn, task::Task{<:BaseTaskMap}, C)
    pne = chart_to_emb(pn, C)
    Ge = metric_emb(pne, task)
    ∂pne_∂pn = chart_to_emb_jacobian(pn, C)
    G = ∂pne_∂pn'*Ge*∂pne_∂pn
end

metric_coord_rep(task::Task{<:Union{BaseTaskMap,TaskMapT}}) = default_coord_rep(task)
metric_chart(pni, task::Task{<:Union{BaseTaskMap,TaskMapT}}, Ci::Chart) = 
    metric_chart(metric_coord_rep(task), pni, task, Ci)
metric_chart(::EmbRep, pni, task, Ci) = metric_from_emb(pni, task, Ci)
metric_chart(::ChartRep, pni, task, Ci) = metric_from_home_chart(pni, task, Ci)
metric_emb(pne, task) = throw("Not defined!")

default_metric(pni, task::Task{<:TaskMap{M,ℝ{n},S}}, CN1::Chart{1,ℝ{n}}) where {M,n,S} =
    SMatrix{n,n,eltype(pni)}(I)
default_metric(pni, task::Task{<:TaskMap{M,S1,S}}, CN1::Chart{<:AngleChart,S1}) where {M,S} =
    SMatrix{1,1,eltype(pni)}(I)
default_metric(pne, task::Task{<:TaskMap{M,𝕊{n},S}}) where {M,n,S} =
    SMatrix{n+1,n+1,eltype(pne)}(I)

## Weight metrics
weight_metric_coord_rep(task::Task{<:Union{BaseTaskMap,TaskMapT}}) = default_coord_rep(task)

function weight_metric_from_home_chart(pni, wni, task::Task{<:BaseTaskMap}, Ci)
    Ch = home_task_chart(task, Ci)
    pnh, wnh = chart_transition_differential(pni, wni, Ci, Ch)
    Wh = weight_metric_chart(pnh, wnh, task, Ch)
    Wi = metric_chart_transition(pnh, Wh, Ch, Ci)
end

function weight_metric_from_emb(pn, wn, task::Task{<:BaseTaskMap}, C)
    pne, wne = chart_to_emb_differential(pn, wn, C)
    We = weight_metric_emb(pne, wne, task)
    ∂pne_∂pn = chart_to_emb_jacobian(pn, C)
    W = ∂pne_∂pn'*We*∂pne_∂pn
end

weight_metric_chart(pn, wn, task::Task{<:Union{BaseTaskMap,TaskMapT}}, C::Chart) =
    weight_metric_chart(weight_metric_coord_rep(task), pn, wn, task, C)
weight_metric_chart(::EmbRep, pn, wn, task, C) = weight_metric_from_emb(pn, wn, task, C)
weight_metric_chart(::ChartRep, pn, wn, task, C) =
    weight_metric_from_home_chart(pn, wn, task, C)
weight_metric_emb(pne, wne, task) = throw("Not defined!")

# function active_weight_position_from_home_chart(pn, task::Task{<:BaseTaskMap}, C)
#     Ch = home_task_chart(task, C)
#     pnh = chart_transition(pni, Ci, Ch)
#     active_weight_position_chart(pnh, task, Ch)
# end

# function active_weight_position_from_emb(pn, task::Task{<:BaseTaskMap}, C)
#     pne = chart_to_emb(pn, C)
#     active_weight_position_emb(pne, task)
# end

# active_weight_position_chart(pn, task::Task{<:Union{BaseTaskMap,TaskMapT}}, C::Chart) =
#     active_weight_position_chart(weight_metric_coord_rep(task), pn, C)
# active_weight_position_chart(::EmbRep, pn, C) = active_weight_position_from_emb(pn, wn, task, C)
# active_weight_position_chart(::ChartRep, pn, C) = 
#     active_weight_position_from_home_chart(pn, task, C)
# active_weight_position_emb(pne, task) = throw("Not defined!")

active_weight_position_chart(pn, task::Task{<:Union{BaseTaskMap,TaskMapT}}, C::Chart) = true

default_weight_metric(pn, wn, task::Task{<:TaskMap{M,ℝ{n},S}}, CN::Chart{1,ℝ{n}}) where {M,n,S} =
    SMatrix{n,n,eltype(pn)}(I)
function default_weight_metric(pn, wn, task::Task{<:TaskMap{M,S1,S}},
        CN::Chart{<:AngleChart,S1}) where {M,S}
    SMatrix{1,1,eltype(pn)}(I)
end
default_weight_metric(pne, vne, task::Task{<:TaskMap{M,𝕊{n},S}}) where {M,n,S} =
    SMatrix{n+1,n+1,eltype(pne)}(I)

## GDS metrics
function metric_from_home_chart(xni, vni, task, Ci::Chart{I,N}) where {N,I}
    Ch = home_task_chart(task, Ci)
    xnh, vnh = chart_transition_differential(xni, vni, Ci, Ch)
    Gh = metric_chart(xnh, vnh, task, Ch)
    
    # Since the block of the metric used for GDS is on the horizontal bundle, cannot perform a chart
    # transition without considering the full metric
    n, nt = dim(N), dim(T{N})
    Gh_full = SMatrix{nt,nt,eltype(Gh)}([Gh zeros(n,n); zeros(n,n) I])
    pnh = [xnh; vnh]
    Gi_full = metric_chart_transition(pnh, Gh, T(Ch), T(Ci))
    Gi = SMatrix{n,n,eltype(Gh)}(Gi_full[1:n,1:n])

    # If projection was onto vertical bundle:
    # Gi = metric_chart_transition(xnh, Gh, Ch, Ci)
end

function metric_from_emb(xn, vn, task, C::Chart{J,N}) where {N,J}
    xne, vne = chart_to_emb_differential(xn, vn, C)
    Ge = metric_emb(xne, vne, task)
    # Since the block of the metric used for GDS is on the horizontal bundle, cannot perform a chart
    # transition without considering the full metric
    n, ne, nte = dim(N), embdim(N), embdim(T{N})
    Ge_full = SMatrix{nte,nte,eltype(Ge)}([Ge zeros(ne,ne); zeros(ne,ne) I])
    pn = [xn; vn]
    ∂pne_∂pn = chart_to_emb_jacobian(pn, T(C))
    G_full = ∂pne_∂pn'*Ge_full*∂pne_∂pn
    G = SMatrix{n,n,eltype(Ge)}(G_full[1:n,1:n])

    # If projection was onto vertical bundle:
    # ∂xne_∂xn = chart_to_emb_jacobian(xn, C)
    # G = ∂xne_∂xn'*Ge*∂xne_∂xn
end

metric_coord_rep(task::TaskGDS) = default_coord_rep(task)
metric_chart(xni, vni, task::TaskGDS, Ci) = metric_chart(metric_coord_rep(task), xni, vni, task, Ci)
metric_chart(::EmbRep, xni, vni, task, Ci) = metric_from_emb(xni, vni, task, Ci)
metric_chart(::ChartRep, xni, vni, task, Ci) = metric_from_home_chart(xni, vni, task, Ci)
metric_emb(xne, vne, task) = throw("Not defined!")

default_metric(xn, vn, task::TaskGDS{<:TaskMapT{M,ℝ{n},S}}, CN::Chart{1,ℝ{n}}) where {M,n,S} =
    SMatrix{n,n,eltype(xn)}(I)