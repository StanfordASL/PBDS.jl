# Note: Not applicable to permuted task map scheme
task_forces_chart(pn, wn, task::TaskX,  C) =
    potential_force_chart(pn, task, C) + dissipative_forces_chart(pn, wn, task, C)

## Potential force
function potential_force_chart(pn, task::TaskX, C::Chart)
    ∇Φ = ForwardDiff.gradient(pn -> potential_chart(pn, task, C), pn)
    -∇Φ
end

function potential_from_home_chart(pni, task, Ci)
    Ch = home_task_chart(task)
    pnh = chart_transition(pni, Ci, Ch)
    Φ = potential_chart(pnh, task, Ch)
end

function potential_from_emb(pn, task, C)
    pne = chart_to_emb(pn, C)
    Φ = potential_emb(pne, task)
end

potential_coord_rep(task::TaskX) = default_coord_rep(task)
potential_chart(pn, task, C) = potential_chart(potential_coord_rep(task), pn, task, C)
potential_chart(::EmbRep, pn, task, C) = potential_from_emb(pn, task, C)
potential_chart(::ChartRep, pn, task, C) = potential_from_home_chart(pn, task, C)
potential_emb(pne, task) = throw("Not defined!")

default_potential(pn, task::TaskX, C::Chart) = zero(eltype(pn))

## Dissipative forces
function dissipative_forces_from_home_chart(pni, wni, task, Ci)
    Ch = home_task_chart(task)
    pnh, wnh = chart_transition_differential(pni, wni, Ci, Ch)
    ∂pnh_∂pni = chart_transition_jacobian(pni, Ci, Ch)
    ℱ_dish = dissipative_forces_chart(pnh, wnh, task, Ch)
    ℱ_disi = ∂pnh_∂pni'*ℱ_dish
end

function dissipative_forces_from_emb(pn, wn, task, C)
    pne, wne = chart_to_emb_differential(pn, wn, C)
    ℱ_dise = dissipative_forces_emb(pne, wne, task)
    ∂pne_∂pn = chart_to_emb_jacobian(pn, C)
    ℱ_dis = ∂pne_∂pn'*ℱ_dise
end

dissipative_forces_coord_rep(task::TaskX) = default_coord_rep(task)
dissipative_forces_chart(pn, wn, task, C) = 
    dissipative_forces_chart(dissipative_forces_coord_rep(task), pn, wn, task, C)
dissipative_forces_chart(::EmbRep, pn, wn, task, C) = dissipative_forces_from_emb(pn, wn, task, C)
dissipative_forces_chart(::ChartRep, pn, wn, task, C) =
    dissipative_forces_from_home_chart(pn, wn, task, C)
dissipative_forces_emb(pne, wne, task) = throw("Not defined!")

default_dissipative_forces(pn, wn, task::TaskX, C) = zero(eltype(pn)) * wn

## GDS forces
function dissipative_term_from_home_chart(xni, vni, task::TaskGDS, Ci::Chart)
    Ch = home_task_chart(task)
    xnh, vnh = chart_transition_differential(xni, vni, Ci, Ch)
    ∂xnh_∂xni = chart_transition_jacobian(xni, Ci, Ch)
    Bh = dissipative_term_chart(xnh, vnh, task, Ch)
    Bi = ∂xnh_∂xni'*Bh*∂xnh_∂xni
end

function dissipative_term_from_emb(xn, vn, task, C)
    xne, vne = chart_to_emb_differential(xn, vn, C)
    Be = dissipative_term_emb(xne, vne, task)
    ∂xne_∂xn = chart_to_emb_jacobian(xn, C)
    B = ∂xne_∂xn'*Be*∂xne_∂xn
end

dissipative_term_coord_rep(task::TaskGDS) = default_coord_rep(task)
dissipative_term_chart(xn, vn, task, C) = 
    dissipative_term_chart(dissipative_term_coord_rep(task), xn, vn, task, C)
dissipative_term_chart(::EmbRep, xn, vn, task, C) = dissipative_term_from_emb(xn, vn, task, C)
dissipative_term_chart(::ChartRep, xn, vn, task, C) =
    dissipative_term_from_home_chart(xn, vn, task, C)
dissipative_term_emb(xne, vne, task) = throw("Not defined!")

default_dissipative_term_emb(xn, vn, task::TaskGDS{<:TaskMapT{M,N}}, C) where {M,N} =    
    @SMatrix zeros(eltype(xn), embdim(N), embdim(N))
default_dissipative_term_chart(xn, vn, task::TaskGDS{<:TaskMapT{M,N}}, C) where {M,N} =    
    @SMatrix zeros(eltype(xn), dim(N), dim(N))