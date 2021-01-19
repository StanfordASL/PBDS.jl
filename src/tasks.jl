include(joinpath("tasks", "task_types.jl"))
include(joinpath("tasks", "product_tasks.jl"))
include(joinpath("tasks", "identity.jl"))
include(joinpath("tasks", "position_around_point.jl"))
include(joinpath("tasks", "distance_from_point.jl"))
include(joinpath("tasks", "distance_from_sphere_surface.jl"))
include(joinpath("tasks", "distance_sphere_to_sphere.jl"))
include(joinpath("tasks", "distance_sphere_to_cup.jl"))
include(joinpath("tasks", "distance_sphere_to_cylinder.jl"))
include(joinpath("tasks", "distance_sphere_to_box.jl"))
include(joinpath("tasks", "distance_sphere_to_line.jl"))
include(joinpath("tasks", "angular_distance.jl"))
include(joinpath("tasks", "angular_position_around_point.jl"))
include(joinpath("tasks", "coordinate.jl"))
include(joinpath("tasks", "coordinate_distance.jl"))
include(joinpath("tasks", "arm", "joint_to_link_position.jl"))
include(joinpath("tasks", "arm", "link_sphere_pair_distance.jl"))
include(joinpath("tasks", "arm", "link_sphere_position.jl"))
include(joinpath("tasks", "arm", "link_frame_transform.jl"))

Base.eltype(::Type{TaskMapX{A,B,S}}) where {A,B,S} = S

domain_manifold(task_map::Union{TaskMap{M,N},TaskMapT{M,N}}) where {M,N} = M
codomain_manifold(task_map::Union{TaskMap{M,N},TaskMapT{M,N}}) where {M,N} = N
domain_manifold(task::TaskX) = domain_manifold(task.task_map)
codomain_manifold(task::TaskX) = codomain_manifold(task.task_map)
domain_manifold(M::Manifold) = M

domain_coord_rep(task_map::TaskMapT) = domain_coord_rep(task_map.base_map)
codomain_coord_rep(task_map::TaskMapT) = codomain_coord_rep(task_map.base_map)

default_coord_rep(::TaskX) = ChartRep()

choose_chart_emb(xme, taskx::Union{TaskX,TaskMapX}, CM, CN) =
    choose_chart_emb(xme, base_task_map(taskx), CM, CN)
choose_chart_chart(xm, taskx::Union{TaskX,TaskMapX}, CM, CN) =
    choose_chart_chart(xm, base_task_map(taskx), CM, CN)

choose_chart_emb(xme, task_map::TaskMap, CM, CN) =
    choose_chart_emb(chart_choice_coord_rep(CN), xme, task_map, CM, CN)
choose_chart_chart(xm, task_map::TaskMap, CM, CN) =
    choose_chart_chart(chart_choice_coord_rep(CN), xm, task_map, CM, CN)

choose_chart_emb(::EmbRep, xme, task_map::TaskMap, CM, CN) =
    choose_chart_emb(EmbRep(), task_map_emb(xme, task_map, CM, CN), CN)
choose_chart_emb(::ChartRep, xme, task_map::TaskMap, CM, CN) =
    choose_chart_chart(ChartRep(), task_map_emb_chart(xme, task_map, CM, CN), CN)
choose_chart_chart(::EmbRep, xm, task_map::TaskMap, CM, CN) =
    choose_chart_emb(EmbRep(), task_map_chart_emb(xm, task_map, CM, CN), CN)
choose_chart_chart(::ChartRep, xm, task_map::TaskMap, CM, CN) =
    choose_chart_chart(ChartRep(), task_map_chart(xm, task_map, CM, CN), CN)

choose_chart(::EmbRep, xme, taskx::Union{TaskX,TaskMapX}, CM, CN) =
    choose_chart_emb(xme, taskx, CM, CN)
choose_chart(::ChartRep, xm, taskx::Union{TaskX,TaskMapX}, CM, CN) =
    choose_chart_chart(xm, taskx, CM, CN)

base_task_map(task::TaskX) = base_task_map(task.task_map)
base_task_map(task_map::TaskMap) = task_map
base_task_map(task_map::TaskMapX) = task_map.base_map

task_map_emb(pme, task_map::TaskMapX, arg...) =
    task_map_emb(domain_coord_rep(task_map), codomain_coord_rep(task_map), pme, task_map, arg...)
task_map_emb_chart(pme, task_map::TaskMapX, arg...) =
    task_map_emb_chart(domain_coord_rep(task_map), codomain_coord_rep(task_map),
    pme, task_map, arg...)
task_map_chart_emb(pm, task_map::TaskMapX, arg...) =
    task_map_chart_emb(domain_coord_rep(task_map), codomain_coord_rep(task_map),
    pm, task_map, arg...)
task_map_chart(pm, task_map::TaskMapX, arg...) =
    task_map_chart(domain_coord_rep(task_map), codomain_coord_rep(task_map), pm, task_map, arg...)

task_map_emb(pme, task::TaskX, arg...) = task_map_emb(pme, task.task_map, arg...)
task_map_emb_chart(pme, task::TaskX, arg...) = task_map_emb_chart(pme, task.task_map, arg...)
task_map_chart_emb(pm, task::TaskX, arg...) = task_map_chart_emb(pm, task.task_map, arg...)
task_map_chart(pm, task::TaskX, arg...) = task_map_chart(pm, task.task_map, arg...)

base_task_map_emb(pme, taskx::Union{TaskX,TaskMapX}, arg...) =
    task_map_emb(pme, base_task_map(taskx), arg...)
base_task_map_emb_chart(pme, taskx::Union{TaskX,TaskMapX}, arg...) =
    task_map_emb_chart(pme, base_task_map(taskx), arg...)
base_task_map_chart_emb(pm, taskx::Union{TaskX,TaskMapX}, arg...) =
    task_map_chart_emb(pm, base_task_map(taskx), arg...)
base_task_map_chart(pm, taskx::Union{TaskX,TaskMapX}, arg...) =
    task_map_chart(pm, base_task_map(taskx), arg...)

## No chart specified
# task_map_emb(::EmbRep, ::EmbRep, pme, task_map::TaskMap) # Should be defined for specific task
# Throw out given charts as needed; base emb map should be defined
task_map_emb(::EmbRep, ::EmbRep, pme, task_map::TaskMap, arg...) =
    task_map_emb(EmbRep(), EmbRep(), pme, task_map)
task_map_emb(::CoordinateRep, ::CoordinateRep, pme, task_map) =
    throw(ArgumentError("Chart(s) not specified"))

task_map_emb_chart(::CoordinateRep, ::CoordinateRep, pme, task_map) =
    throw(ArgumentError("Chart not specified"))
task_map_chart_emb(::CoordinateRep, ::CoordinateRep, pm, task_map) =
    throw(ArgumentError("Chart not specified"))
task_map_chart(::CoordinateRep, ::CoordinateRep, pm, task_map) =
    throw(ArgumentError("Charts not specified"))

## One chart specified
function task_map_emb(::EmbRep, ::ChartRep, pme, task_map::TaskMap, CN)
    pn = task_map_emb_chart(pme, task_map, CN)
    pne = chart_to_emb(pn, CN)
end
function task_map_emb(::ChartRep, ::EmbRep, pme, task_map::TaskMap, CM)
    pm = emb_to_chart(pme, CM)
    pne = task_map_chart_emb(pm, task_map, CM)
end
task_map_emb(::ChartRep, ::ChartRep, pme, task_map, C) =
    throw(ArgumentError("Second chart not specified"))

function task_map_emb_chart(::EmbRep, ::EmbRep, pme, task_map::TaskMap{M,N,S},
        CN::Chart{J,N}) where {M,N,S,J}
    pne = task_map_emb(pme, task_map)
    pn = emb_to_chart(pne, CN)
end
# Should be defined for specific task:
# task_map_emb_chart(::EmbRep, ::ChartRep, pm, task_map::TaskMap, CN) 
task_map_emb_chart(::ChartRep, ::EmbRep, pme, task_map, C) =
    throw(ArgumentError("Second chart not specified"))
task_map_emb_chart(::ChartRep, ::ChartRep, pme, task_map, C) =
    throw(ArgumentError("Second chart not specified"))

function task_map_chart_emb(::EmbRep, ::EmbRep, pm, task_map::TaskMap{M,N,S},
        CM::Chart{I,M}) where {M,N,S,I}
    pme = chart_to_emb(pm, CM)
    pne = task_map_emb(pme, task_map)
end
task_map_chart_emb(::EmbRep, ::ChartRep, pm, task_map, C) =
    throw(ArgumentError("Second chart not specified"))
# Should be defined for specific task:
# task_map_chart_emb(::ChartRep, ::EmbRep, pm, task_map::TaskMap, C) 
task_map_chart_emb(::ChartRep, ::ChartRep, pm, task_map, C) =
    throw(ArgumentError("Second chart not specified"))

task_map_chart(::CoordinateRep, ::CoordinateRep, pm, task_map, C) =
    throw(ArgumentError("Second chart not specified"))

## Both charts specified
task_map_emb(::EmbRep, ::ChartRep, pme, task_map::TaskMap, CM, CN) =
    task_map_emb(pme, task_map, CN)
task_map_emb(::ChartRep, ::EmbRep, pme, task_map::TaskMap, CM, CN) =
    task_map_emb(pme, task_map, CM)
function task_map_emb(::ChartRep, ::ChartRep, pme, task_map::TaskMap, CM, CN)
    pm = emb_to_chart(pme, CM)
    pn = task_map_chart(pm, task_map, CM, CN)
    pne = chart_to_emb(pn, CN)
end

task_map_emb_chart(::EmbRep, ::EmbRep, pme, task_map::TaskMap, CM, CN) =
    task_map_emb_chart(pme, task_map, CN)
task_map_emb_chart(::EmbRep, ::ChartRep, pme, task_map::TaskMap, CM, CN) =
    task_map_emb_chart(pme, task_map, CN)
function task_map_emb_chart(::ChartRep, ::EmbRep, pme, task_map::TaskMap{M,N,S}, CM,
        CN::Chart{J,N}) where {M,N,S,J}
    pm = emb_to_chart(pme, CM)
    pne = task_map_chart_emb(pm, task_map, CM)
    pn = emb_to_chart(pne, CN)
end
function task_map_emb_chart(::ChartRep, ::ChartRep, pme, task_map::TaskMap, CM, CN)
    pm = emb_to_chart(pme, CM)
    pn = task_map_chart(pm, task_map, CM, CN)
end

task_map_chart_emb(::EmbRep, ::EmbRep, pm, task_map::TaskMap, CM, CN) =
    task_map_chart_emb(pm, task_map, CM)
function task_map_chart_emb(::EmbRep, ::ChartRep, pm, task_map::TaskMap{M,N,S}, CM::Chart{I,M},
        CN) where {M,N,S,I}
    pme = chart_to_emb(pm, CM)
    pn = task_map_emb_chart(pme, task_map, CN)
    pne = chart_to_emb(pn, CN)
end
task_map_chart_emb(::ChartRep, ::EmbRep, pm, task_map::TaskMap, CM, CN) =
    task_map_chart_emb(pm, task_map, CM)
function task_map_chart_emb(::ChartRep, ::ChartRep, pm, task_map::TaskMap, CM, CN)
    pn = task_map_chart(pm, task_map, CM, CN)
    pne = chart_to_emb(pn, CN)
end

function task_map_chart(::EmbRep, ::EmbRep, pm, task_map::TaskMap{M,N,S}, CM::Chart{I,M},
        CN::Chart{J,N}) where {M,N,S,I,J}
    pme = chart_to_emb(pm, CM)
    pne = task_map_emb(pme, task_map)
    pn = emb_to_chart(pne, CN)
end
function task_map_chart(::EmbRep, ::ChartRep, pm, task_map::TaskMap{M,N,S}, CM::Chart{I,M},
        CN) where {M,N,S,I}
    pme = chart_to_emb(pm, CM)
    pn = task_map_emb_chart(pme, task_map, CN)
end
function task_map_chart(::ChartRep, ::EmbRep, pm, task_map::TaskMap{M,N,S}, CM,
        CN::Chart{J,N}) where {M,N,S,J}
    pne = task_map_chart_emb(pm, task_map, CM)
    pn = emb_to_chart(pne, CN)
end
# Should be defined for specific task
# task_map_chart(::ChartRep, ::ChartRep, pm, task_map::TaskMap, CM, CN)

function task_differential_map_emb(pme, wme, task_map::TaskMap, arg...)
    pne = task_map_emb(pme, task_map, arg...)
    ∂pne_∂pme = task_jacobian_emb(pme, task_map, arg...)
    wne = ∂pne_∂pme*wme
    pne, wne
end

function task_differential_map_emb_chart(pme, wme, task_map::TaskMap, arg...)
    pn = task_map_emb_chart(pme, task_map, arg...)
    ∂pn_∂pme = task_jacobian_emb_chart(pme, task_map, arg...)
    wn = ∂pn_∂pme*wme
    pn, wn
end

function task_differential_map_chart_emb(pm, wm, task_map::TaskMap, arg...)
    pne = task_map_chart_emb(pm, task_map, arg...)
    ∂pne_∂pm = task_jacobian_chart_emb(pm, task_map, arg...)
    wne = ∂pne_∂pm*wm
    pne, wne
end

function task_differential_map_chart(pm, wm, task_map::TaskMap, arg...)
    pn = task_map_chart(pm, task_map, arg...)
    ∂pn_∂pm = task_jacobian_chart(pm, task_map, arg...)
    wn = ∂pn_∂pm*wm
    pn, wn
end

task_differential_map_emb(pme, wme, task::TaskX, arg...) =
    task_differential_map_emb(pme, wme, base_task_map(task), arg...)
task_differential_map_emb_chart(pme, wme, task::TaskX, arg...) =
    task_differential_map_emb_chart(pme, wme, base_task_map(task), arg...) 
task_differential_map_chart_emb(pm, wme, task::TaskX, arg...) =
    task_differential_map_chart_emb(pm, wme, base_task_map(task), arg...)
task_differential_map_chart(pm, wme, task::TaskX, arg...) =
    task_differential_map_chart(pm, wme, base_task_map(task), arg...) 

## Tangent bundle task maps
function task_map_emb(pme, task_map::TaskMapT{M,N,S}, arg...) where {M,N,S}
    xme, vme = tangent_vector_emb_splitview(pme, T{M})
    xne, vne = task_differential_map_emb(xme, vme, task_map.base_map, arg...)
    pne = [xne; vne]
end

function task_map_emb_chart(pme, task_map::TaskMapT{M,N,S}, arg...) where {M,N,S}
    xme, vme = tangent_vector_emb_splitview(pme, T{M})
    xn, vn = task_differential_map_emb_chart(xme, vme, task_map.base_map, arg...)
    pn = [xn; vn]
end

function task_map_chart_emb(pm, task_map::TaskMapT{M,N,S}, arg...) where {M,N,S}
    xm, vm = tangent_vector_chart_splitview(pm, T{M})
    xne, vne = task_differential_map_chart_emb(xm, vm, task_map.base_map, arg...)
    pne = [xne; vne]
end

function task_map_chart(pm, task_map::TaskMapT{M,N,S}, arg...) where {M,N,S}
    xm, vm = tangent_vector_chart_splitview(pm, T{M})
    xn, vn = task_differential_map_chart(xm, vm, task_map.base_map, arg...)
    pn = [xn; vn]
end

## Jacobians
task_jacobian_emb(pme, task_map::TaskMapX, arg...) =
    ForwardDiff.jacobian(pme -> task_map_emb(pme, task_map, arg...), pme)
task_jacobian_emb_chart(pme, task_map::TaskMapX, arg...) =
    ForwardDiff.jacobian(pme -> task_map_emb_chart(pme, task_map, arg...), pme)
task_jacobian_chart_emb(pm, task_map::TaskMapX, arg...) =
    ForwardDiff.jacobian(pm -> task_map_chart_emb(pm, task_map, arg...), pm)
task_jacobian_chart(pm, task_map::TaskMapX, arg...) =
    ForwardDiff.jacobian(pm -> task_map_chart(pm, task_map, arg...), pm)

base_task_jacobian_emb(pme, taskx::Union{TaskX,TaskMapX}, arg...) = 
    task_jacobian_emb(pme, base_task_map(taskx), arg...)
base_task_jacobian_emb_chart(pme, taskx::Union{TaskX,TaskMapX}, arg...) = 
    task_jacobian_emb_chart(pme, base_task_map(taskx), arg...)
base_task_jacobian_chart_emb(pm, taskx::Union{TaskX,TaskMapX}, arg...) = 
    task_jacobian_chart_emb(pm, base_task_map(taskx), arg...)
base_task_jacobian_chart(pm, taskx::Union{TaskX,TaskMapX}, arg...) = 
    task_jacobian_chart(pm, base_task_map(taskx), arg...)

task_jacobian_emb(pme, task::TaskX, arg...) = task_jacobian_emb(pme, task.task_map, arg...)
task_jacobian_emb_chart(pme, task::TaskX, arg...) =
    task_jacobian_emb_chart(pme, task.task_map, arg...)
task_jacobian_chart_emb(pm, task::TaskX, arg...) =
    task_jacobian_chart_emb(pm, task.task_map, arg...)
task_jacobian_chart(pm, task::TaskX, arg...) = task_jacobian_chart(pm, task.task_map, arg...)

## Jacobian Time Derivatives
base_task_jacobian_emb_dot(xme, vme, taskx::Union{TaskX,TaskMapX}, arg...) =
    task_jacobian_emb_dot(xme, vme, base_task_map(taskx), arg...)
base_task_jacobian_chart_dot(xm, vm, taskx::Union{TaskX,TaskMapX}, arg...) =
    task_jacobian_chart_dot(xm, vm, base_task_map(taskx), arg...)

task_jacobian_emb_dot(xme, vme, task::TaskX, arg...) =
    task_jacobian_emb_dot(xme, vme, task.task_map, arg...)
task_jacobian_chart_dot(xm, vm, task::TaskX, arg...) =
    task_jacobian_chart_dot(xm, vm, task.task_map, arg...)

function task_jacobian_emb_dot(xme, vme, task_map::BaseTaskMap, arg...)
    m, n = embdim(domain_manifold(task_map)), embdim(codomain_manifold(task_map))
    ∇xf = ForwardDiff.jacobian(xme -> SVector{m*n,eltype(xme)}(reshape(task_jacobian_emb(xme,
        task_map, arg...), m*n)), xme)::SArray{Tuple{m*n,m},eltype(xme),2,m*m*n}
    Jfdot_emb = reshape(∇xf*vme, Size(n,m))
end

function task_jacobian_chart_dot(xm, vm, task_map::BaseTaskMap, arg...)
    m, n = dim(domain_manifold(task_map)), dim(codomain_manifold(task_map))
    ∇xf = ForwardDiff.jacobian(xm -> SVector{m*n,eltype(xm)}(reshape(task_jacobian_chart(xm,
        task_map, arg...), m*n)), xm)::SArray{Tuple{m*n,m},eltype(xm),2,m*m*n}
    Jfdot = reshape(∇xf*vm, Size(n,m))
end

# Default coordinate representations for tasks on Rn
domain_coord_rep(::TaskMap{ℝ{m},N,S}) where {m,N,S} = EmbRep()
codomain_coord_rep(::TaskMap{M,ℝ{n},S}) where {M,n,S} = EmbRep()

# Default home task charts for tasks on Rn
home_task_chart(::Task{<:TaskMap{M,ℝ{n},S}}) where {M,n,S} = Chart{1,ℝ{n}}()
home_task_chart(::Task{<:TaskMapT{M,ℝ{n},S}}) where {M,n,S} = Chart{1,ℝ{n}}()
home_task_chart(::TaskGDS{<:TaskMapT{M,ℝ{n},S}}) where {M,n,S} = Chart{1,ℝ{n}}()
home_task_chart(::Task{<:PTM{<:TaskMap{M1,ℝ{n1},S},<:TaskMap{M2,ℝ{n2},S}}}
    ) where {M1,M2,n1,n2,S} = Chart{Tuple{1,1},PM{ℝ{n1},ℝ{n2}}}()
home_task_chart(::Task{<:PTMT{<:PTM{<:TaskMap{M1,ℝ{n1},S},<:TaskMap{M2,ℝ{n2},S}}}}
    ) where {M1,M2,n1,n2,S} = Chart{Tuple{1,1},PM{ℝ{n1},ℝ{n2}}}()
# Default home task charts for S1
home_task_chart(::Task{<:TaskMap{M,S1,S}}) where {M,S} = Chart{Angleπ,S1}()
home_task_chart(::Task{<:PTM{<:TaskMap{M1,ℝ{n},S},<:TaskMap{M2,S1,S}}}
    ) where {M1,M2,n,S} = Chart{Tuple{1,Angleπ},PM{ℝ{n},S1}}()

## Task Weighting
task_weight(pn1, task::TaskX, C1) = 1. # Default weighting is 1 everywhere
task_weight(xn1, vn1, task::TaskX, C1) = 1.
