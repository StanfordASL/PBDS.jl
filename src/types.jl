abstract type Manifold end
struct TangentBundle{M<:Manifold} <: Manifold end
const T = TangentBundle
struct ProductManifold{M1<:Manifold,M2<:Manifold} <: Manifold end
const PM = ProductManifold

abstract type ChartID end
# Usually for regular manifolds, I <: ChartID or type(I) = Int
# For product manifolds, I <: Tuple{Int,Int}
struct Chart{I, M<:Manifold} end
ChartList = Vector{Chart}

abstract type TaskMapX{A,B,S} end
abstract type TaskMap{M<:Manifold, N<:Manifold, S} <: TaskMapX{M,N,S} end
abstract type TaskMapT{M<:Manifold, N<:Manifold, S} <: TaskMapX{M,N,S} end
Base.eltype(::TaskMapX{A,B,S}) where {A,B,S} = S

abstract type TaskX{F} end
abstract type Task{F<:TaskMapX} <: TaskX{F} end
abstract type TaskGDS{F<:TaskMapT} <: TaskX{F} end
Base.eltype(::Task{F}) where F = eltype(F)

struct ProductTaskMap{T1<:TaskMap,T2<:TaskMap,S} <: TaskMapX{T1,T2,S}
    map1::T1
    map2::T2
end
const PTM = ProductTaskMap
PTM(map1::T1,map2::T2) where {T1,T2} = PTM{T1,T2,eltype(T1)}(map1,map2)
Base.eltype(::Type{PTM{T1,T2,S}}) where {T1,T2,S} = S
BaseTaskMap = Union{TaskMap, ProductTaskMap}

struct ProductTaskMapT{T1<:TaskMap,T2<:TaskMap,S} <: TaskMapX{T1,T2,S}
    base_map::ProductTaskMap{T1,T2,S}
end
ProductTaskMapT(T1,T2) = ProductTaskMapT(ProductTaskMap(T1,T2))
const PTMT = ProductTaskMapT

ProductTaskMapX = Union{PTM, PTMT}
PTMX = ProductTaskMapX

TaskList = Vector{Task{<:BaseTaskMap}}
TaskTList = Vector{Task{<:TaskMapT}}
TaskGDSList = Vector{TaskGDS{<:TaskMapT}}
TaskXList = Union{TaskList, TaskTList, TaskGDSList}

# Traits for specifying base definitions of task maps
# Also for specifying coordinate representation used for manifold chart selection,
# task trajectory logging
abstract type CoordinateRep end
struct EmbRep <: CoordinateRep end
struct ChartRep <: CoordinateRep end