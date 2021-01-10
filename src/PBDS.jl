module PBDS
# __precompile__(false)

using ComputedFieldTypes
using LinearAlgebra, ForwardDiff
using StaticArrays, StaticNumbers, Tullio
using RigidBodyDynamics, MeshCatMechanisms
using GeometryBasics, Rotations
using Base64

import RigidBodyDynamics: CustomCollections, CacheElement, isdirty, setdirty!
import GeometryBasics: min_euclidean, min_euclideansq, min_dist_dim

# types.jl
export TangentBundle, T, ProductManifold, PM, ProductTaskMap, PTM, Chart, ChartList, TaskList,
       TaskGDSList, EmbRep, ChartRep
# manifolds.jl
export dim, embdim, basedim, default_chart, choose_chart_emb, choose_chart_chart, base,
       emb_to_chart, chart_to_emb
# manifolds/
export ℝ, 𝕊, SterProj, SterProjSouth, SterProjNorth, Angleπ, Angle2π
# tasks.jl
export task_map_emb, task_map_emb_chart, task_map_chart_emb, task_map_chart,
       base_task_map_emb, base_task_map_emb_chart, base_task_map_chart_emb, base_task_map_chart,
       task_jacobian_emb, task_jacobian_emb_chart, task_jacobian_chart_emb, task_jacobian_chart,
       base_task_jacobian_emb, base_task_jacobian_emb_chart,
       base_task_jacobian_chart_emb, base_task_jacobian_chart,
       task_jacobian_emb_dot, task_jacobian_chart_dot,
       base_task_jacobian_emb_dot, base_task_jacobian_chart_dot
# tasks/
export AngularDistance, AngularDistanceT,
       AngularPositionAroundPoint, AngularPositionAroundPointT,
       DistanceFromPoint, DistanceFromPointT, DistanceFromPointDynamic,
       DistanceFromSphereSurface, DistanceFromSphereSurfaceT,
       DistanceSphereToSphere, DistanceSphereToSphereT,
       DistanceSphereToCylinder, DistanceSphereToCylinderT,
       DistanceSphereToCup, DistanceSphereToCupT, DistanceSphereToCupDynamic,
       DistanceSphereToBox, DistanceSphereToBoxT, DistanceSphereToBoxDynamic,
       DistanceSphereToLine, DistanceSphereToLineT,
       Identity, IdentityT,
       PositionAroundPoint, PositionAroundPointT,
       Coordinate, CoordinateT,
       XCoordinateDistance, XCoordinateDistanceT,
       YCoordinateDistance, YCoordinateDistanceT,
       LinkFrameTransform, LinkFrameTransformT,
       LinkSpherePosition, LinkSpherePositionT
# tasks/arm/
export @ext_standard_task_type, LinkSpherePairDistance, LinkSpherePairDistanceT,
       JointToLinkPosition, JointToLinkPositionT
# tasks/task_types.jl
export BlankTask, Damping, Attractor, CollisionAvoidance, VelocityLimiter,
       BlankTaskGDS, DampingGDS, AttractorGDS, CollisionAvoidanceGDS, VelocityLimiterGDS
# caches.jl
export ControllerCache
# caches/
export FramePointCache, LinkFrameCache
# trees.jl
export TreeNode, add_child!, check_tree, set_charts_emb!, set_charts_chart!
# forces.jl
export default_potential, default_dissipative_forces,
       potential_from_home_chart, dissipative_forces_from_home_chart
# metrics.jl
export default_metric, metric_from_home_chart, default_weight_metric, weight_metric_from_home_chart
# zero_order_tasks.jl
export single_task_acceleration, multiple_task_acceleration, task_acceleration
# horizontal_bundle_tasks.jl
export dissipative_term_emb, dissipative_term_chart
# trajectory.jl
export propagate_task, propagate_tasks, propagate_task_emb, propagate_tasks_emb
# utils.jl
export smooth_normalization, smooth_norm, exp_barrier_value, exp_barrier_slope, softplus_slope, 
       smooth_max, smooth_min, smooth_clamp, smooth_abs, logistic_val, wrapTo2Pi, wrapToPi,
       line_segment_distance, aligned_box_distance, html_video

include("types.jl")
include("manifolds.jl")
include("metrics.jl")
include("connections.jl")
include("caches.jl")
include("tasks.jl")
include("forces.jl")
include("trees.jl")
include("pbds_tasks.jl")
include("pbds_tb_tasks.jl")
include("rmpflow_tasks.jl")
include("trajectory.jl")
include("utils.jl")

end #module