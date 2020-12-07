struct DistanceSphereToCup{R3,R1,S} <: TaskMap{R3,R1,S}
    center_bottom::SVector{3,S}
    axis::SVector{3,S}
    cup_radius::S
    sph_radius
    height::S
end
function DistanceSphereToCup{R3,R1}(cb::SVector{3,S}, a::SVector{3,S}, cr::S,
        cs::S, h::S) where {R3,R1,S}
    DistanceSphereToCup{R3,R1,S}(cb, a./norm(a), cr, cs, h)
end

struct DistanceSphereToCupT{R3,R1,S} <: TaskMapT{R3,R1,S}
    base_map::DistanceSphereToCup{R3,R1,S}
end
function DistanceSphereToCupT{R3,R1}(cb::SVector{3,S}, a::SVector{3,S}, cr::S,
        cs::S, h::S) where {R3,R1,S}
    DistanceSphereToCupT{R3,R1,S,m}(DistanceSphereToCup{R3,R1}(cb, a, cr, cs, h))
end

function task_map_emb(::EmbRep, ::EmbRep, xme, task_map::DistanceSphereToCup{R3,R1})
    center_bottom, axis = task_map.center_bottom, task_map.axis
    cup_radius, sph_radius, height = task_map.cup_radius, task_map.sph_radius, task_map.height
    α = 1.e-20                                  # Smoothing factor
    p_world = xme - center_bottom               # Position from origin in world frame
    y_dist_cup = dot(axis, p_world)             # Axial distance from origin
    y_world = y_dist_cup*axis                   # Axial position in world frame
    x_world = p_world - y_world                 # Transverse position in world frame
    x_dist_cup = smooth_norm(x_world, α)        # Transverse distance from origin
    p_cup = SA[x_dist_cup, y_dist_cup]

    p_corner1 = SA[-cup_radius, 0.]
    p_corner2 = SA[cup_radius, 0.]
    p_corner3 = SA[cup_radius, height]
    d1 = line_segment_distance(p_cup, p_corner1, p_corner2)
    d2 = line_segment_distance(p_cup, p_corner2, p_corner3)

    SA[min(max(0., d1 - sph_radius), max(0, d2 - sph_radius))]
end

# Updateable Cup Pose
struct DistanceSphereToCupDynamic{R3,R1,S} <: TaskMap{R3,R1,S}
    center_bottom::Vector{SVector{3,S}}
    axis::Vector{SVector{3,S}}
    cup_radius::S
    sph_radius
    height::S
end
function DistanceSphereToCupDynamic{R3,R1}(cb::Vector{SVector{3,S}}, a::Vector{SVector{3,S}}, cr::S,
        cs::S, h::S) where {R3,R1,S}
    DistanceSphereToCupDynamic{R3,R1,S}(cb, [a[1]./norm(a[1])], cr, cs, h)
end

struct DistanceSphereToCupDynamicT{R3,R1,S} <: TaskMapT{R3,R1,S}
    base_map::DistanceSphereToCupDynamic{R3,R1,S}
end
function DistanceSphereToCupDynamicT{R3,R1}(cb::SVector{3,S}, a::SVector{3,S}, cr::S,
        cs::S, h::S) where {R3,R1,S}
    DistanceSphereToCupDynamicT{R3,R1,S,m}(DistanceSphereToCupDynamic{R3,R1}(cb, a, cr, cs, h))
end

function task_map_emb(::EmbRep, ::EmbRep, xme, task_map::DistanceSphereToCupDynamic{R3,R1})
    center_bottom, axis = task_map.center_bottom[1], task_map.axis[1]
    cup_radius, sph_radius, height = task_map.cup_radius, task_map.sph_radius, task_map.height
    α = 1.e-20                                  # Smoothing factor
    p_world = xme - center_bottom               # Position from origin in world frame
    y_dist_cup = dot(axis, p_world)             # Axial distance from origin
    y_world = y_dist_cup*axis                   # Axial position in world frame
    x_world = p_world - y_world                 # Transverse position in world frame
    x_dist_cup = smooth_norm(x_world, α)        # Transverse distance from origin
    p_cup = SA[x_dist_cup, y_dist_cup]

    p_corner1 = SA[-cup_radius, 0.]
    p_corner2 = SA[cup_radius, 0.]
    p_corner3 = SA[cup_radius, height]
    d1 = line_segment_distance(p_cup, p_corner1, p_corner2)
    d2 = line_segment_distance(p_cup, p_corner2, p_corner3)

    SA[min(max(0., d1 - sph_radius), max(0, d2 - sph_radius))]
end
