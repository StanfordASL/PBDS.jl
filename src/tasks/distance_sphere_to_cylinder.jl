struct DistanceSphereToCylinder{R3,R1,S} <: TaskMap{R3,R1,S}
    center_bottom::SVector{3,S}
    axis::SVector{3,S}
    cyl_radius::S
    sph_radius
    height::S
end
function DistanceSphereToCylinder{R3,R1}(cb::SVector{3,S}, a::SVector{3,S}, cr::S,
        cs::S, h::S) where {R3,R1,S}
    DistanceSphereToCylinder{R3,R1,S}(cb, a./norm(a), cr, cs, h)
end

struct DistanceSphereToCylinderT{R3,R1,S} <: TaskMapT{R3,R1,S}
    base_map::DistanceSphereToCylinder{R3,R1,S}
end
function DistanceSphereToCylinderT{R3,R1}(cb::SVector{3,S}, a::SVector{3,S}, cr::S,
        cs::S, h::S) where {R3,R1,S}
    DistanceSphereToCylinderT{R3,R1,S,m}(DistanceSphereToCylinder{R3,R1}(cb, a, cr, cs, h))
end
    
function task_map_emb(::EmbRep, ::EmbRep, xme, task_map::DistanceSphereToCylinder{R3,R1})
    center_bottom, axis = task_map.center_bottom, task_map.axis
    cyl_radius, sph_radius, height = task_map.cyl_radius, task_map.sph_radius, task_map.height
    α = 1.e-20                                  # Smoothing factor
    p_world = xme - center_bottom               # Position from origin in world frame
    y_dist_cyl = dot(axis, p_world)             # Axial distance from origin
    y_world = y_dist_cyl*axis                   # Axial position in world frame
    x_world = p_world - y_world                 # Transverse position in world frame
    x_dist_cyl = smooth_norm(x_world, α)        # Transverse distance from origin
    p_cyl = SA[x_dist_cyl, y_dist_cyl]

    p_corner1 = SA[cyl_radius, 0.]
    p_corner2 = SA[cyl_radius, height]
    d = line_segment_distance(p_cyl, p_corner1, p_corner2)
    SA[max(0., d - sph_radius)]
end