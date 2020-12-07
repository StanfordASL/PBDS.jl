struct XCoordinateDistance{M,R1,S} <: TaskMap{M,R1,S}
    x_position::S
end
XCoordinateDistance{M,R1}(c::S) where {M,R1,S} = XCoordinateDistance{M,R1,S}(c)

struct XCoordinateDistanceT{M,R1,S} <: TaskMapT{M,R1,S}
    base_map::XCoordinateDistance{M,R1,S}
end
XCoordinateDistanceT{M,R1}(c::S) where {M,R1,S} =
    XCoordinateDistanceT{M,R1,S}(XCoordinateDistance{M,R1,S}(c))

task_map_emb(::EmbRep, ::EmbRep, xm, task_map::XCoordinateDistance{ℝ{m}, R1}) where m =
    SA[abs(xm[1] - task_map.x_position)]

struct YCoordinateDistance{M,R1,S} <: TaskMap{M,R1,S}
    y_position::S
end
YCoordinateDistance{M,R1}(c::S) where {M,R1,S} = YCoordinateDistance{M,R1,S}(c)

struct YCoordinateDistanceT{M,R1,S} <: TaskMapT{M,R1,S}
    base_map::YCoordinateDistance{M,R1,S}
end
YCoordinateDistanceT{M,R1}(c::S) where {M,R1,S} =
    YCoordinateDistanceT{M,R1,S}(YCoordinateDistance{M,R1,S}(c))

task_map_emb(::EmbRep, ::EmbRep, xm, task_map::YCoordinateDistance{ℝ{m}, R1}) where m =
    SA[abs(xm[2] .- task_map.y_position)]