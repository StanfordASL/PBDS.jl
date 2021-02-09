struct Coordinate{M,R1,S} <: TaskMap{M,R1,S}
    ind::Int
end
struct CoordinateT{M,R1,S} <: TaskMapT{M,R1,S}
    base_map::Coordinate{M,R1,S}
end
CoordinateT{M,R1,S}(i::Int) where {M,R1,S} = CoordinateT{M,R1,S}(Coordinate{M,R1,S}(i))

task_map_emb(::EmbRep, ::EmbRep, xm, task_map::Coordinate{ℝ{m}, R1}) where m =
    SA[xm[task_map.ind]]