struct Identity{M,N,S} <: TaskMap{M,N,S} end
struct IdentityT{M,N,S} <: TaskMapT{M,N,S}
    base_map::Identity{M,N,S}
end
IdentityT{M,N,S}() where {M,N,S} = IdentityT{M,N,S}(Identity{M,N,S}())

task_map_emb(::EmbRep, ::EmbRep, xm, task_map::Identity{ℝ{m},ℝ{m}}) where m = xm

task_map_emb(::EmbRep, ::EmbRep, xme, task_map::Identity{𝕊{m},𝕊{m}}) where m = xme
domain_coord_rep(::Identity{𝕊{n},𝕊{n}}) where n = EmbRep()
codomain_coord_rep(::Identity{𝕊{n},𝕊{n}}) where n = EmbRep()