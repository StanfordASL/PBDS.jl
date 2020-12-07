struct Identity{M,N,S} <: TaskMap{M,N,S} end
struct IdentityT{M,N,S} <: TaskMapT{M,N,S}
    base_map::Identity{M,N,S}
end
IdentityT{M,N,S}() where {M,N,S} = IdentityT{M,N,S}(Identity{M,N,S}())

task_map_emb(::EmbRep, ::EmbRep, xm, task_map::Identity{ℝ{m},ℝ{m},S}) where {m,S} = xm

task_map_emb(::EmbRep, ::EmbRep, xme, task_map::Identity{𝕊{n},𝕊{n},S}) where {n,S} = xme
domain_coord_rep(::Identity{𝕊{n},𝕊{n},S}) where {n,S} = EmbRep()
codomain_coord_rep(::Identity{𝕊{n},𝕊{n},S}) where {n,S} = EmbRep()