struct 𝕊{n} <: Manifold end

dim(::Type{𝕊{n}}) where n = n
embdim(::Type{𝕊{n}}) where n = n+1

# Creates aliases S1, TS1, etc.
for n = 1:3
    @eval const $(Symbol("S$n")) = 𝕊{$n}
    @eval const $(Symbol("TS$n")) = T{𝕊{$n}}
    @eval export $(Symbol("S$n"))
    @eval export $(Symbol("TS$n"))
end
TS{n} = T{𝕊{n}}

# Sterographic Projection
abstract type SterProj <: ChartID end
isglobal(::Chart{<:SterProj,𝕊{n}}) where n = false
struct SterProjSouth <: SterProj end  # South pole
struct SterProjNorth <: SterProj end  # North pole
emb_to_chart(pe, ::Chart{SterProjSouth,𝕊{n}}) where n = SVector{n}(pe[1:n])./(1+pe[n+1])
emb_to_chart(pe, ::Chart{SterProjNorth,𝕊{n}}) where n = SVector{n}(pe[1:n])./(1-pe[n+1])
function chart_to_emb(pn, ::Chart{SterProjSouth,S2})
    den = 1 + pn[1]^2 + pn[2]^2
    SA[2*pn[1]/den, 2*pn[2]/den, (2-den)/den]
end
function chart_to_emb(ps, ::Chart{SterProjNorth,S2})
    den = 1 + ps[1]^2 + ps[2]^2
    SA[2*ps[1]/den, 2*ps[2]/den, (den-2)/den]
end
chart_choice_coord_rep(::Chart{<:SterProj,S2}) = EmbRep()
choose_chart_emb(::EmbRep, pe, ::Chart{<:SterProj,S2}) =
    (pe[3] < 0)[1] ? Chart{SterProjSouth,S2}() : Chart{SterProjNorth,S2}()

# Same both ways
chart_transition(ps, ::Chart{SterProjSouth,𝕊{n}}, ::Chart{SterProjNorth,𝕊{n}}) where n =
    ps/norm(ps)^2
chart_transition(pn, ::Chart{SterProjNorth,𝕊{n}}, ::Chart{SterProjSouth,𝕊{n}}) where n =
    pn/norm(pn)^2

# Angle representation for S1
abstract type AngleChart <: ChartID end
isglobal(::Chart{<:AngleChart,S1}) = false
struct Angleπ <: AngleChart end    # Angles in (-π, π)
struct Angle2π <: AngleChart end   # Angles in (0, 2π)
default_chart(::Type{S1}) = Chart{Angleπ,S1}()
emb_to_chart(pe, ::Chart{Angleπ,S1}) = SA[atan(pe[2],pe[1])]
function emb_to_chart(pe, C2π::Chart{Angle2π,S1})
    Cπ = Chart{Angleπ,S1}()
    pπ = emb_to_chart(pe, Cπ) # Between -π and π
    chart_transition(pπ, Cπ, C2π)
end

chart_to_emb(p, ::Chart{<:AngleChart,S1}) = [cos.(p); sin.(p)]

chart_transition(p2π, ::Chart{Angle2π,S1}, ::Chart{Angleπ,S1}) = (p2π .< π)[1] ? p2π : p2π .- 2π
chart_transition(pπ, ::Chart{Angleπ,S1}, ::Chart{Angle2π,S1}) = (pπ .> 0.)[1] ? pπ : pπ .+ 2π

chart_choice_coord_rep(::Chart{<:AngleChart,S1}) = EmbRep()
choose_chart_emb(::EmbRep, pe, ::Chart{<:AngleChart,S1}) =
    (pe[1] .< 0.)[1] ? Chart{Angle2π,S1}() : Chart{Angleπ,S1}()