struct ℝ{n} <: Manifold end

dim(::Type{ℝ{n}}) where n = n

# Creates aliases R1, TR1, etc.
for n = 1:12
    @eval const $(Symbol("R$n")) = ℝ{$n}
    @eval const $(Symbol("TR$n")) = T{ℝ{$n}}
    @eval export $(Symbol("R$n"))
    @eval export $(Symbol("TR$n"))
end
TR{n} = T{ℝ{n}}

# Assume all charts are global for now
isglobal(::Chart{k,ℝ{n}}) where {k,n} = true

emb_to_chart(p, C::Chart{1,ℝ{n}}) where n = p
chart_to_emb(p, C::Chart{1,ℝ{n}}) where n = p

# Not a diffeomorphism for all of R^n, but useful for manual checks
emb_to_chart(p, C::Chart{2,ℝ{n}}) where n = p.^2
chart_to_emb(p, C::Chart{2,ℝ{n}}) where n = p.^(1/2)

# Not a diffeomorphism for all of R^n, but useful for manual checks
emb_to_chart(p, C::Chart{3,ℝ{n}}) where n = p.^3
chart_to_emb(p, C::Chart{3,ℝ{n}}) where n = p.^(1/3)

emb_to_chart(p, C::Chart{4,ℝ{n}}) where n = p.^3 + p
function chart_to_emb(p, C::Chart{4,ℝ{n}}) where n
    a = sqrt(3)*sqrt.(27*p.^2 .+ 4) .+ 9*p
    (a/18).^(1/3) .- ((2/3)./a).^(1/3)
end

global_frame_coeffs(p, C::Chart{1,ℝ{n}}) where n = SMatrix{n,n,eltype(p)}(I)
global_coframe_coeffs(p, C::Chart{1,ℝ{n}}) where n = SMatrix{n,n,eltype(p)}(I)

function global_frame_coeffs(p, C::Chart{I,ℝ{n}}) where {n,I}
    C1 = Chart{1,ℝ{n}}()
    p1 = chart_transition(p, C, C1)
    B1 = global_frame_coeffs(p, C1)
    global_frame_coeffs_chart_transition(p1, B1, C1, C)
end

function global_coframe_coeffs(p, C::Chart{I,ℝ{n}}) where {n,I}
    C1 = Chart{1,ℝ{n}}()
    p1 = chart_transition(p, C, C1)
    β1 = global_frame_coeffs(p, C1)
    global_coframe_coeffs_chart_transition(p1, β1, C1, C)
end

inv_global_frame_coeffs(p, C) = inv(global_coframe_coeffs(p, C))
inv_global_coframe_coeffs(p, C) = inv(global_frame_coeffs(p, C))