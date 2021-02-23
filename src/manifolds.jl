include(joinpath("manifolds", "SN.jl"))
include(joinpath("manifolds", "SO3.jl"))
include(joinpath("manifolds", "RN.jl"))
include(joinpath("manifolds", "product_manifolds.jl"))

embdim(::Type{M}) where M <: Manifold = dim(M)
basedim(::Type{M}) where M <: Manifold = dim(M)

dim(::M) where M <: Manifold = dim(M)
embdim(::M) where M <: Manifold = embdim(M)
basedim(::M) where M <: Manifold = dim(M)

isglobal(::Chart) = true    # Charts are considered global by default
default_chart(::Type{M}) where M <: Manifold = Chart{1,M}()
default_chart(::M) where M <: Manifold = default_chart(M)

function chart_transition(p1, C1::Chart{I,M}, C2::Chart{J,M}) where {M,I,J}
    p_emb = chart_to_emb(p1, C1)
    p2 = emb_to_chart(p_emb, C2)
end

choose_chart_emb(pe, C::Chart) =
    isglobal(C) ? C : choose_chart_emb(chart_choice_coord_rep(C), pe, C)
choose_chart_chart(p, C::Chart) =
    isglobal(C) ? C : choose_chart_chart(chart_choice_coord_rep(C), p, C)

choose_chart_emb(::ChartRep, pe, C) = choose_chart_chart(ChartRep(), emb_to_chart(pe, C), C)
choose_chart_chart(::EmbRep, p, C) = choose_chart_emb(EmbRep(), chart_to_emb(p, C), C)

choose_chart(::EmbRep, pe, C) = choose_chart_emb(p, C)
choose_chart(::ChartRep, p, C) = choose_chart_chart(p, C)

### Tangent bundles
T(::Chart{I,M}) where {M<:Manifold,I} = Chart{I,T{M}}()
base(::Chart{I,T{M}}) where {M<:Manifold,I} = Chart{I,M}()

dim(::Type{TangentBundle{M}}) where M <: Manifold = 2*dim(M)
embdim(::Type{TangentBundle{M}}) where M <: Manifold = 2*embdim(M)
basedim(::Type{TangentBundle{M}}) where M <: Manifold = dim(M)
baseembdim(::Type{TangentBundle{M}}) where M <: Manifold = embdim(M)

function tangent_vector_emb_splitview(pe, ::Type{TM}) where TM <: TangentBundle
    m, m_tb = baseembdim(TM), embdim(TM)
    xe, ve = pe[SVector{m}(1:m)], pe[SVector{m}(m+1:m_tb)]
end

function tangent_vector_chart_splitview(p, ::Type{TM}) where TM <: TangentBundle
    m, m_tb = basedim(TM), dim(TM)
    x, v = p[SVector{m}(1:m)], p[SVector{m}(m+1:m_tb)]
end

function chart_to_emb_differential(x, v, C::Chart{I,M}) where {M,I}
    xe = chart_to_emb(x, C)
    ∂xe_∂x = chart_to_emb_jacobian(x, C)
    ve = ∂xe_∂x*v
    xe, ve
end

function chart_to_emb_differential(x, v, a, C::Chart{I,M}) where {M,I}
    xe = chart_to_emb(x, C)
    p = [x; v]
    w = [v; a]
    ∂pe_∂p = chart_to_emb_jacobian(p, T(C))
    we = ∂pe_∂p*w
    ve, ae = tangent_vector_emb_splitview(we, T{M})
    xe, ve, ae
end

function emb_to_chart_differential(xe, ve, C::Chart{I,M}) where {M,I}
    x = emb_to_chart(xe, C)
    ∂x_∂xe = emb_to_chart_jacobian(xe, C)
    v = ∂x_∂xe*ve
    x, v
end

function emb_to_chart_differential(xe, ve, ae, C::Chart{I,M}) where {M,I}
    x = emb_to_chart(xe, C)
    pe = [xe; ve]
    we = [ve; ae]
    ∂p_∂pe = emb_to_chart_jacobian(pe, T(C))
    w = ∂p_∂pe*we
    v, a = tangent_vector_chart_splitview(w, T{M})
    x, v, a
end

function emb_to_chart(pe, C::Chart{I,T{M}}) where {M,I}
    xe, ve = tangent_vector_emb_splitview(pe, T{M})
    x, v = emb_to_chart_differential(xe, ve, base(C))
    p = [x; v]
end

function chart_to_emb(p, C::Chart{I,T{M}}) where {M,I}
    x, v = tangent_vector_chart_splitview(p, T{M})
    xe, ve = chart_to_emb_differential(x, v, Chart{I,M}())
    [xe; ve]
end

emb_to_chart_jacobian(pe, C::Chart{I,M}) where {M,I} =
    ForwardDiff.jacobian(pe -> emb_to_chart(pe, C), pe)
chart_to_emb_jacobian(pe, C::Chart{I,M}) where {M,I} =
    ForwardDiff.jacobian(pe -> chart_to_emb(pe, C), pe)
chart_transition_jacobian(p1, C1::Chart{I,M}, C2::Chart{J,M}) where {M,I,J} = 
    ForwardDiff.jacobian(p1 -> chart_transition(p1, C1, C2), p1)

function chart_transition_jacobian(p1, C1::Chart{I,T{M}}, C2::Chart{J,T{M}}) where {M,I,J}
    m = dim(M)
    x1, v1 = tangent_vector_chart_splitview(p1, T{M})
    ∂x2_∂x1 = chart_transition_jacobian(x1, base(C1), base(C2))
    Hx2_x1 =  chart_transition_hessian(x1, base(C1), base(C2))
    ∂v2_∂v1 = SMatrix{dim(T{M}),dim(T{M})}([∂x2_∂x1 @SMatrix zeros(m,m)
               sum(Hx2_x1[k,:,:]*v1[k] for k in 1:m) ∂x2_∂x1])
end

function chart_transition_differential(x1, v1, C1::Chart{I,M}, C2::Chart{J,M}) where {M,I,J}
    x2 = chart_transition(x1, C1, C2)
    ∂x2_∂x1 = chart_transition_jacobian(x1, C1, C2)
    v2 = ∂x2_∂x1*v1
    x2, v2
end

function chart_transition_differential(x1, v1, a1, C1::Chart{I,M}, C2::Chart{J,M}) where {M,I,J}
    x2 = chart_transition(x1, C1, C2)
    p1 = [x1; v1]
    w1 = [v1; a1]
    ∂p2_∂p1 = chart_transition_jacobian(p1, T(C1), T(C2))
    w2 = ∂p2_∂p1*w1
    v2, a2 = tangent_vector_chart_splitview(w2, T{M})
    x2, v2, a2
end

function chart_transition(p1, C1::Chart{I,T{M}}, C2::Chart{J,T{M}}) where {M,I,J}
    x1, v1 = tangent_vector_chart_splitview(p1, T{M})
    x2, v2 = chart_transition_differential(x1, v1, base(C1), base(C2))
    p2 = [x2; v2]
end

# Print H[i,:,:] to see symmetric matrices
# index convention: ∂p2k_∂p1i∂p1j -> H[k,i,j]
function chart_transition_hessian(p1, C1::Chart{I,M}, C2::Chart{J,M}) where {M,I,J}
    m, S = dim(M), eltype(p1)
    inds = static(1):static(m)
    SArray{Tuple{m,m,m},S}([ForwardDiff.hessian(p1 -> chart_transition(p1, C1, C2)[k], p1)[i,j]
        for k in inds, i in inds, j in inds])
end

# Frames and coframes
function global_frame_coeffs_chart_transition(p1, B1, C1::Chart{I,M}, C2::Chart{J,M}) where {M,I,J}
    ∂p2_∂p1 = chart_transition_jacobian(p1, C1, C2)
    B2 = B1*∂p2_∂p1'
end

function global_coframe_coeffs_chart_transition(p1, β1, C1::Chart{I,M},
        C2::Chart{J,M}) where {M,I,J}
    ∂p2_∂p1 = chart_transition_jacobian(p1, C1, C2)
    β2 = inv(∂p2_∂p1)'*β1
end
