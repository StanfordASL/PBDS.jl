product_manifold_split(::Type{PM{M1,M2}}) where {M1,M2} = M1, M2 

dim(::Union{PM{M1,M2},Type{PM{M1,M2}}}) where {M1,M2} = dim(M1) + dim(M2)
embdim(::Union{PM{M1,M2},Type{PM{M1,M2}}}) where {M1,M2} = embdim(M1) + embdim(M2)

dims(::Union{PM{M1,M2},Type{PM{M1,M2}}}) where {M1,M2} = dim(M1), dim(M2)
embdims(::Union{PM{M1,M2},Type{PM{M1,M2}}}) where {M1,M2} = embdim(M1), embdim(M2)

product_chart_split(C::Chart{Tuple{I,J},PM{M1,M2}}) where {M1,M2,I,J} =
    Chart{I,M1}(), Chart{J,M2}()
product_chart_merge(C1::Chart{I,M1}, C2::Chart{J,M2}) where {M1,M2,I,J} =
    Chart{Tuple{I,J},PM{M1,M2}}()

isglobal(::Chart{Tuple{I,J},PM{M1,M2}}) where {M1,M2,I,J} =
    isglobal(Chart{I,M1}()) && isglobal(Chart{J,M2}())

function product_manifold_emb_splitview(pe, ::Type{M}) where M <: ProductManifold
    m = embdim(M)
    m1, m2 = embdims(M)
    x1e, x2e = pe[SVector{m1}(1:m1)], pe[SVector{m2}(m1+1:m)]
end

function product_manifold_chart_splitview(p, ::Type{M}) where M <: ProductManifold
    m = dim(M)
    m1, m2 = dims(M)
    x1, x2 = p[SVector{m1}(1:m1)], p[SVector{m2}(m1+1:m)]
end

function emb_to_chart(xe, C::Chart{I,M}) where {M<:PM,I}
    C1, C2 = product_chart_split(C)
    x1e, x2e = product_manifold_emb_splitview(xe, M)
    x1, x2 = emb_to_chart(x1e, C1), emb_to_chart(x2e, C2)
    x = [x1; x2]
end

function chart_transition(xa, Ca::Chart{Tuple{Ia,Ja},M},
        Cb::Chart{Tuple{Ib,Jb},M}) where {M<:PM,Ia,Ja,Ib,Jb}
    C1a, C2a = product_chart_split(Ca)
    C1b, C2b = product_chart_split(Cb)
    x1a, x2a = product_manifold_chart_splitview(xa, M)
    x1b, x2b = chart_transition(x1a, C1a, C1b), chart_transition(x2a, C2a, C2b)
    xb = [x1b; x2b]
end

function choose_chart_emb(xe, C::Chart{Tuple{I,J},M}) where {M<:PM,I,J}
    C1, C2 = product_chart_split(C)
    x1e, x2e = product_manifold_emb_splitview(xe, M)
    C1c, C2c = choose_chart_emb(x1e, C1), choose_chart_emb(x2e, C2)
    product_chart_merge(C1c, C2c)
end

function choose_chart_chart(x, C::Chart{Tuple{I,J},M}) where {M<:PM,I,J}
    C1, C2 = product_chart_split(C)
    x1, x2 = product_manifold_chart_splitview(x, M)
    C1c, C2c = choose_chart_chart(x1, C1), choose_chart_chart(x2, C2)
    product_chart_merge(C1c, C2c)
end