# indexing convention for ∂G_ij/∂p^k: (i, j, k)
# indexing convention Γ^k_ij: (i, j, k)
function christoffel_symbols(pn, task::Task{<:TaskMapT{M,N}}, C::Chart{I,N}) where {M,N,I}
    n, S = dim(T{N}), eltype(pn)
    Ginv = inv(metric_chart(pn, task, C))
    ∂G_∂xn_matrix = ForwardDiff.jacobian(pn -> reshape(metric_chart(pn, task, C), Size(n^2)), pn)
    ∂G_∂xn_matrix == (@SMatrix zeros(n^2,n)) && return @SArray zeros(n,n,n)
    ∂G_∂xn = reshape(∂G_∂xn_matrix, Size(n,n,n))
    inds = static(1):static(n)
    Γ = SArray{Tuple{n,n,n},}([0.5*sum(Tuple(Ginv[k,l]*(∂G_∂xn[j,l,i]+∂G_∂xn[i,l,j]-∂G_∂xn[i,j,l])
        for l in inds)) for i in inds, j in inds, k in inds])
end

function christoffel_symbols(pn, task::Task{<:BaseTaskMap}, CN)
    n, S = dim(codomain_manifold(task)), eltype(pn)
    Ginv = inv(metric_chart(pn, task, CN))
    ∂G_∂xn_matrix = ForwardDiff.jacobian(pn -> reshape(metric_chart(pn, task, CN), Size(n^2)), pn)
    ∂G_∂xn_matrix == (@SMatrix zeros(n^2,n)) && return @SArray zeros(n,n,n)
    ∂G_∂xn = reshape(∂G_∂xn_matrix, Size(n,n,n))
    inds = static(1):static(n)
    Γ = SArray{Tuple{n,n,n},S}([0.5*sum(Tuple(Ginv[k,l]*(∂G_∂xn[j,l,i]+∂G_∂xn[i,l,j]-∂G_∂xn[i,j,l])
        for l in inds)) for i in inds, j in inds, k in inds])
end

function christoffel_symbols_chart_transition(pn1, Γ1, C1::Chart{I,N},
        C2::Chart{J,N}) where {N,I,J}
    n, S = dim(N), eltype(pn1)
    Γ1 == (@SArray zeros(n,n,n)) && return @SArray zeros(n,n,n)
    pn2 = chart_transition(pn1, C1, C2)
    ∂pn1_∂pn2 = chart_transition_jacobian(pn2, C2, C1)
    ∂pn2_∂pn1 = inv(∂pn1_∂pn2)
    Hpn2_pn1 = chart_transition_hessian(pn2, C2, C1)
    inds = static(1):static(n)
    Γ2 = SArray{Tuple{n,n,n},S}([sum(Tuple(∂pn2_∂pn1[k,l]*∂pn1_∂pn2[r,i]*∂pn1_∂pn2[s,j]*Γ1[r,s,l] 
        for l in inds, r in inds, s in inds)) + sum(Hpn2_pn1[l,i,j]*∂pn2_∂pn1[k,l] for l in inds) 
        for i in inds, j in inds, k in inds])
end