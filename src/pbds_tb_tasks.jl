function permute_chart(pn, ::Type{N}) where N <: Manifold
    n = dim(N)
    pn_perm = [pn[SVector{n}(n+1:2n)]; pn[SVector{n}(1:n)]]
end

function permuted_task_map_chart(pm, task_map::TaskMapT{M,N,S}, CM::Chart{I,M}, CN::Chart{J,N}) where {M,N,S,I,J}
    xm, vm = tangent_vector_chart_splitview(pm, T{M})
    xn = base_task_map_chart(xm, task_map, CM, CN)
    ∂xn_∂xm = base_task_jacobian_chart(xm, task_map, CM, CN)
    vn = ∂xn_∂xm*vm  # Chart representation cancels out effect of global frame
    pn_perm = [vn; xn]
end

function permuted_metric_chart(pn_perm, task::Task{<:TaskMapT{M,N,S}},
        C::Chart{I,N}) where {M,N,S,I}
    n = dim(N)
    pn = permute_chart(pn_perm, N)
    g = metric_chart(pn, task, C)
    indsx, indsv = SVector{n}(1:n), SVector{n}(n+1:2n)
    gxx = g[indsx,indsx]
    gxv = g[indsx,indsv]
    gvx = g[indsv,indsx]
    gvv = g[indsv,indsv]
    
    g_perm = SMatrix{2n,2n,eltype(pn_perm)}([gvv gvx; gxv gxx])
end

# indexing convention Γ^k_ij: (i, j, k)
function permuted_christoffel_symbols(pn_perm, task::Task{<:TaskMapT{M,N,S}},
        C::Chart{I,N}) where {M,N,I,S}
    n = dim(T{N})
    Ginv = inv(permuted_metric_chart(pn_perm, task, C))
    ∂G_∂xn_matrix = ForwardDiff.jacobian(
        pn_perm -> reshape(permuted_metric_chart(pn_perm, task, C), Size(n^2)), pn_perm)
    ∂G_∂xn = SArray{Tuple{n,n,n},S}(reshape(∂G_∂xn_matrix, Size(n,n,n)))
    inds = static(1):static(n)
    Γ = SArray{Tuple{n,n,n},S}([0.5*sum(Tuple(Ginv[k,l]*(∂G_∂xn[j,l,i]+∂G_∂xn[i,l,j]-∂G_∂xn[i,j,l])
        for l in inds)) for i in inds, j in inds, k in inds])
end

function permuted_task_jacobian_chart(pm, task_map::TaskMapT{M,N,S}, arg...)::
        SMatrix{dim(T{N}),dim(T{M}),S,dim(T{M})*dim(T{N})} where {M,N,S}
    ForwardDiff.jacobian(pm -> permuted_task_map_chart(pm, task_map, arg...), pm)
end

function permuted_potential_forces(pn, task::Task{<:TaskMapT{M,N,S}},
        CN::Chart{J,N}) where {M,N,S,J}
    n = dim(N)
    xn, vn = tangent_vector_chart_splitview(pn, T{N})
    ∂Φ_∂xn = ForwardDiff.gradient(xn -> potential(xn, task, CN), xn)
    [@SVector zeros(n); -∂Φ_∂xn]
end

function permuted_dissipative_forces(pn, task::Task{<:TaskMapT{M,N,S}},
        CN::Chart{J,N}) where {M,N,S,J}
    n = dim(N)
    xn, vn = tangent_vector_chart_splitview(pn, T{N})
    [@SVector zeros(n); dissipative_forces(xn, vn, task, CN)]
end

function single_task_components(xm, vm, task::Task{<:TaskMapT{M,N,S}}, CM::Chart{I,M},
        CN::Chart{J,N}) where {M,N,S,I,J}
    m, n = dim(M), dim(N)
    mt, nt = dim(T{M}), dim(T{N})
    
    σx, σxdot = xm, vm
    pm = [xm; vm]
    Jf = base_task_jacobian_chart(σx, task, CM, CN)
    Jfdot = base_task_jacobian_chart_dot(σx, σxdot, task, CM, CN)
    JF_perm = permuted_task_jacobian_chart(pm, task.task_map, CM, CN)
    pn = task_map_chart(pm, task, CM, CN)
    pn_perm = [pn[SVector{n}(n+1:nt)]; pn[SVector{n}(1:n)]]
    
    Γ_perm = permuted_christoffel_symbols(pn_perm, task, CN)
    g_perm = permuted_metric_chart(pn_perm, task, CN)
    inv_g_perm = inv(g_perm)
    inv_gv_perm = inv_g_perm[SVector{n}(n+1:nt),:]

    ℱ_pot = permuted_potential_forces(pn, task, CN)
    ℱ_dis = permuted_dissipative_forces(pn, task, CN)
    ℱ = ℱ_pot + ℱ_dis
    
    nt_inds, n_inds, m_inds = static(1):static(nt), static(1):static(n), static(1):static(m)
    Ξx = SMatrix{n,m,S}([sum(Tuple(JF_perm[l,j]*Γ_perm[l,h+n,k+n]*Jf[h,r]*σxdot[r]
        for l=nt_inds, h=n_inds, r=m_inds)) for k=n_inds, j=m_inds])
    Ξv = SMatrix{n,m,S}([sum(Tuple(JF_perm[l,j+m]*Γ_perm[l,h+n,k+n]*Jf[h,r]*σxdot[r]
        for l=nt_inds, h=n_inds, r=m_inds)) for k=n_inds, j=m_inds])
    λ = task_weight(pn, task, CN)

    JfΞv = Jf + Ξv
    JfΞvt_JfΞv = λ*JfΞv'*JfΞv
    JfΞvt_A = λ*JfΞv'*(inv_gv_perm*ℱ - (Jfdot + Ξx)*σxdot)
    
    JfΞvt_JfΞv, JfΞvt_A
end

function single_task_acceleration(xm, vm, task::Task{<:TaskMapT{M,N,S}}, CM::Chart{I,M},
        CN::Chart{J,N}) where {M,N,S,I,J}
    m, n = dim(M), dim(N)
    CN = isglobal(CN) ? CN : choose_chart_emb(xm, task, CM, CN)
    Jft_JfΞv, Jft_A = single_task_components(xm, vm, task, CM, CN)
    σxddot = SMatrix{m,m,S}(pinv(Matrix(Jft_JfΞv)))*Jft_A
    σxddot, CN
end

function multiple_task_acceleration(xm, vm, tasks::TaskTList, CM::Chart{I,M}, CNs::ChartList;
        log_task_chart=false) where {M,I}
    m = dim(M)
    Jft_JfΞv_sum = zeros(m,m)
    Jft_A_sum = zeros(m)
    CNs_out = ChartList()
    for i in 1:length(tasks)
        CN = isglobal(CNs[i]) ? CNs[i] : choose_chart_emb(xm, tasks[i], CM, CNs[i])
        Jft_JfΞv, Jft_A = single_task_components(xm, vm, tasks[i], CM, CN)
        Jft_JfΞv_sum += Jft_JfΞv
        Jft_A_sum += Jft_A
        log_task_chart && push!(CNs_out, CN)
    end
    σxddot = SMatrix{m,m,eltype(xm)}(pinv(Matrix(Jft_JfΞv_sum)))*Jft_A_sum
    σxddot, CNs_out
end