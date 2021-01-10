function single_task_components(xm, vm, task::Task{<:BaseTaskMap}, CM::Chart{I,M},
        CN::Chart{J,N}) where {M,N,I,J}
    m, n = dim(M), dim(N)
    
    σ, σdot = xm, vm
    Jf = task_jacobian_chart(σ, task, CM, CN)
    Jfdot = task_jacobian_chart_dot(σ, σdot, task, CM, CN)
    xn = task_map_chart(σ, task, CM, CN)
    vn = Jf*vm
    
    Γ = christoffel_symbols(xn, task, CN)
    if any(isinf.(Γ)) || any(isnan.(Γ))
        Γ = eltype(xm).(christoffel_symbols(BigFloat.(xn), task, CN))
    end
    g = metric_chart(xn, task, CN)
    ginv = inv(g)

    ℱ_pot = potential_force_chart(xn, task, CN)
    ℱ_dis = dissipative_forces_chart(xn, vn, task, CN)
    ℱ = ℱ_pot + ℱ_dis
    
    m_inds, n_inds = static(1):static(m), static(1):static(n)
    Ξ = SMatrix{n,m,eltype(xm)}([sum(Tuple(Jf[l,j]*Γ[l,h,k]*Jf[h,r]*σdot[r]
        for l=n_inds, h=n_inds, r=m_inds)) for k=n_inds, j=m_inds])
    W = weight_metric_chart(xn, vn, task, CN)

    JftWJf = Jf'*W*Jf
    JftWA = Jf'*W*(ginv*ℱ - (Jfdot + Ξ)*σdot)

    JftWJf, JftWA
end

function single_task_acceleration(xm, vm, task::Task{<:BaseTaskMap}, CM::Chart{I,M},
        CN::Chart{J,N}, robot_coord_rep=ChartRep()) where {M,N,I,J}
    m, n = dim(M), dim(N)
    (xm, vm) = robot_coord_rep == EmbRep() ? emb_to_chart_differential(xm, vm, CM) : (xm, vm)
    CN = isglobal(CN) ? CN : choose_chart_chart(xm, task, CM, CN)
    JftWJf, JftWA = single_task_components(xm, vm, task, CM, CN)
    σxddot = SMatrix{m,m,eltype(xm)}(pinv(Matrix(JftWJf)))*JftWA
    if robot_coord_rep == EmbRep()
        xme, vme, σxddot = chart_to_emb_differential(xm, vm, σxddot, CM)
    end
    σxddot, CN
end

function multiple_task_acceleration(xm, vm, tasks::TaskList, CM::Chart{I,M}, CNs::ChartList,
        robot_coord_rep=ChartRep(); log_task_chart=false) where {M,I}
    m = dim(M)
    JftWJf_sum = zeros(m,m)
    JftWA_sum = zeros(m)
    CNs_out = ChartList()

    (xm, vm) = robot_coord_rep == EmbRep() ? emb_to_chart_differential(xm, vm, CM) : (xm, vm)
    for i in 1:length(tasks)
        CN = isglobal(CNs[i]) ? CNs[i] : choose_chart_chart(xm, tasks[i], CM, CNs[i])
        JftWJf, JftWA = single_task_components(xm, vm, tasks[i], CM, CN)
        JftWJf_sum += JftWJf
        JftWA_sum += JftWA
        log_task_chart && push!(CNs_out, CN)
    end

    σxddot = SMatrix{m,m,eltype(xm)}(pinv(Matrix(JftWJf_sum)))*JftWA_sum
    if robot_coord_rep == EmbRep()
        xme, vme, σxddot = chart_to_emb_differential(xm, vm, σxddot, CM)
    end
    σxddot, CNs_out
end

function task_map_pullback!(JftWJfm_sum, Am_sum, Bm_sum, JftWginvℱm_sum, Γm_sum, xm, vm, 
        node::TreeNode{<:BaseTaskMap}, CM::Chart{I,M}, CN::Chart{J,N};
        log_tasks=false) where {M,N,I,J}
    m, n = dim(M), dim(N)
    JftWJfn_sum, An_sum, Bn_sum = zeros(n,n), zeros(n,n), zeros(n,n)
    JftWginvℱn_sum = zeros(n)
    Γn_sum = zeros(n,n,n)
        
    task_map = node.data
        
    Jf = task_jacobian_chart(xm, task_map, CM, CN)
    xn = task_map_chart(xm, task_map, CM, CN)
    vn = Jf*vm
    
    for child in children(node)
        if !isglobal(child.chart)
            CN_child = choose_chart_chart(xn, child.data, CN, child.chart)
        else
            CN_child = child.chart
        end
        task_map_pullback!(JftWJfn_sum, An_sum, Bn_sum, JftWginvℱn_sum, Γn_sum,
            xn, vn, child, CN, CN_child; log_tasks)
    end
    
    if Bn_sum != @SMatrix zeros(n,n)
        Jfdot = task_jacobian_chart_dot(xm, vm, task_map, CM, CN)
        Am = Jf'*(An_sum*Jf + Bn_sum*Jfdot)
        Bm = Jf'*Bn_sum*Jf
        Bm_sum .+= Bm
        log_tasks && push!(node.traj_log.B, Bm)
    else
        Am = Jf'*An_sum*Jf
    end
    Am_sum .+= Am

    if JftWJfn_sum != @SMatrix zeros(n,n)
        JftWJfm = Jf'*JftWJfn_sum*Jf
        JftWJfm_sum .+= JftWJfm
        log_tasks && push!(node.traj_log.JftWJf, JftWJfm)
        if JftWginvℱn_sum != @SVector zeros(n)
            JftWginvℱm = Jf'*JftWginvℱn_sum
            JftWginvℱm_sum .+= JftWginvℱm
            log_tasks && push!(node.traj_log.JftWginvℱ, JftWginvℱm)
        end
    end

    if Γn_sum != @SArray zeros(n,n,n)
        @tullio Γmnn[i,h,k] := Jf[l,i]*Γn_sum[l,h,k]
        @tullio Γmmn[l,s,k] := Jf[h,s]*Γmnn[l,h,k]
        @tullio Γm[l,h,q] := Jf[k,q]*Γmmn[l,h,k]
        Γm_sum .+= Γm
        log_tasks && push!(node.traj_log.Γ, Γm)
    end

    if log_tasks
        if node.coord_rep == ChartRep()
            push!(node.traj_log.x, xn)
            push!(node.traj_log.v, vn)
        else
            xne, vne = chart_to_emb_differential(xn, vn, CN)
            push!(node.traj_log.x, xne)
            push!(node.traj_log.v, vne)
        end
        push!(node.traj_log.chart, CN)
        push!(node.traj_log.A, Am)
    end

    nothing
end

function task_components!(JftWJfm_sum, Am_sum, Bm_sum, JftWginvℱm_sum, Γm_sum, xm, vm,
        node::TreeNode{<:Task{<:BaseTaskMap}}, CM::Chart{I,M}, CN::Chart{J,N};
        log_tasks=false) where {M,N,I,J}
    m, n = dim(M), dim(N)
    task = node.data
    
    Jf = task_jacobian_chart(xm, task, CM, CN)
    xn = task_map_chart(xm, task, CM, CN)
    vn = Jf*vm

    if log_tasks
        if node.coord_rep == ChartRep()
            push!(node.traj_log.x, xn)
            push!(node.traj_log.v, vn)
        else
            xne, vne = chart_to_emb_differential(xn, vn, CN)
            push!(node.traj_log.x, xne)
            push!(node.traj_log.v, vne)
        end
        push!(node.traj_log.chart, CN)
    end

    W = weight_metric_chart(xn, vn, task, CN)

    if W != @SMatrix zeros(n,n)
        JftW = Jf'*W
        Jfdot = task_jacobian_chart_dot(xm, vm, task, CM, CN)
        Γn = christoffel_symbols(xn, task, CN)
        if any(isinf.(Γn)) || any(isnan.(Γn))
            Γn = eltype(xm).(christoffel_symbols(BigFloat.(xn), task, CN))
        end
        g = metric_chart(xn, task, CN)
        
        ℱ_pot = potential_force_chart(xn, task, CN)
        ℱ_dis = dissipative_forces_chart(xn, vn, task, CN)
        ℱ = ℱ_pot + ℱ_dis
        JftWJf = Jf'*W*Jf
        JftWJfm_sum .+= JftWJf
        Bm = JftWJf
        Bm_sum .+= Bm

        if ℱ != @SVector zeros(n)
            JftWginvℱm = Jf'*W*inv(g)*ℱ
            JftWginvℱm_sum .+= JftWginvℱm
            log_tasks && push!(node.traj_log.JftWginvℱ, JftWginvℱm)
        end
        if Γn != @SArray zeros(n,n,n)
            @tullio Γmnn[i,h,k] := Jf[l,i]*Γn[l,h,k]
            @tullio Γmmn[l,s,k] := Jf[h,s]*Γmnn[l,h,k]
            @tullio Γm[l,h,q] := JftW[q,k]*Γmmn[l,h,k]
            Γm_sum .+= Γm
            log_tasks && push!(node.traj_log.Γ, Γm)
        end

        Am = Jf'*W*Jfdot
        Am_sum .+= Am
        
        if log_tasks
            push!(node.traj_log.g, g)
            push!(node.traj_log.ginv, inv(g))
            push!(node.traj_log.JftWJf, JftWJf)
            push!(node.traj_log.A, Am)
            push!(node.traj_log.B, Bm)
        end

    elseif log_tasks
        push!(node.traj_log.g, 0)
        push!(node.traj_log.ginv, 0)
        push!(node.traj_log.Γ, 0)
        push!(node.traj_log.JftWJf, 0)
        push!(node.traj_log.JftWginvℱ, 0)
        push!(node.traj_log.A, 0)
        push!(node.traj_log.B, 0)
    end

    nothing
end

function task_map_pullback!(JftWJfm_sum, Am_sum, Bm_sum, JftWginvℱm_sum, Γm_sum, xm, vm,
        node::TreeNode{<:TaskX}, CM, CN; log_tasks=false)
    task_components!(JftWJfm_sum, Am_sum, Bm_sum, JftWginvℱm_sum, Γm_sum, xm, vm,
        node, CM, CN; log_tasks)
end

function task_acceleration(xm, vm, node::TreeNode{M}, CM::Chart{I,M}, robot_coord_rep=ChartRep();
        log_tasks=false) where {M<:Manifold,I}
    m = dim(M)
    SVm, SMmm = SVector{m,eltype(xm)}, SMatrix{m,m,eltype(xm)}
    JftWJf_sum, A_sum, B_sum = zeros(m,m), zeros(m,m), zeros(m,m)
    JftWginvℱ_sum = zeros(m)
    Γ_sum = zeros(m,m,m)

    (xm, vm) = robot_coord_rep == EmbRep() ? emb_to_chart_differential(xm, vm, CM) : (xm, vm)
    for child in children(node)
        if !isglobal(child.chart)
            CN_child = choose_chart_chart(xm, child.data, CM, child.chart)
        else
            CN_child = child.chart
        end
        task_map_pullback!(JftWJf_sum, A_sum, B_sum, JftWginvℱ_sum, Γ_sum,
            xm, vm, child, CM, CN_child; log_tasks)
    end
    
    m_inds = static(1):static(m)
    vmΓm = SMmm([sum(Tuple(vm[l]*Γ_sum[l,h,k] for l=m_inds)) for k=m_inds, h=m_inds])
    vmΓmvm = vmΓm*vm

    JftWginvℱ_sum_st = SVm(JftWginvℱ_sum)
    A_sum_st = SMmm(A_sum)
    
    if log_tasks
        push!(node.traj_log.x, xm)
        push!(node.traj_log.v, vm)
        push!(node.traj_log.chart, CM)
        push!(node.traj_log.JftWJf, JftWJf_sum)
        push!(node.traj_log.JftWginvℱ, JftWginvℱ_sum_st)
        push!(node.traj_log.A, A_sum_st)
        push!(node.traj_log.vmΓmvm, vmΓmvm)
    end

    σxddot = SMmm(pinv(JftWJf_sum))*(JftWginvℱ_sum_st - A_sum_st*vm - vmΓmvm)
    if robot_coord_rep == EmbRep()
        xme, vme, σxddot = chart_to_emb_differential(xm, vm, σxddot, CM)
    end
    σxddot
end

# Version where mechanism state needs to be updated
function task_acceleration(xm, vm, node::TreeNode{M}, CM::Chart{I,M}, state::MechanismState, 
    robot_coord_rep=ChartRep(); log_tasks=false) where {M<:Manifold,I}
    set_configuration!(state, xm)
    task_acceleration(xm, vm, node, CM, robot_coord_rep; log_tasks)
end

# Version where cache is also present
function task_acceleration(xm, vm, node::TreeNode{M}, CM::Chart{I,M}, state::MechanismState, 
        cache::ControllerCache, robot_coord_rep=ChartRep(); log_tasks=false) where {M<:Manifold,I}
    set_configuration!(state, xm)
    setdirty!(cache)
    task_acceleration(xm, vm, node, CM, robot_coord_rep; log_tasks)
end