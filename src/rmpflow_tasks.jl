function metric_curvature(xn, vn, task::TaskGDS{<:TaskMapT{M,N}}, C::Chart{I,N}) where {M,N,I}
    n = dim(N)
    ∂g_∂v = ForwardDiff.jacobian(vn -> reshape(metric_chart(xn, vn, task, C), Size(n^2)), vn)
    Ξ = SMatrix{n,n,eltype(xn)}(0.5*sum(∂g_∂v[(i-1)*n+1:i*n,:].*vn[i] for i in 1:n)')
end

function force_curvature(xn, vn, task::TaskGDS{<:TaskMapT{M,N}}, C::Chart{I,N}) where {M,N,I}
    n = dim(N)
    ∂g_∂x = SMatrix{n^2,n,eltype(xn)}(ForwardDiff.jacobian(xn -> reshape(metric_chart(xn, vn, task, C), n^2), xn))
    ξ1 = reshape(∂g_∂x*vn, Size(n,n))*vn
    ξ2 = ForwardDiff.gradient(xn -> (vn'*metric_chart(xn, vn, task, C)*vn)[1], xn)
    ξ = ξ1 - 0.5*ξ2
end

function single_task_components(xm, vm, task::TaskGDS{<:TaskMapT{M,N}}, CM::Chart{I,M},
        CN::Chart{J,N}) where {M,N,I,J}
    n = dim(N)
    xn = base_task_map_chart(xm, task, CM, CN)
    Jf = base_task_jacobian_chart(xm, task, CM, CN)
    vn = Jf*vm
    Jfdot = base_task_jacobian_chart_dot(xm, vm, task, CM, CN)
    
    G = metric_chart(xn, vn, task, CN)
    Ξ = metric_curvature(xn, vn, task, CN)
    ξ = force_curvature(xn, vn, task, CN)
    ℱ_pot = potential_force_chart(xn, task, CN)
    B = dissipative_term_chart(xn, vn, task, CN)
    A = Ξ + G

    xnddot_des = SVector{n,eltype(xm)}(pinv(A)*(ℱ_pot - B*vn - ξ))

    JftAJf = Jf'*A*Jf
    JftA_a = Jf'*A*(xnddot_des - Jfdot*vm)
    
    JftAJf, JftA_a
end

function single_task_acceleration(xm, vm, task::TaskGDS{<:TaskMapT{M,N}}, CM::Chart{I,M},
        CN::Chart{J,N}, robot_coord_rep=ChartRep()) where {M,N,I,J}
    m, n = dim(M), dim(N)
    (xm, vm) = robot_coord_rep == EmbRep() ? emb_to_chart_differential(xm, vm, CM) : (xm, vm)
    CN = isglobal(CN) ? CN : choose_chart_chart(xm, task, CM, CN)
    JftAJf, JftA_a = single_task_components(xm, vm, task, CM, CN)
    σxddot = SMatrix{m,m,eltype(xm)}(pinv(Matrix(JftAJf)))*JftA_a
    if robot_coord_rep == EmbRep()
        xme, vme, σxddot = chart_to_emb_differential(xm, vm, σxddot, CM)
    end
    σxddot, CN
end

function multiple_task_acceleration(xm, vm, tasks::TaskGDSList, CM::Chart{I,M}, CNs::ChartList,
        robot_coord_rep=ChartRep(); log_task_chart=false) where {M,I}
    m = dim(M)
    JftAJf_sum = zeros(m,m)
    JftA_a_sum = zeros(m)
    CNs_out = ChartList()

    (xm, vm) = robot_coord_rep == EmbRep() ? emb_to_chart_differential(xm, vm, CM) : (xm, vm)
    for i in 1:length(tasks)
        CN = isglobal(CNs[i]) ? CNs[i] : choose_chart_chart(xm, tasks[i], CM, CNs[i])
        JftAJf, JftA_a = single_task_components(xm, vm, tasks[i], CM, CN)
        JftAJf_sum += JftAJf
        JftA_a_sum += JftA_a
        log_task_chart && push!(CNs_out, CN)
    end
    σxddot = SMatrix{m,m,eltype(xm)}(pinv(Matrix(JftAJf_sum)))*JftA_a_sum
    if robot_coord_rep == EmbRep()
        xme, vme, σxddot = chart_to_emb_differential(xm, vm, σxddot, CM)
    end
    σxddot, CNs_out
end