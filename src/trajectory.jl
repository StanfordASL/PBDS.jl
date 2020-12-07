struct RMPTrajectory{M,S,m}
    xm::Vector{SVector{m,S}}
    vm::Vector{SVector{m,S}}
    am::Vector{SVector{m,S}}
    robot_coord_rep::CoordinateRep
    CM::Vector{Chart}
    
    xns::Vector{Vector}
    vns::Vector{Vector}
    task_coord_reps::Vector{CoordinateRep}
    CNs::Vector{Vector{Chart}}

    Time::S
    dt::S
end

function RMPTrajectory(xm_in, vm_in, Time, dt, task::TaskX, CM::Chart{I,M}, CN,
        robot_coord_rep::ChartRep, task_coord_rep=ChartRep(), log_tasks=false) where {M,I}
    
    am, CN = single_task_acceleration(xm_in, vm_in, task, CM, CN, robot_coord_rep)
    xm0, vm0, am0 = [xm_in], [vm_in], [am]

    xns, vns, CNs = [], [], [CN]
    if log_tasks
        if task_coord_rep == ChartRep()
            xn, vn = task_differential_map_chart(xm_in, vm_in, task, CM, CN)
            push!(xns, xn)
            push!(vns, vn)
        else
            xne, vne = task_differential_map_chart_emb(xm_in, vm_in, task, CM, CN)
            push!(xns, xne)
            push!(vns, vne)
        end
    end

    RMPTrajectory{M,eltype(xm_in),dim(M)}(xm0, vm0, am0, robot_coord_rep, [CM],
        xns, vns, [task_coord_rep], [CNs], Time, dt)
end

function RMPTrajectory(xme_in, vme_in, Time, dt, task::TaskX, CM::Chart{I,M}, CN,
    robot_coord_rep::EmbRep, task_coord_rep=ChartRep(), log_tasks=false) where {M,I}

    CM = choose_chart_emb(xme_in, CM)
    ame, CN = single_task_acceleration(xme_in, vme_in, task, CM, CN, robot_coord_rep)
    xm0, vm0, am0 = [xme_in], [vme_in], [ame]

    xns, vns, CNs = [], [], [CN]
    if log_tasks
        if task_coord_rep == ChartRep()
            xn, vn = task_differential_map_emb_chart(xme_in, vme_in, task, CM, CN)
            push!(xns, xn)
            push!(vns, vn)
        else
            xne, vne = task_differential_map_emb(xme_in, vme_in, task, CM, CN)
            push!(xns, xne)
            push!(vns, vne)
        end
    end

    RMPTrajectory{M,eltype(xme_in),embdim(M)}(xm0, vm0, am0, robot_coord_rep, [CM],
        xns, vns, [task_coord_rep], [CNs], Time, dt)
end

function RMPTrajectory(xm_in, vm_in, Time, dt, tasks::TaskXList, CM::Chart{I,M}, CNs::ChartList,
        robot_coord_rep::ChartRep, task_coord_reps=fill(ChartRep(), length(tasks)),
        log_tasks=false) where {M,I}

    am, CNs = multiple_task_acceleration(xm_in, vm_in, tasks, CM, CNs,
        robot_coord_rep, log_task_chart=true)
    xm0, vm0, am0 = [xm_in], [vm_in], [am]

    xns, vns = [], []
    if log_tasks
        for i = 1:length(tasks)
            if task_coord_reps[i] == ChartRep()
                xn, vn = task_differential_map_chart(xm_in, vm_in, tasks[i], CM, CNs[i])
                push!(xns, [xn])
                push!(vns, [vn])
            else
                xne, vne = task_differential_map_chart_emb(xm_in, vm_in, tasks[i], CM, CNs[i])
                push!(xns, [xne])
                push!(vns, [vne])
            end
        end
    end

    RMPTrajectory{M,eltype(xm_in),dim(M)}(xm0, vm0, am0, robot_coord_rep, [CM],
        xns, vns, task_coord_reps, [CNs], Time, dt)
end

function RMPTrajectory(xme_in, vme_in, Time, dt, tasks::TaskXList, CM::Chart{I,M}, CNs::ChartList,
    robot_coord_rep::EmbRep, task_coord_reps=fill(ChartRep(), length(tasks)),
    log_tasks=false) where {M,I}

    CM = choose_chart_emb(xme_in, CM)
    ame, CNs = multiple_task_acceleration(xme_in, vme_in, tasks, CM, CNs,
        robot_coord_rep, log_task_chart=true)
    xm0, vm0, am0 = [xme_in], [vme_in], [ame]

    xns, vns = [], []
    if log_tasks
        for i = 1:length(tasks)
            if task_coord_reps[i] == ChartRep()
                xn, vn = task_differential_map_emb_chart(xme_in, vme_in, tasks[i], CM, CNs[i])
                push!(xns, [xn])
                push!(vns, [vn])
            else
                xne, vne = task_differential_map_emb(xme_in, vme_in, tasks[i], CM, CNs[i])
                push!(xns, [xne])
                push!(vns, [vne])
            end
        end
    end

    RMPTrajectory{M,eltype(xme_in),embdim(M)}(xm0, vm0, am0, robot_coord_rep, [CM],
        xns, vns, task_coord_reps, [CNs], Time, dt)
end

function rk4step(x, v, a, dt, f, args...)
    k1x = x
    k1v = v
    k1a = a

    k2x = k1x + dt*k1v/2
    k2v = k1v + dt*k1a/2
    k2a = f(k2x, k2v, args...)

    k3x = k1x + dt*k2v/2
    k3v = k1v + dt*k2a/2
    k3a = f(k3x, k3v, args...)

    k4x = k1x + dt*k3v
    k4v = k1v + dt*k3a
    k4a = f(k4x, k4v, args...)

    x_next = k1x + dt*(k1v + 2*(k2v + k3v) + k4v)/6
    v_next = k1v + dt*(k1a + 2*(k2a + k3a) + k4a)/6

    x_next, v_next
end

function rk4step_timedep(x, v, a, t, dt, f, args...)
    k1x = x
    k1v = v
    k1a = a

    k2x = k1x + dt*k1v/2
    k2v = k1v + dt*k1a/2
    k2a = f(k2x, k2v, t + dt/2, args...)

    k3x = k1x + dt*k2v/2
    k3v = k1v + dt*k2a/2
    k3a = f(k3x, k3v, t + dt/2, args...)

    k4x = k1x + dt*k3v
    k4v = k1v + dt*k3a
    k4a = f(k4x, k4v, t + dt, args...)

    x_next = k1x + dt*(k1v + 2*(k2v + k3v) + k4v)/6
    v_next = k1v + dt*(k1a + 2*(k2a + k3a) + k4a)/6

    x_next, v_next
end

propagate_task(xm, vm, task, CM, CNs, Time, dt) =
    propagate_task(xm, vm, task, CM, CNs, Time, dt, ChartRep())

function propagate_task(xm, vm, task::TaskX, CM, CN, Time, dt,
        robot_coord_rep::ChartRep; task_coord_rep=ChartRep(), log_tasks=false)
    Nsteps = Int(cld(Time, dt))
    traj = RMPTrajectory(xm, vm, Time, dt, task, CM, CN, robot_coord_rep, task_coord_rep, log_tasks)
    am, CN = traj.am[1], traj.CNs[1][1]

    for i = 1:Nsteps-1
        xm, vm = rk4step(xm, vm, am, dt,
            (x, v, args...) -> single_task_acceleration(x, v, args...)[1],
            task, CM, CN, robot_coord_rep)

        CM_next = choose_chart_chart(xm, CM)
        if CM_next != CM
            xm, vm = chart_transition_differential(xm, vm, CM, CM_next)
            CM = CM_next
        end
        am, CN = single_task_acceleration(xm, vm, task, CM, CN, robot_coord_rep)
        push!(traj.xm, xm)
        push!(traj.vm, vm)
        push!(traj.am, am)
        push!(traj.CM, CM)
        
        if log_tasks
            if task_coord_rep == ChartRep()
                xn, vn = task_differential_map_chart(xm, vm, task, CM, CN)
                push!(traj.xns, xn)
                push!(traj.vns, vn)
            else
                xne, vne = task_differential_map_chart_emb(xm, vm, task, CM, CN)
                push!(traj.xns, xne)
                push!(traj.vns, vne)
            end
        end
        push!(traj.CNs, [CN])
    end

    return traj
end

function propagate_task(xme, vme, task::TaskX, CM, CN, Time, dt,
        robot_coord_rep::EmbRep; task_coord_rep=ChartRep(), log_tasks=false)
    Nsteps = Int(cld(Time, dt))
    traj = RMPTrajectory(xme, vme, Time, dt, task, CM, CN, robot_coord_rep, task_coord_rep, log_tasks)
    ame, CM, CN = traj.am[1], traj.CM[1], traj.CNs[1][1]
    xm, vm, am = emb_to_chart_differential(xme, vme, ame, CM)

    for i = 1:Nsteps-1
        xm, vm = rk4step(xm, vm, am, dt,
            (x, v, args...) -> single_task_acceleration(x, v, args...)[1],
            task, CM, CN, ChartRep())
        CM_next = choose_chart_chart(xm, CM)
        if CM_next != CM
            xm, vm = chart_transition_differential(xm, vm, CM, CM_next)
            CM = CM_next
        end
        am, CN = single_task_acceleration(xme, vme, task, CM, CN, ChartRep())
        xme, vme, ame = chart_to_emb_differential(xm, vm, am, CM)
        push!(traj.xm, xme)
        push!(traj.vm, vme)
        push!(traj.am, ame)
        push!(traj.CM, CM)
        
        if log_tasks
            if task_coord_rep == ChartRep()
                xn, vn = task_differential_map_emb_chart(xme, vme, task, CM, CN)
                push!(traj.xns, xn)
                push!(traj.vns, vn)
            else
                xne, vne = task_differential_map_emb(xme, vme, task, CM, CN)
                push!(traj.xns, xne)
                push!(traj.vns, vne)
            end
        end
        push!(traj.CNs, [CN])
    end

    return traj
end

propagate_tasks(xm, vm, tasks, CM, CNs, Time, dt) =
    propagate_tasks(xm, vm, tasks, CM, CNs, Time, dt, ChartRep())

function propagate_tasks(xm, vm, tasks::TaskXList, CM, CNs::ChartList, Time, dt,
        robot_coord_rep::ChartRep; task_coord_reps=fill(ChartRep(), length(tasks)), log_tasks=false)
    Nsteps = Int(cld(Time, dt))
    Ntasks = length(tasks)
    traj = RMPTrajectory(xm, vm, Time, dt, tasks, CM, CNs,
        robot_coord_rep, task_coord_reps, log_tasks)
    xm, vm, am, CNs = traj.xm[1], traj.vm[1], traj.am[1], traj.CNs[1]
    
    for i = 1:Nsteps-1
        xm, vm = rk4step(xm, vm, am, dt,
            (x, v, args...) -> multiple_task_acceleration(x, v, args...)[1],
            tasks, CM, CNs, robot_coord_rep)

        CM_next = choose_chart_chart(xm, CM)
        if CM_next != CM
            xm, vm = chart_transition_differential(xm, vm, CM, CM_next)
            CM = CM_next
        end
        am, CNs = multiple_task_acceleration(xm, vm, tasks, CM, CNs,
            robot_coord_rep, log_task_chart=true)
        push!(traj.xm, xm)
        push!(traj.vm, vm)
        push!(traj.am, am)
        push!(traj.CM, CM)
        
        if log_tasks
            for j in 1:Ntasks
                if task_coord_reps[j] == ChartRep()
                    xn, vn = task_differential_map_chart(xm, vm, tasks[j], CM, CNs[j])
                    push!(traj.xns[j], xn)
                    push!(traj.vns[j], vn)
                else
                    xne, vne = task_differential_map_chart_emb(xm, vm, tasks[j], CM, CNs[j])
                    push!(traj.xns[j], xne)
                    push!(traj.vns[j], vne)
                end
            end
            push!(traj.CNs, CNs)
        end
    end

    return traj
end

function propagate_tasks(xme, vme, tasks::TaskXList, CM, CNs::ChartList, Time, dt,
        robot_coord_rep::EmbRep; task_coord_reps=fill(ChartRep(), length(tasks)), log_tasks=false)
    Nsteps = Int(cld(Time, dt))
    Ntasks = length(tasks)
    traj = RMPTrajectory(xme, vme, Time, dt, tasks, CM, CNs,
        robot_coord_rep, task_coord_reps, log_tasks)
    ame, CM, CN = traj.am[1], traj.CM[1], traj.CNs[1][1]
    xm, vm, am = emb_to_chart_differential(xme, vme, ame, CM)

    for i = 1:Nsteps-1
        xm, vm = rk4step(xm, vm, am, dt,
            (x, v, args...) -> multiple_task_acceleration(x, v, args...)[1],
            tasks, CM, CNs, ChartRep())
        CM_next = choose_chart_chart(xm, CM)
        if CM_next != CM
            xm, vm = chart_transition_differential(xm, vm, CM, CM_next)
            CM = CM_next
        end
        am, CN = multiple_task_acceleration(xm, vm, tasks, CM, CNs, ChartRep())
        xme, vme, ame = chart_to_emb_differential(xm, vm, am, CM)
        push!(traj.xm, xme)
        push!(traj.vm, vme)
        push!(traj.am, ame)
        push!(traj.CM, CM)
        
        if log_tasks
            for j in 1:Ntasks
                if task_coord_reps[j] == ChartRep()
                    xn, vn = task_differential_map_emb_chart(xme, vme, tasks[j], CM, CNs[j])
                    push!(traj.xns[j], xn)
                    push!(traj.vns[j], vn)
                else
                    xne, vne = task_differential_map_emb(xme, vme, tasks[j], CM, CNs[j])
                    push!(traj.xns[j], xne)
                    push!(traj.vns[j], vne)
                end
            end
            push!(traj.CNs, CNs)
        end
    end

    return traj
end

struct RMPTreeTrajectory{M,S,m}
    xm::Vector{SVector{m,S}}
    vm::Vector{SVector{m,S}}
    am::Vector{SVector{m,S}}
    robot_coord_rep::CoordinateRep
    CM::Vector{Chart}
    root::TreeNode{M}

    Time::S
    dt::S
end

function RMPTreeTrajectory(xm_in, vm_in, Time, dt, root::TreeNode{M}, CM,
        robot_coord_rep::ChartRep, args...; log_tasks=false, time_dep=false) where M <: Manifold
    
    log_tasks && clear_logs!(root)
    if !time_dep
        am = task_acceleration(xm_in, vm_in, root, CM, args..., robot_coord_rep; log_tasks)
    else
        am = task_acceleration(xm_in, vm_in, 0., root, CM, args..., robot_coord_rep; log_tasks)
    end
    xm0, vm0, am0 = [xm_in], [vm_in], [am]

    RMPTreeTrajectory{M,eltype(xm_in),dim(M)}(xm0, vm0, am0, robot_coord_rep, [CM], root,
        Time, dt)
end

function RMPTreeTrajectory(xme_in, vme_in, Time, dt, root::TreeNode{M}, CM,
        robot_coord_rep::EmbRep, args...; log_tasks=false, time_dep=false) where M <: Manifold

    log_tasks && clear_logs!(root)
    CM = choose_chart_emb(xme_in, CM)
    if !time_dep
        ame = task_acceleration(xme_in, vme_in, root, CM, args..., robot_coord_rep; log_tasks)
    else
        ame = task_acceleration(xme_in, vme_in, 0., root, CM, args..., robot_coord_rep; log_tasks)
    end
    xm0, vm0, am0 = [xme_in], [vme_in], [ame]

    RMPTreeTrajectory{M,eltype(xme_in),embdim(M)}(xm0, vm0, am0, robot_coord_rep, [CM], root,
        Time, dt)
end

function propagate_tasks(xm, vm, root::TreeNode{M}, CM, Time, dt,
        robot_coord_rep::ChartRep, args...; log_tasks=false, time_dep=false) where M <: Manifold
    Nsteps = Int(cld(Time, dt))
    set_charts_chart!(root, xm, CM)
    traj = RMPTreeTrajectory(xm, vm, Time, dt, root, CM, robot_coord_rep, args...;
        log_tasks, time_dep)
    am = traj.am[1]
    t = 0.
    
    for i = 1:Nsteps-1
        if !time_dep
            xm, vm = rk4step(xm, vm, am, dt,
                (x, v, args...) -> task_acceleration(x, v, args...),
                root, CM, args..., robot_coord_rep)
        else    
            xm, vm = rk4step_timedep(xm, vm, am, t, dt,
                (x, v, args...) -> task_acceleration(x, v, args...),
                root, CM, args..., robot_coord_rep)
        end

        CM_next = choose_chart_chart(xm, CM)
        if CM_next != CM
            xm, vm = chart_transition_differential(xm, vm, CM, CM_next)
            CM = CM_next
        end
        if !time_dep
            am = task_acceleration(xm, vm, root, CM, args..., robot_coord_rep; log_tasks)
        else
            am = task_acceleration(xm, vm, t, root, CM, args..., robot_coord_rep; log_tasks)
        end
        push!(traj.xm, xm)
        push!(traj.vm, vm)
        push!(traj.am, am)
        push!(traj.CM, CM)
        t += dt
    end

    return traj
end

function propagate_tasks(xme, vme, root::TreeNode{M}, CM, Time, dt,
        robot_coord_rep::EmbRep, args...; log_tasks=false, time_dep=false) where M <: Manifold
    Nsteps = Int(cld(Time, dt))
    set_charts_emb!(root, xme, CM)
    traj = RMPTreeTrajectory(xme, vme, Time, dt, root, CM, robot_coord_rep, args...;
        log_tasks, time_dep)
    ame, CM = traj.am[1], traj.CM[1]
    xm, vm, am = emb_to_chart_differential(xme, vme, ame, CM)
    t = 0.

    for i = 1:Nsteps-1
        if !time_dep
            xm, vm = rk4step(xm, vm, am, dt,
                (x, v, args...) -> task_acceleration(x, v, args...),
                root, CM, args..., ChartRep())
        else
            xm, vm = rk4step_timedep(xm, vm, am, t, dt,
                (x, v, args...) -> task_acceleration(x, v, args...),
                root, CM, args..., ChartRep())
        end
        CM_next = choose_chart_chart(xm, CM)
        if CM_next != CM
            xm, vm = chart_transition_differential(xm, vm, CM, CM_next)
            CM = CM_next
        end
        if !time_dep
            am = task_acceleration(xm, vm, root, CM, args..., ChartRep(); log_tasks)
        else
            am = task_acceleration(xm, vm, t, root, CM, args..., ChartRep(); log_tasks)
        end
        xme, vme, ame = chart_to_emb_differential(xm, vm, am, CM)
        push!(traj.xm, xme)
        push!(traj.vm, vme)
        push!(traj.am, ame)
        push!(traj.CM, CM)
        t += dt
    end

    return traj
end