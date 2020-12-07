using AbstractTrees

struct TrajectoryLog
    x
    v
    chart

    g
    ginv
    Γ

    JftWJf
    JftWginvℱ
    A
    B
    vmΓmvm
    
    TrajectoryLog() = new([],[],[],[],[],[],[],[],[],[],[])
end

function clear!(traj_log::TrajectoryLog)
    empty!(traj_log.x)
    empty!(traj_log.v)
    empty!(traj_log.chart)
    empty!(traj_log.g)
    empty!(traj_log.ginv)
    empty!(traj_log.Γ)
    empty!(traj_log.JftWJf)
    empty!(traj_log.JftWginvℱ)
    empty!(traj_log.A)
    empty!(traj_log.B)
    empty!(traj_log.vmΓmvm)
end

mutable struct TreeNode{T}
    data::T # Manifold for root, Task for leaves, TaskMap otherwise
    children::Vector{TreeNode}
    chart::Chart
    coord_rep::CoordinateRep
    traj_log::TrajectoryLog
    parent::TreeNode

    function TreeNode(data::M; CM::Chart{I,M}=default_chart(M),
            robot_coord_rep=ChartRep()) where {M<:Manifold,I}
        new{M}(data, [], CM, robot_coord_rep, TrajectoryLog())
    end
    function TreeNode(data::T; CN=codomain_manifold(data), 
            task_coord_rep=ChartRep()) where T<:Union{TaskMapX,TaskX}
        new{T}(data, [], CN, task_coord_rep, TrajectoryLog())
    end
end

AbstractTrees.children(node::TreeNode) = node.children
AbstractTrees.printnode(io::IO, node::TreeNode{T}) where T = print(io, T)
function AbstractTrees.printnode(io::IO, node::TreeNode{T}, inds) where T
    if !isempty(inds)
        print(io, '(')
        show(IOContext(io, :compact => true), inds[end])
        print(io, ") ")
    end
    AbstractTrees.printnode(io, T)
end

function add_child!(node::TreeNode{M}, taskx::Union{TaskMapX,TaskX};
        CN=default_chart(codomain_manifold(taskx)), task_coord_rep=ChartRep()) where M <: Manifold
    if domain_manifold(taskx) == M
        new_node = TreeNode(taskx; CN, task_coord_rep)
        push!(node.children, new_node)
        return new_node
    else
        throw(ArgumentError("Domain manifold mismatch"))
    end
    nothing
end

function add_child!(node::TreeNode{<:TaskMapX}, taskx::Union{TaskMapX,TaskX};
        CN=default_chart(codomain_manifold(taskx)), task_coord_rep=ChartRep()) where M <: Manifold
    if codomain_manifold(node.data) == domain_manifold(taskx)
        new_node = TreeNode(taskx; CN, task_coord_rep)
        push!(node.children, new_node)
        return new_node
    else
        throw(ArgumentError("Manifold mismatch"))
    end
    nothing
end

add_child!(node::TreeNode{<:TaskX}, args...) = throw(ArgumentError("All tasks must be leaf nodes"))

function check_tree(root::TreeNode)
    for node in Leaves(root)
        if typeof(node.data) <: Union{Manifold, TaskMapX}
            error("Check failed: All leaf nodes must be tasks")
        end
    end
    println("Check passed")
end

function set_charts_emb!(node::TreeNode{M}, xme, CM) where M <: Manifold
    !isglobal(node.chart) ? node.chart = choose_chart_emb(xme, node.chart) : nothing
    for child in children(node)
        set_charts_emb!(child, xme, node.chart)
    end
end
function set_charts_chart!(node::TreeNode{M}, xm, CM) where M <: Manifold
    !isglobal(CM) ? node.chart = choose_chart_chart(xm, CM) : nothing
    for child in children(node)
        set_charts_chart!(child, xm, node.chart)
    end
end
function set_charts_emb!(node::TreeNode, xme, CM)
    taskx = node.data
    !isglobal(node.chart) ? node.chart = choose_chart_emb(xme, taskx, CM, node.chart) : nothing
    for child in children(node)
        set_charts_emb!(child, task_map_emb(xme, taskx, CM, node.chart), node.chart)
    end
end
function set_charts_chart!(node::TreeNode, xm, CM)
    taskx = node.data
    !isglobal(node.chart) ? node.chart = choose_chart_chart(xm, taskx, CM, node.chart) : nothing
    for child in children(node)
        set_charts_chart!(child, task_map_chart(xm, taskx, CM, node.chart), node.chart)
    end
end

function clear_logs!(node::TreeNode)
    clear!(node.traj_log)
    for child in node.children
        clear_logs!(child)
    end
    nothing
end