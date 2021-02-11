using PBDS
using Test
using NBInclude

@testset "Examples" begin
    "no_plots" in ARGS && (global const no_plots = true)
    PBDS_dir = joinpath(@__DIR__, "..", "examples", "PBDS")

    test_notebook = "R7Arm_DynamicMugGrasping"
    @testset "$test_notebook" begin
        file = string(test_notebook, ".ipynb")
        @nbinclude(joinpath(PBDS_dir, file))
        @test root.children[end-1].children[3].traj_log.x[end][1] .< 5e-3
        robot_coord_rep = ChartRep()
        traj = propagate_tasks(xm, vm, root, CM, Time, dt, robot_coord_rep, state, cache, 
            mugparams; time_dep, log_tasks=true)
        @test root.children[end-1].children[3].traj_log.x[end][1] .< 5e-3
        println("Finished example ", test_notebook)
    end

    test_notebooks = ["R1_To_R1PointPositionAttractor",
                      "R2_To_R1PointDistanceAttractor_S1Damping",
                      "R2_To_R1PointDistanceAttractor_S1Damping_R1BoxAvoidance",
                      "R2_To_R1PointDistanceAttractor_S1Damping_R1SphereAvoidance",
                      "R2_To_R2PointPositionAttractor",
                      "R3_To_R1PointDistanceAttractor_S2Damping_R1BoxAvoidance",
                      "R3_To_R1PointDistanceAttractor_S2Damping_R1SphereAvoidance",
                      "S2_To_R1Attractor_S2Damping_R1ObstacleAvoidance"]

    for test_notebook in test_notebooks
        @testset "$test_notebook" begin
            file = string(test_notebook, ".ipynb")
            @nbinclude(joinpath(PBDS_dir, file))
            ε = 5e-2
            @test all([norm(traj.xm[end] - xm_goal) for traj in trajs] .< ε)

            if test_notebook == "R3_To_R1PointDistanceAttractor_S2Damping_R1SphereAvoidance"
                log_tasks = true
                robot_coord_rep = ChartRep()
                traj = propagate_tasks(xm, vm, tasks, CM, CNs, Time, dt, robot_coord_rep; log_tasks)
                @test norm(traj.xm[end] - xm_goal) < ε
                robot_coord_rep = EmbRep()
                traj = propagate_tasks(xm, vm, tasks, CM, CNs, Time, dt, robot_coord_rep; log_tasks)
                @test norm(traj.xm[end] - xm_goal) < ε
            end
            println("Finished example ", test_notebook)
        end
    end

    test_notebook = "S2_To_R1Attractor_S2Damping_ConsistencyTest"
    @testset "$test_notebook" begin
        file = string(test_notebook, ".ipynb")
        @nbinclude(joinpath(PBDS_dir, file))
        Δx = 0.
        Δx += sum(@. norm(traj_north.xm - traj_south.xm))
        Δx += sum(@. norm(traj_south.xm - traj_switching.xm))
        Δx += sum(@. norm(traj_switching.xm - traj_north.xm))
        @test Δx < 1e-4
        ε = 5e-3
        @test norm(traj_north.xm[end] - xm_goal) < ε

        RMP_dir = joinpath(@__DIR__, "..", "examples", "RMPflow")
        @nbinclude(joinpath(RMP_dir, file))
        Δx = 0.
        Δx += sum(@. norm(traj_north.xm - traj_south.xm))
        Δx += sum(@. norm(traj_south.xm - traj_switching.xm))
        Δx += sum(@. norm(traj_switching.xm - traj_north.xm))
        @test Δx > 1e3
        Δx = 0.
        Δx += norm(traj_north.xm[end] - xm_goal)
        Δx += norm(traj_south.xm[end] - xm_goal)
        Δx += norm(traj_switching.xm[end] - xm_goal)
        @test Δx < ε
        println("Finished example ", test_notebook)
    end

end