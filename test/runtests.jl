using PBDS
using Test
using NBInclude

@testset "Examples" begin
    dir = joinpath(@__DIR__, "..", "examples", "PBDS")

    test_notebook = "R7Arm_DynamicMugGrasping"
    @testset "$test_notebook" begin
        file = string(test_notebook, ".ipynb")
        @nbinclude(joinpath(dir, file))
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
            @nbinclude(joinpath(dir, file))
            @test all([norm(traj.xm[end] - xm_goal) for traj in trajs] .< 5e-2)
            println("Finished example ", test_notebook)
        end
    end
end