"""

julia --project test/PySDMCall/runtests.jl will run the experiment from the main ClimateMachine.jl directory.

"""

using Test

@testset "PySDMCall tests" begin
    @testset "PyCall invocation" begin
        include(joinpath("PyCall_invocation.jl"))
    end
    @testset "PySDM presence" begin
        include(joinpath("PySDM_presence.jl"))
    end
    @testset "PySDMCallback invocation" begin
        include(joinpath("PySDMCallback_invocation.jl"))
    end
    @testset "PySDM cloud flatness" begin
        include(joinpath("KM_cloud_flatness.jl"))
    end
end
