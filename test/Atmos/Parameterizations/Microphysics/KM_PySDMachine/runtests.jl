"""

julia --project="test/Atmos/Parameterizations/Microphysics/KM_PySDMachine" test/Atmos/Parameterizations/Microphysics/KM_PySDMachine/runtests.jl will run tests from the main ClimateMachine.jl directory.

"""

using Test, Pkg


begin
    root_folder_index = findlast("ClimateMachine.jl", pwd())
    tmp_path = pwd()[root_folder_index[1]:end]
    n = length(splitpath(tmp_path)) -1

    if n==0
        path = "."
    else
        path = repeat("../", n)
    end
    
    Pkg.add(url=path)
end

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
