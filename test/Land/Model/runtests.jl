using Test, Pkg
@testset "Land" begin
    include("test_heat_parameterizations.jl")
    include("test_water_parameterizations.jl")
    include("prescribed_twice.jl")
    include("freeze_thaw_alone.jl")
    include("test_runoff_functions.jl")
    include("test_no_river.jl")
end
