# Test the branching of the SoilWaterParameterizations functions,
# test the instantiation of them works as expected, and that they
# return what we expect.
using MPI
using OrderedCollections
using StaticArrays
using Test

using CLIMAParameters
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

using ClimateMachine
using ClimateMachine.Land.SoilWaterParameterizations
using ClimateMachine.VariableTemplates

@testset "Land water parameterizations" begin
    FT = Float64
    test_array = [0.5, 1.0]
    vg_model = vanGenuchten(FT;)
    mm = MoistureDependent{FT}()
    bc_model = BrooksCorey(FT;)
    hk_model = Haverkamp(FT;)
    struct SoilWaterTestModel end
    function state(m::SoilWaterTestModel, T)
        @vars begin
            z::T
        end
    end
    model = SoilWaterTestModel()
    
    st = state(model, FT)
    v1 = Vars{st}(zeros(MVector{varsize(st), FT}))# want this to be an array of pts
    v2 = Vars{st}(ones(MVector{varsize(st), FT}))# want this to be an array of pts
    
    #Use an array to confirm that extra arguments are unused.
    @test viscosity_factor.(Ref(ConstantViscosity{FT}()), test_array) ≈
          [1.0, 1.0]
    @test impedance_factor.(Ref(NoImpedance{FT}()), test_array, test_array) ≈
          [1.0, 1.0]
    @test moisture_factor.(
        Ref(MoistureIndependent{FT}()),
        Ref(vg_model),
        test_array,[v1,v2]
    ) ≈ [1.0, 1.0]

    viscosity_model = TemperatureDependentViscosity{FT}(; T_ref = FT(1.0))
    @test viscosity_factor(viscosity_model, FT(1.0)) == 1
    impedance_model = IceImpedance{FT}(; Ω = 2.0)
    @test impedance_factor(impedance_model, 0.2, 0.2) == FT(0.1)

    @test moisture_factor(mm, vg_model, FT(1),v1) == 1
    @test moisture_factor(mm, vg_model, FT(0),v1) == 0


    @test moisture_factor(mm, bc_model, FT(1),v1) == 1
    @test moisture_factor(mm, bc_model, FT(0),v1) == 0


    @test moisture_factor(mm, hk_model, FT(1),v1) == 1
    @test moisture_factor(mm, hk_model, FT(0),v1) == 0


    @test hydraulic_conductivity(
        impedance_model,
        viscosity_model,
        MoistureDependent{FT}(),
        vanGenuchten(FT;),
        0.5,
        0.5,
        1.0,
        1.0,
        v1
    ) == FT(0.1)

    @test_throws Exception effective_saturation(0.5, -1.0)
    @test effective_saturation(0.5, 0.25) == 0.5

    test_array = [0.5, 1.0]
    n = FT(1.43)
    m = 1.0 - 1.0 / n
    α = FT(2.6)

    @test pressure_head.(
        Ref(vg_model),
        Ref(1.0),
        Ref(0.001),
        test_array,
        Ref(0.0),
        Ref(v)
    ) ≈ .-((-1 .+ test_array .^ (-1 / m)) .* α^(-n)) .^ (1 / n)
    #test branching in pressure head
    @test pressure_head(vg_model, 1.0, 0.001, 1.5, 0.0, v) == 500

    @test pressure_head.(
        Ref(hk_model),
        Ref(1.0),
        Ref(0.001),
        test_array,
        Ref(0.0),
    ) ≈ .-((-1 .+ test_array .^ (-1 / m)) .* α^(-n)) .^ (1 / n)


    m = FT(0.5)
    ψb = FT(0.1656)
    @test pressure_head(bc_model, 1.0, 0.001, 0.5, 0.0) ≈ -ψb * 0.5^(-1 / m)
    @test volumetric_liquid_fraction.([0.5, 1.5], Ref(0.75)) ≈ [0.5, 0.75]

    @test inverse_matric_potential(
        bc_model,
        matric_potential(bc_model, FT(0.5)),
    ) == FT(0.5)
    @test inverse_matric_potential(
        vg_model,
        matric_potential(vg_model, FT(0.5)),
    ) == FT(0.5)
    @test inverse_matric_potential(
        hk_model,
        matric_potential(hk_model, FT(0.5)),
    ) == FT(0.5)

    @test_throws Exception inverse_matric_potential(bc_model, 1)
    @test_throws Exception inverse_matric_potential(hk_model, 1)
    @test_throws Exception inverse_matric_potential(vg_model, 1)


end
