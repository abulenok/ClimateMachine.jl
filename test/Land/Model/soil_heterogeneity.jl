using MPI
using OrderedCollections
using StaticArrays
using Statistics
using Dierckx
using Test
using Pkg.Artifacts
using DelimitedFiles

using CLIMAParameters
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

using ClimateMachine
using ClimateMachine.Land
using ClimateMachine.Land.SoilWaterParameterizations
using ClimateMachine.Mesh.Topologies
using ClimateMachine.Mesh.Grids
using ClimateMachine.DGMethods
using ClimateMachine.DGMethods.NumericalFluxes
using ClimateMachine.DGMethods: BalanceLaw, LocalGeometry
using ClimateMachine.MPIStateArrays
using ClimateMachine.GenericCallbacks
using ClimateMachine.ODESolvers
using ClimateMachine.VariableTemplates
using ClimateMachine.SingleStackUtils
using ClimateMachine.BalanceLaws:
    BalanceLaw, Prognostic, Auxiliary, Gradient, GradientFlux, vars_state
#=
using ClimateMachine.ArtifactWrappers

haverkamp_dataset = ArtifactWrapper(
    @__DIR__,
    isempty(get(ENV, "CI", "")),
    "richards",
    ArtifactFile[ArtifactFile(
        url = "https://caltech.box.com/shared/static/dfijf07io7h5dk1k87saaewgsg9apq8d.csv",
        filename = "bonan_haverkamp_data.csv",
    ),],
)
haverkamp_dataset_path = get_data_folder(haverkamp_dataset)

@testset "Richard's equation - Haverkamp test" begin
    ClimateMachine.init()
    FT = Float64

    function init_soil_water!(land, state, aux, localgeo, time)
        myfloat = eltype(aux)
        state.soil.water.ϑ_l = myfloat(land.soil.water.initialϑ_l(aux))
        state.soil.water.θ_i = myfloat(land.soil.water.initialθ_i(aux))
    end

    soil_heat_model = PrescribedTemperatureModel()

    soil_param_functions = SoilParamFunctions{FT}(
        porosity = 0.495,
    )
    # Keep in mind that what is passed is aux⁻.
    # Scalar fluxes are multiplied by ẑ (normal to the surface, -normal to the bottom,
    # where normal point outs of the domain.)
    surface_state = (aux, t) -> eltype(aux)(0.494)
    bottom_flux = (aux, t) -> aux.soil.water.K * eltype(aux)(-1)
    ϑ_l0 = (aux) -> eltype(aux)(0.24)
    sigmoid(x, offset, width) = typeof(x)(exp((x-offset)/width)/(1+exp((x-offset)/width)))
    bc = LandDomainBC(
        bottom_bc = LandComponentBC(soil_water = Neumann(bottom_flux)),
        surface_bc = LandComponentBC(soil_water = Dirichlet(surface_state)),
    )
    Ksat = (aux) -> eltype(aux)(0.0443 / (3600 * 100)*(9*sigmoid(aux.z,-1.0,0.02)+1)/10)
    S_s = FT(1e-3)
    soil_water_model = SoilWaterModel(
        FT;
        moisture_factor = MoistureDependent{FT}(),
        hydraulics = Haverkamp(FT;),
        Ksat = Ksat,
        S_s = S_s,
        initialϑ_l = ϑ_l0,
    )

    m_soil = SoilModel(soil_param_functions, soil_water_model, soil_heat_model)
    sources = ()
    m = LandModel(
        param_set,
        m_soil;
        boundary_conditions = bc,
        source = sources,
        init_state_prognostic = init_soil_water!,
    )


    N_poly = 5
    nelem_vert = 20

    # Specify the domain boundaries
    zmax = FT(0)
    zmin = FT(-2)

    driver_config = ClimateMachine.SingleStackConfiguration(
        "LandModel",
        N_poly,
        nelem_vert,
        zmax,
        param_set,
        m;
        zmin = zmin,
        numerical_flux_first_order = CentralNumericalFluxFirstOrder(),
    )

    t0 = FT(0)
    timeend = FT(60 * 60 * 24)

    dt = FT(6)

    solver_config = ClimateMachine.SolverConfiguration(
        t0,
        timeend,
        driver_config,
        ode_dt = dt,
    )
    mygrid = solver_config.dg.grid
    Q = solver_config.Q
    aux = solver_config.dg.state_auxiliary

    ClimateMachine.invoke!(solver_config)
    ϑ_l_ind = varsindex(vars_state(m, Prognostic(), FT), :soil, :water, :ϑ_l)
    ϑ_l = Array(Q[:, ϑ_l_ind, :][:])
    z_ind = varsindex(vars_state(m, Auxiliary(), FT), :z)
    z = Array(aux[:, z_ind, :][:])

    # Compare with Bonan simulation data at 1 day.
    data = joinpath(haverkamp_dataset_path, "bonan_haverkamp_data.csv")
    ds_bonan = readdlm(data, ',')
    bonan_moisture = reverse(ds_bonan[:, 1])
    bonan_z = reverse(ds_bonan[:, 2]) ./ 100.0


    # Create an interpolation from the Bonan data
    bonan_moisture_continuous = Spline1D(bonan_z, bonan_moisture)
    bonan_at_clima_z = bonan_moisture_continuous.(z)
    MSE = mean((bonan_at_clima_z .- ϑ_l) .^ 2.0)
    @test MSE < 1e-5
end
=#


@testset "hydrostatic test 1" begin
    ClimateMachine.init()
    FT = Float64

    function init_soil_water!(land, state, aux, localgeo, time)
        myfloat = eltype(aux)
        state.soil.water.ϑ_l = myfloat(land.soil.water.initialϑ_l(aux))
        state.soil.water.θ_i = myfloat(land.soil.water.initialθ_i(aux))
    end

    soil_heat_model = PrescribedTemperatureModel()

    soil_param_functions = SoilParamFunctions{FT}(
        porosity = 0.495,
    )
    # Keep in mind that what is passed is aux⁻.
    # Scalar fluxes are multiplied by ẑ (normal to the surface, -normal to the bottom,
    # where normal point outs of the domain.)
    bottom_flux = (aux, t) -> eltype(aux)(0.0)
    surface_flux = bottom_flux
    ϑ_l0 = (aux) -> eltype(aux)(0.494)
    bc = LandDomainBC(
        bottom_bc = LandComponentBC(soil_water = Neumann(bottom_flux)),
        surface_bc = LandComponentBC(soil_water = Neumann(surface_flux)),
    )
    Ksat = (aux) -> eltype(aux)(0.0443 / (3600 * 100))
    S_s = (aux) -> eltype(aux)((1e-3)*exp(-0.2*aux.z))
    vgn = FT(2)
    soil_water_model = SoilWaterModel(
        FT;
        moisture_factor = MoistureDependent{FT}(),
        hydraulics = vanGenuchten(FT;n = vgn),
        Ksat = Ksat,
        S_s = S_s,
        initialϑ_l = ϑ_l0,
    )

    m_soil = SoilModel(soil_param_functions, soil_water_model, soil_heat_model)
    sources = ()
    m = LandModel(
        param_set,
        m_soil;
        boundary_conditions = bc,
        source = sources,
        init_state_prognostic = init_soil_water!,
    )


    N_poly = 2
    nelem_vert = 20

    # Specify the domain boundaries
    zmax = FT(0)
    zmin = FT(-10)

    driver_config = ClimateMachine.SingleStackConfiguration(
        "LandModel",
        N_poly,
        nelem_vert,
        zmax,
        param_set,
        m;
        zmin = zmin,
        numerical_flux_first_order = CentralNumericalFluxFirstOrder(),
    )

    t0 = FT(0)
    timeend = FT(60 * 60 * 24*200)

    dt = FT(500)
    n_outputs = 3
    every_x_simulation_time = ceil(Int, timeend / n_outputs)
    solver_config = ClimateMachine.SolverConfiguration(
        t0,
        timeend,
        driver_config,
        ode_dt = dt,
    )
    aux = solver_config.dg.state_auxiliary
    state_types = (Prognostic(),Auxiliary())
    dons_arr =
        Dict[dict_of_nodal_states(solver_config, state_types; interp = true)]
    time_data = FT[0]
    callback = GenericCallbacks.EveryXSimulationTime(every_x_simulation_time) do
        dons = dict_of_nodal_states(solver_config, state_types; interp = true)
        push!(dons_arr, dons)
        push!(time_data, gettime(solver_config.solver))
        nothing
    end
    ClimateMachine.invoke!(solver_config; user_callbacks = (callback,))
    z = dons_arr[1]["z"]
    interface_z = -1.0395
    function hydrostatic_profile(z, zm, porosity, n, α)
        myf = eltype(z)
        m = FT(1 - 1 / n)
        S = FT((FT(1) + (α * (z - zm))^n)^(-m))
        return FT(S * porosity)
    end
    function soln(z, interface,porosity, n, α, δ, S_s)
        if z< interface
            return porosity +S_s*(interface-z)*exp(-δ*z)
        else
            return hydrostatic_profile(z,interface,porosity,n,α)
        end
    end
    
  #  plot(dons_arr[1]["soil.water.ϑ_l"],z, label = "initial state")
  #  plot!(dons_arr[2]["soil.water.ϑ_l"],z, label = "67 d")
  #  plot!(dons_arr[3]["soil.water.ϑ_l"],z, label = "133 d")
  #  plot!(dons_arr[4]["soil.water.ϑ_l"],z, label = "200 d")
  #  plot!(soln.(z,interface_z,0.495,vgn,2.6,0.2,1e-3),z, label = "steady state soln")
  #  plot!(legend = :bottomleft)
  #  plot!(xlabel = "ϑ_l")
  #  plot!(ylabel = "z")

    MSE = mean((soln.(z,interface_z,0.495,vgn,2.6,0.2,1e-3) .- dons_arr[4]["soil.water.ϑ_l"]) .^ 2.0)
    @test MSE < 1e-4
end


@testset "hydrostatic test 2" begin
    ClimateMachine.init()
    FT = Float64

    function init_soil_water!(land, state, aux, localgeo, time)
        myfloat = eltype(aux)
        state.soil.water.ϑ_l = myfloat(land.soil.water.initialϑ_l(aux))
        state.soil.water.θ_i = myfloat(land.soil.water.initialθ_i(aux))
    end

    soil_heat_model = PrescribedTemperatureModel()

    soil_param_functions = SoilParamFunctions{FT}(
        porosity = 0.41,
    )
    # Keep in mind that what is passed is aux⁻.
    # Scalar fluxes are multiplied by ẑ (normal to the surface, -normal to the bottom,
    # where normal point outs of the domain.)
    bottom_flux = (aux, t) -> eltype(aux)(0.0)
    surface_flux = bottom_flux
    ϑ_l0 = (aux) -> eltype(aux)(0.15)
    bc = LandDomainBC(
        bottom_bc = LandComponentBC(soil_water = Neumann(bottom_flux)),
        surface_bc = LandComponentBC(soil_water = Neumann(surface_flux)),
    )
    Ksat = (aux) -> eltype(aux)(0.443 / (3600 * 100))
    S_s = (aux) -> eltype(aux)(1e-3)
    sigmoid(x, offset, width) = typeof(x)(exp((x-offset)/width)/(1+exp((x-offset)/width)))

    vgα = (aux) -> eltype(aux)(sigmoid(aux.z, -0.5, 0.02)*(14.5-0.8)+0.8)
    vgn = (aux) -> eltype(aux)(sigmoid(aux.z, -0.5, 0.02)*(2.68-1.09)+1.09)
    
    soil_water_model = SoilWaterModel(
        FT;
        moisture_factor = MoistureDependent{FT}(),
        hydraulics = vanGenuchten(FT;n = vgn, α = vgα),
        Ksat = Ksat,
        S_s = S_s,
        initialϑ_l = ϑ_l0,
    )

    m_soil = SoilModel(soil_param_functions, soil_water_model, soil_heat_model)
    sources = ()
    m = LandModel(
        param_set,
        m_soil;
        boundary_conditions = bc,
        source = sources,
        init_state_prognostic = init_soil_water!,
    )


    N_poly = 1
    nelem_vert = 20

    # Specify the domain boundaries
    zmax = FT(0)
    zmin = FT(-1)

    driver_config = ClimateMachine.SingleStackConfiguration(
        "LandModel",
        N_poly,
        nelem_vert,
        zmax,
        param_set,
        m;
        zmin = zmin,
        numerical_flux_first_order = CentralNumericalFluxFirstOrder(),
    )

    t0 = FT(0)
    timeend = FT(60 * 60*24*50)

    dt = FT(100)
    n_outputs = 3
    every_x_simulation_time = ceil(Int, timeend / n_outputs)
    solver_config = ClimateMachine.SolverConfiguration(
        t0,
        timeend,
        driver_config,
        ode_dt = dt,
    )
    aux = solver_config.dg.state_auxiliary
    state_types = (Prognostic(),Auxiliary())
    dons_arr =
        Dict[dict_of_nodal_states(solver_config, state_types; interp = true)]
    time_data = FT[0]
    callback = GenericCallbacks.EveryXSimulationTime(every_x_simulation_time) do
        dons = dict_of_nodal_states(solver_config, state_types; interp = true)
        push!(dons_arr, dons)
        push!(time_data, gettime(solver_config.solver))
        nothing
    end
    ClimateMachine.invoke!(solver_config; user_callbacks = (callback,))
    z = dons_arr[1]["z"]

    function hydrostatic_profile(z, porosity, n, α)
        myf = eltype(z)
        m = FT(1 - 1 / n)
        S = FT((FT(1) + (α * abs(z))^n)^(-m))
        return FT(S * porosity)
    end
    function soln(z,porosity)
        if z< -0.5
            return hydrostatic_profile(z,porosity,1.09,0.8)
        else
            return hydrostatic_profile(z,porosity,2.68,14.5)
        end
    end
    
    plot(dons_arr[1]["soil.water.ϑ_l"],z, label = "initial state")
    plot!(dons_arr[2]["soil.water.ϑ_l"],z, label = "67 d")
    plot!(dons_arr[3]["soil.water.ϑ_l"],z, label = "133 d")
    plot!(dons_arr[4]["soil.water.ϑ_l"],z, label = "200 d")
    plot!(soln.(z,0.41),z, label = "steady state soln")
    plot!(legend = :bottomleft)
    plot!(xlabel = "ϑ_l")
    plot!(ylabel = "z")

    MSE = mean((soln.(z,interface_z,0.495,vgn,2.6,0.2,1e-3) .- dons_arr[4]["soil.water.ϑ_l"]) .^ 2.0)
    @test MSE < 1e-4
end
