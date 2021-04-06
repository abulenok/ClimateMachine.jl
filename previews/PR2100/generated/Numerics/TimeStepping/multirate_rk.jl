using ClimateMachine
const clima_dir = dirname(dirname(pathof(ClimateMachine)));
include(joinpath(
    clima_dir,
    "tutorials",
    "Numerics",
    "TimeStepping",
    "tutorial_acousticwave_config.jl",
))

FT = Float64;

ode_solver = ClimateMachine.MultirateSolverType(
    splitting_type = ClimateMachine.HEVISplitting(),
    slow_method = LSRK54CarpenterKennedy,
    fast_method = ARK2GiraldoKellyConstantinescu,
    implicit_solver_adjustable = true,
    timestep_ratio = 100,
)

timeend = FT(3600)
CFL = FT(5)
cfl_direction = HorizontalDirection()
run_acousticwave(ode_solver, CFL, cfl_direction, timeend);

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

