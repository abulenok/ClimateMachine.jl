using ClimateMachine
const clima_dir = dirname(dirname(pathof(ClimateMachine)));
include(joinpath(
    clima_dir,
    "tutorials",
    "Numerics",
    "TimeStepping",
    "tutorial_risingbubble_config.jl",
))

FT = Float64;

timeend = FT(100)
ode_solver =
    ClimateMachine.ExplicitSolverType(solver_method = LSRK54CarpenterKennedy)
CFL = FT(0.4)
run_simulation(ode_solver, CFL, timeend);

ode_solver = ClimateMachine.ExplicitSolverType(
    solver_method = LSRK144NiegemannDiehlBusch,
)
CFL = FT(1.7)
run_simulation(ode_solver, CFL, timeend);

ode_solver = ClimateMachine.ExplicitSolverType(solver_method = SSPRK33ShuOsher)
CFL = FT(0.2)
run_simulation(ode_solver, CFL, timeend);

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

