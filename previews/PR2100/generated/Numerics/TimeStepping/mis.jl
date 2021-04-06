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

ode_solver = ClimateMachine.MISSolverType(;
    mis_method = MIS2,
    fast_method = LSRK144NiegemannDiehlBusch,
    nsubsteps = (40,),
)

timeend = FT(500)
CFL = FT(20)
run_simulation(ode_solver, CFL, timeend);

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

