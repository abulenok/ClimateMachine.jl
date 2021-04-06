using ClimateMachine
const clima_dir = dirname(dirname(pathof(ClimateMachine)));
include(joinpath(
    clima_dir,
    "tutorials",
    "Numerics",
    "TimeStepping",
    "tutorial_acousticwave_config.jl",
));

FT = Float64
timeend = FT(100)

ode_solver = ClimateMachine.ExplicitSolverType(
    solver_method = LSRK144NiegemannDiehlBusch,
);

CFL = FT(0.002)
cfl_direction = HorizontalDirection()
run_acousticwave(ode_solver, CFL, cfl_direction, timeend);

timeend = FT(3600)
ode_solver = ClimateMachine.IMEXSolverType(
    solver_method = ARK2GiraldoKellyConstantinescu,
)
CFL = FT(0.5)
cfl_direction = HorizontalDirection()
run_acousticwave(ode_solver, CFL, cfl_direction, timeend);

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

