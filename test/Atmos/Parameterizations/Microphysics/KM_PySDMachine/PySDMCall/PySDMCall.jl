
module PySDMCall

using PyCall
import Thermodynamics
const THDS = Thermodynamics
include("../../KinematicModel.jl")

export PySDM, PySDMConfig, PySDMCallWrapper, __init__


mutable struct PySDMConfig
    grid::Tuple{Int64, Int64}
    size::Tuple{Int64, Int64}
    dxdz::Tuple{Float64, Float64}
    simtime::Float64
    dt::Float64
    n_sd::Int64
    kappa::Int64 # hygroscopicity
    kernel::PyCall.PyObject # from PySDM
    spectrum_per_mass_of_dry_air::PyCall.PyObject
end

function PySDMConfig(
    size::Tuple{Int64, Int64},
    dxdz::Tuple{Float64, Float64},
    simtime::Float64,
    dt::Float64,
    n_sd_per_gridbox::Int64,
    kappa::Int64,
    kernel::Any,
    spectrum_per_mass_of_dry_air::Any,
)
    grid = (Int(size[1] / dxdz[1]), Int(size[2] / dxdz[2]))

    n_sd = grid[1] * grid[2] * n_sd_per_gridbox

    PySDMConfig(
        grid,
        size,
        dxdz,
        simtime,
        dt,
        n_sd,
        kappa,
        kernel,
        spectrum_per_mass_of_dry_air,
    )
end

mutable struct PySDM
    config::PySDMConfig
    particulator::Any
    rhod::Any #TODO: Clima.ρ != PySDM.rhod
    exporter::Any
end

mutable struct PySDMCallWrapper
    pysdm::PySDM
    init!::Any
    do_step!::Any
    fini!::Any

    function PySDMCallWrapper(pysdm_conf::PySDMConfig, init!, do_step!, fini!)

        return new(
            PySDM(pysdm_conf, nothing, nothing, nothing),
            init!,
            do_step!,
            fini!,
        )
    end
end


# CMStepper # coś co ma metodę wait & step # https://github.com/atmos-cloud-sim-uj/PySDM-examples/blob/main/PySDM_examples/Arabas_et_al_2015/mpdata.py

function __init__()
    # adds directories to the Python search path.
    pushfirst!(PyVector(pyimport("sys")."path"), "")
    pushfirst!(PyVector(pyimport("sys")."path"), "PySDMCall/")
    pushfirst!(
        PyVector(pyimport("sys")."path"),
        "test/Atmos/Parameterizations/Microphysics/KM_PySDMachine/PySDMCall/",
    )
end

end # module PySDMCall
