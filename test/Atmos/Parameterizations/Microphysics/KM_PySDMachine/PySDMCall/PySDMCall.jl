
module PySDMCall

using PyCall
import Thermodynamics
const THDS = Thermodynamics
include("../../KinematicModel.jl")

export PySDM, PySDMConfig, PySDMCallWrapper, __init__


mutable struct PySDM
    config
    particulator
    rhod #TODO: Clima.ρ != PySDM.rhod
    exporter
end

mutable struct PySDMCallWrapper
    pysdm::PySDM
    init!
    do_step!
    fini!
end


function PySDMCallWrapper(pysdm_conf, init!, do_step!, fini!)
    return PySDMCallWrapper(
                PySDM(pysdm_conf, nothing, nothing, nothing),
                init!,
                do_step!,
                fini!
            )
end



mutable struct PySDMConfig
    grid
    size
    dxdz
    simtime
    dt
    n_sd
    kappa # hygroscopicity
    kernel # Geometric from PySDM
    spectrum_per_mass_of_dry_air
end

function PySDMConfig(
    size,
    dxdz,
    simtime,
    dt,
    n_sd_per_gridbox,
    kappa,
    kernel,
    spectrum_per_mass_of_dry_air
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
        spectrum_per_mass_of_dry_air
    )
end


# CMStepper # coś co ma metodę wait & step # https://github.com/atmos-cloud-sim-uj/PySDM-examples/blob/main/PySDM_examples/Arabas_et_al_2015/mpdata.py



function __init__()
    # adds directories to the Python search path.
    pushfirst!(PyVector(pyimport("sys")."path"), "")
    pushfirst!(PyVector(pyimport("sys")."path"), "PySDMCall/")
    pushfirst!(PyVector(pyimport("sys")."path"), "test/Atmos/Parameterizations/Microphysics/KM_PySDMachine/PySDMCall/")
end


end # module PySDMCall