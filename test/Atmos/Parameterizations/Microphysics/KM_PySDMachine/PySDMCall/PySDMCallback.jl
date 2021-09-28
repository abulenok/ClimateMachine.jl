module PySDMCallbacks

export PySDMCallback

include("PySDMCall.jl")
include("../../KinematicModel.jl")

using .PySDMCall

using MPI
using OrderedCollections
using PyCall

using ClimateMachine.BalanceLaws
using ClimateMachine.Mesh.Interpolation
using ClimateMachine.DGMethods: SpaceDiscretization

import Thermodynamics
const THDS = Thermodynamics
import ClimateMachine.GenericCallbacks


mutable struct PySDMCallback
    name::String
    dg::SpaceDiscretization
    interpol::InterpolationTopology
    mpicomm::MPI.Comm
    pysdmcw::PySDMCallWrapper
end


function GenericCallbacks.init!(cb::PySDMCallback, solver, Q, param, t)
    println()
    println("PySDMCallback init")
    println(t)
    println()

    varvals = vals_interpol(cb, Q)

    if !isnothing(cb.pysdmcw.init!)
        cb.pysdmcw.init!(cb.pysdmcw.pysdm, varvals)
    end

    return nothing
end

function GenericCallbacks.call!(cb::PySDMCallback, solver, Q, param, t)
    println()
    println("PySDMCallback call")
    println(t)
    println()

    vals = vals_interpol(cb, Q)

    if !isnothing(cb.pysdmcw.do_step!)
        cb.pysdmcw.do_step!(cb.pysdmcw.pysdm, vals, t)
    end

    return nothing
end

function GenericCallbacks.fini!(cb::PySDMCallback, solver, Q, param, t)
    println()
    println("PySDMCallback fini")
    println(t)
    println()

    if !isnothing(cb.pysdmcw.fini!)
        cb.pysdmcw.fini!(cb.pysdmcw.pysdm, vals, t)
    end

    return nothing
end



function vals_interpol(cb::PySDMCallback, Q)

    interpol = cb.interpol
    mpicomm = cb.mpicomm
    dg = cb.dg
    FT = eltype(Q.data)
    bl = dg.balance_law
    mpirank = MPI.Comm_rank(mpicomm)

    istate = similar(Q.data, interpol.Npl, number_states(bl, Prognostic()))

    interpolate_local!(interpol, Q.data, istate)

    if interpol isa InterpolationCubedSphere
        # TODO: get indices here without hard-coding them
        _ρu, _ρv, _ρw = 2, 3, 4
        project_cubed_sphere!(interpol, istate, (_ρu, _ρv, _ρw))
    end

    iaux = similar(
        dg.state_auxiliary.data,
        interpol.Npl,
        number_states(bl, Auxiliary()),
    )

    interpolate_local!(interpol, dg.state_auxiliary.data, iaux)

    all_state_data = accumulate_interpolated_data(mpicomm, interpol, istate)
    all_aux_data = accumulate_interpolated_data(mpicomm, interpol, iaux)

    pysdm_vars = ["ρ", "ρu[1]", "ρu[3]", "q_vap", "theta_dry", "q_tot", "e_int"]


    varvals = nothing

    if mpirank == 0
        statenames = flattenednames(vars_state(bl, Prognostic(), FT))
        auxnames = flattenednames(vars_state(bl, Auxiliary(), FT))

        varvals = OrderedDict()
        for (vari, varname) in enumerate(statenames)
            if varname in pysdm_vars

                varvals[varname] = all_state_data[:, :, :, vari]

            end
        end

        for (vari, varname) in enumerate(auxnames)
            if varname in pysdm_vars

                varvals[varname] = all_aux_data[:, :, :, vari]

            end
        end

    end

    MPI.Barrier(mpicomm)

    return varvals
end

end #module