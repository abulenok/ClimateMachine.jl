include("../PySDMCall/PySDMCall.jl")
using .PySDMCall


"""
    setup_dump_pysdm_vars(
        ::ClimateMachineConfigType,
        interval::String,
        out_prefix::String;
        writer = NetCDFWriter(),
        interpol = nothing,
    )

Create the "DumpPysdm" `DiagnosticsGroup` which contains all the
 state variables needed for pysdm. These are output on `x`, `y`, `z` or `lat`,
`long`, `level` dimensions of an interpolated grid (`interpol` _must_
be specified) as well as a `time` dimension at the specified `interval`.
"""
function setup_dump_pysdm_vars(
    ::ClimateMachineConfigType,
    interval::String,
    out_prefix::String;
    writer = NetCDFWriter(),
    interpol = nothing,
)
    # TODO: remove this
    @assert !isnothing(interpol)

    return DiagnosticsGroup(
        "DumpPysdm",
        Diagnostics.dump_pysdm_init,
        Diagnostics.dump_pysdm_fini,
        Diagnostics.dump_pysdm_collect,
        interval,
        out_prefix,
        writer,
        interpol,
    )
end

function dump_pysdm_init(dgngrp, currtime)
    FT = eltype(Settings.Q)
    bl = Settings.dg.balance_law
    mpicomm = Settings.mpicomm
    mpirank = MPI.Comm_rank(mpicomm)

    if mpirank == 0
        # get dimensions for the interpolated grid
        dims = dimensions(dgngrp.interpol)

        # set up the variables we're going to be writing
        vars = OrderedDict()
        #statenames = flattenednames(vars_state(bl, Prognostic(), FT))
        statenames = ["ρ", "ρu[1]", "ρu[2]", "ρu[3]", "moisture.θ_v"]
        
        for varname in statenames
            vars[varname] = (tuple(collect(keys(dims))...), FT, Dict())
        end


        dprefix = @sprintf("%s_%s", dgngrp.out_prefix, dgngrp.name)
        dfilename = joinpath(Settings.output_dir, dprefix)
        noov = Settings.no_overwrite
        init_data(dgngrp.writer, dfilename, noov, dims, vars)
    end

    return nothing
end

function dump_pysdm_collect(dgngrp, currtime) #TODO: make interpolation only for selected variables
    interpol = dgngrp.interpol
    mpicomm = Settings.mpicomm
    dg = Settings.dg
    Q = Settings.Q
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


    pysdm_vars = ["ρ", "ρu[1]", "ρu[3]", "moisture.θ_v"] # no "ρu[2]"


    if mpirank == 0
        statenames = flattenednames(vars_state(bl, Prognostic(), FT))
        auxnames = flattenednames(vars_state(bl, Auxiliary(), FT))
        
        println("PYSDM STATENAMES")
        println(statenames)

        println("PYSDM AUXNAMES")
        println(auxnames)

        varvals = OrderedDict()
        for (vari, varname) in enumerate(statenames)
            if varname in pysdm_vars
                println(varname)
                varvals[varname] = all_state_data[:, :, :, vari]
            end    
        end

        for (vari, varname) in enumerate(auxnames)
            if varname in pysdm_vars
                println(varname)
                varvals[varname] = all_aux_data[:, :, :, vari]
            end    
        end

        #ptest(varvals) # probably varvals should be converted to an array

        append_data(dgngrp.writer, varvals, currtime)
    end

    MPI.Barrier(mpicomm)
    return nothing
end

function dump_pysdm_fini(dgngrp, currtime) end
