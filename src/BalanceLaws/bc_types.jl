##### Boundary condition types

using DispatchedTuples

export BCDef, DefaultBC, set_bcs!, dispatched_tuple

"""
    BCDef

Subtypes are used for specifying
each boundary condition.
"""
abstract type BCDef{PV <: PrognosticVariable} end

"""
    DefaultBC

The default boundary condition
"""
struct DefaultBC{PV} <: BCDef{PV} end

"""
    bc_val(::BCDef{PV}, bl, nf, args)

Return the value of the boundary condition, given
 - `bcd::BCDef` the boundary condition definition type
 - `bl` the balance law
 - `nf` the numerical flux
 - `args` top-level arguments, packed into a NamedTuple
"""
bc_val(::BCDef{PV}, bl, nf, args) where {PV} = DefaultBC{PV}()

"""
    default_bcs(::BalanceLaw)

A tuple of BalanceLaw-specific default boundary conditions
"""
function default_bcs end

"""
    set_bcs!(state⁺, bl, nf, bc, ntargs, prog_vars = prognostic_vars(bl))

A convenience method for setting `state⁺` inside `boundary_state!`.

Arguments:
 - `state⁺` the exterior state
 - `bl` the balance law
 - `nf::Union{NumericalFluxFirstOrder,NumericalFluxGradient}` the numerical flux
 - `ntargs` the top-level arguments
 - `prog_vars` (optional) the balance law's prognostic variables
"""
function set_bcs!(state⁺, bl, nf, bc, ntargs, prog_vars = prognostic_vars(bl))
    state⁻ = ntargs.state
    bcs_default = dispatched_tuple(default_bcs(bl))
    map(prog_vars) do prog
        var⁺, name = get_prog_state(state⁺, prog)
        var⁻, name = get_prog_state(state⁻, prog)
        tup = used_bcs(bl, prog, bc, bcs_default)
        bcvals = map(tup) do bc_pv
            bc_val(bc_pv, bl, nf, ntargs)
        end
        set_bc!(var⁺, name, bcvals)
    end
end

"""
    dispatched_tuple(bc_defs::Tuple)

A DispatchedTuple, given a tuple of [`BCDef`](@ref)'s.
"""
function dispatched_tuple(bc_defs::Tuple)
    tup = map(bc_defs) do e
        pvi = prog_var_instance(e)
        @assert pvi isa PrognosticVariable
        (pvi, e)
    end
    return DispatchedTuple(tup, nothing)
end

#####
##### Internal methods for setting/diagonalizing boundary conditions
#####

driver_prescribed(bc) = bc.tup
# TODO: Simplify logic
function used_bcs(bl::BalanceLaw, prog::PrognosticVariable, bc, bcs_default)
    bcs_driver_prescribed = driver_prescribed(bc)
    bcs_prescribed_prog = dispatch(bcs_driver_prescribed, prog)
    # TODO: Do we need to check bcs_prescribed_prog == ((),) ?
    if bcs_prescribed_prog == ((),) || bcs_prescribed_prog == (nothing,)
        # use bl-specific BCs:
        tup = dispatch(bcs_default, prog)
    else # use driver-prescribed BCs:
        tup = bcs_prescribed_prog
    end

    # Error checking:
    isa_prog(pv::Type{PV}) where {PV <: PrognosticVariable} = true
    isa_prog(pv::Type{PV}) where {PV} = false
    @assert all(isa_prog.(type_param.(tup)))
    return tup
end

diag_bc(::PVA, ::BCDef{PVB}) where {PVA, PVB} = nothing
diag_bc(::PV, bcd::BCDef{PV}) where {PV} = bcd

filter_bcs(t::Tuple) = filter(x -> !(x == nothing), t)
diag_bc(tup, pv) = filter_bcs(map(bc -> diag_bc(pv, bc), tup))

prog_var_bcs(::Tuple{}, ::PV) where {PV} = (DefaultBC{PV}(),)
prog_var_bcs(t::Tuple, ::PV) where {PV} = t

bcs_per_prog_var(tup, pv::PrognosticVariable) =
    prog_var_bcs(diag_bc(tup, pv), pv)

set_bc!(var⁺, name, bcvals::Tuple{DefaultBC}) = nothing
set_bc!(var⁺, name, bcvals::NTuple{N, DefaultBC}) where {N} = nothing
set_bc!(var⁺, name, bcvals) = setproperty!(var⁺, name, sum(bcvals))

prog_var_instance(bc::BCDef{PV}) where {PV} = PV()
