##### Boundary condition types

import Base
using DispatchedTuples

export BCDef, DefaultBC, set_bcs!, dispatched_tuple
export DefaultBCValue

"""
    BCDef

Subtypes are used for specifying
each boundary condition.
"""
abstract type BCDef{PV <: PrognosticVariable} end

"""
    DefaultBCValue

Default BC value, which results in the default
boundary condition behavior.
"""
struct DefaultBCValue end

"""
    DefaultBC{PV} <: BCDef{PV}

The default boundary condition definition,
which results in yielding the default boundary
condition value [`DefaultBCValue`](@ref).
"""
struct DefaultBC{PV} <: BCDef{PV} end


# If DefaultBCValue's are mixed with real values,
# then let DefaultBCValue's be zero:
Base.:+(::DefaultBCValue, x::FT) where {FT <: AbstractFloat} = x
Base.:+(x::FT, ::DefaultBCValue) where {FT <: AbstractFloat} = x

"""
    bc_val(::BCDef{PV}, bl, nf, args)

Return the value of the boundary condition, given
 - `bcd::BCDef` the boundary condition definition type
 - `bl` the balance law
 - `nf` the numerical flux
 - `args` top-level arguments, packed into a NamedTuple
"""
bc_val(::DefaultBC{PV}, bl, nf, args) where {PV} = DefaultBCValue()

"""
    default_bcs(::BalanceLaw)

A tuple of BalanceLaw-specific default boundary conditions
"""
function default_bcs end

"""
    set_bcs!(state⁺, bl, nf, bc, args, prog_vars = prognostic_vars(bl))

A convenience method for setting `state⁺` inside `boundary_state!`.

Arguments:
 - `state⁺` the exterior state
 - `bl` the balance law
 - `nf::Union{NumericalFluxFirstOrder,NumericalFluxGradient}` the numerical flux
 - `args` the top-level arguments
 - `prog_vars` (optional) the balance law's prognostic variables
"""
function set_bcs!(state⁺, bl, nf, bc, args, prog_vars = prognostic_vars(bl))
    state⁻ = args.state
    bcs_default = dispatched_tuple(default_bcs(bl))
    map(prog_vars) do prog
        var⁺, name = get_prog_state(state⁺, prog)
        var⁻, name = get_prog_state(state⁻, prog)
        tup = used_bcs(bl, prog, bc, bcs_default)
        bcvals = map(tup) do bc_pv
            bc_val(bc_pv, bl, nf, args)
        end
        Σbcvals = sum_bc_vals(bcvals)
        setproperty!(var⁺, name, Σbcvals)
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
    # isa_prog(pv::Type{PV}) where {PV <: PrognosticVariable} = true
    # isa_prog(pv::Type{PV}) where {PV} = false
    # @assert all(isa_prog.(type_param.(tup)))
    return tup
end

sum_bc_vals(bcvals::NTuple{N, DefaultBCValue}) where {N} = 0
sum_bc_vals(bcvals) = sum(bcvals)

prog_var_instance(bc::BCDef{PV}) where {PV} = PV()

"""
    bc_precompute(bc, bl, args, dispatch_helper)

A nested NamedTuple of precomputed (cached) values
and or objects. This is useful for "expensive"
point-wise quantities that are used in multiple
boundary condition functions. For example, computing
a quantity that requires iteration.

This is a separated from [`precompute`](@ref), as
simply because there are many `precompute` definitions,
and splitting the methods helps search-ability.
"""
bc_precompute(bc, bl, args, dispatch_helper) = NamedTuple()
