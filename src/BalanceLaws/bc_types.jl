##### Boundary condition types

import Base
using DispatchedTuples

export DefaultBC, DefaultBCValue, DefaultBCFlux
export set_boundary_values!
export set_boundary_fluxes!

# The following BC naming convention is used:

# `all_bcs` - all boundary conditions for a given balance
#             law (on all faces of the geometry) This is
#             typically a tuple of individual BCs.
#
# `bc_set` - a set of boundary conditions for a single point
#            in space, or wall, for all variables.
#
# `bc`     - a single boundary condition for a single point
#            in space, or wall, for a single variable.

# Additional suffixes are used to indicate one of several things:
# - `driver` - BCs described in the driver
# - `default` - balance law-specific BCs

"""
    DefaultBCValue

Default BC value, which results in the default
boundary condition behavior.
"""
struct DefaultBCValue end

"""
    DefaultBCFlux

Default BC flux, which results in the default
boundary condition behavior.
"""
struct DefaultBCFlux end

"""
    DefaultBC

The default boundary condition definition,
which results in yielding the default boundary
condition value [`DefaultBCValue`](@ref).
"""
struct DefaultBC end

# If DefaultBCValue's are mixed with real values,
# then let DefaultBCValue's be zero:
Base.:+(::DefaultBCValue, x::FT) where {FT <: AbstractFloat} = x
Base.:+(x::FT, ::DefaultBCValue) where {FT <: AbstractFloat} = x
Base.:+(::DefaultBCValue, x::SArray{Tuple{3}, FT, 1, 3}) where {FT} = x
Base.:+(x::SArray{Tuple{3}, FT, 1, 3}, ::DefaultBCValue) where {FT} = x

Base.:+(::DefaultBCFlux, x::FT) where {FT <: AbstractFloat} = x
Base.:+(x::FT, ::DefaultBCFlux) where {FT <: AbstractFloat} = x
Base.:+(::DefaultBCFlux, x::SArray{Tuple{3}, FT, 1, 3}) where {FT} = x
Base.:+(x::SArray{Tuple{3}, FT, 1, 3}, ::DefaultBCFlux) where {FT} = x

"""
    boundary_value(::PrognosticVariable, bc, bl, args, nf)

Return the value of the boundary condition, given
 - `::PrognosticVariable` the prognostic variable
 - `bc` the individual boundary condition type
 - `bl` the balance law
 - `args` top-level arguments, packed into a NamedTuple
 - `nf` the numerical flux
"""
boundary_value(::PrognosticVariable, ::DefaultBC, bl, args, nf) =
    DefaultBCValue()

"""
    boundary_flux(::PrognosticVariable, bc, bl, args, nf)

Return the prescribed boundary flux, given
 - `::PrognosticVariable` the prognostic variable
 - `bc` the individual boundary condition type
 - `bl` the balance law
 - `args` top-level arguments, packed into a NamedTuple
 - `nf` the numerical flux
"""
boundary_flux(::PrognosticVariable, ::DefaultBC, bl, args, nf) = DefaultBCFlux()

"""
    default_bcs(::PrognosticVariable)

A tuple of default boundary condition definitions
for a given prognostic variable.
"""
function default_bcs end

"""
    set_boundary_values!(
        state⁺,
        bl,
        nf,
        bc_set_driver,
        args,
        prog_vars = prognostic_vars(bl)
    )

A convenience method for setting `state⁺` inside `boundary_state!`.

Arguments:
 - `state⁺` the exterior state
 - `bl` the balance law
 - `bc_set_driver` the balance law boundary conditions on a given boundary
 - `nf::Union{NumericalFluxFirstOrder,NumericalFluxGradient}` the numerical flux
 - `args` the top-level arguments
 - `prog_vars` (optional) the balance law's prognostic variables
"""
function set_boundary_values!(
    state⁺,
    bl,
    nf,
    bc_set_driver,
    args,
    prog_vars = prognostic_vars(bl),
)
    state⁻ = args.state⁻
    bc_set_default = default_bcs(bl)
    map(prog_vars) do pv
        var⁺, name = get_prog_state(state⁺, pv)
        var⁻, name = get_prog_state(state⁻, pv)
        tup = used_bcs(bl, pv, bc_set_driver, bc_set_default)
        bcvals = map(tup) do bc
            boundary_value(pv, bc, bl, args, nf)
        end
        set_bc!(var⁺, name, bcvals)
    end
end

"""
    set_boundary_fluxes!(
        fluxᵀn,
        bl,
        nf,
        bc_set_driver,
        args,
        prog_vars = prognostic_vars(bl)
    )

A convenience method for setting the
numerical (boundary) flux `fluxᵀn` in
`normal_boundary_flux_second_order!`.

Arguments:
 - `fluxᵀn` the numerical (boundary) flux
 - `bl` the balance law
 - `bc_set_driver` the balance law boundary conditions on a given boundary
 - `nf::NumericalFluxSecondOrder` the numerical flux
 - `args` the top-level arguments
 - `prog_vars` (optional) the balance law's prognostic variables
"""
function set_boundary_fluxes!(
    fluxᵀn,
    bl,
    nf,
    bc_set_driver,
    args,
    prog_vars = prognostic_vars(bl),
)
    bc_set_default = default_bcs(bl)
    map(prog_vars) do pv
        varflux, name = get_prog_state(fluxᵀn, pv)
        tup = used_bcs(bl, pv, bc_set_driver, bc_set_default)
        bfluxes = map(tup) do bc
            boundary_flux(pv, bc, bl, args, nf)
        end
        Σbfluxes = sum_boundary_fluxes(bfluxes)
        varprop = getproperty(varflux, name)
        varflux_assigned = Σbfluxes + varprop
        setproperty!(varflux, name, varflux_assigned)
    end
end

#####
##### Internal methods for setting/diagonalizing boundary conditions
#####

sum_boundary_fluxes(bcvals::NTuple{N, DefaultBCFlux}) where {N} =
    DefaultBCFlux()
sum_boundary_fluxes(bcvals) = sum(bcvals)

# Internal method: call default_bcs for bl-specific
# prognostic variables:
function default_bcs(bl::BalanceLaw)
    tup = map(prognostic_vars(bl)) do pv
        map(default_bcs(pv)) do bc
            @assert pv isa PrognosticVariable
            (pv, bc)
        end
    end
    tup = tuple_of_tuples(tup)
    return DispatchedTuple(tup, DefaultBC())
end

# Flatten "tuple of tuple of tuples" to "tuple of tuples"
tuple_of_tuples(a::Tuple{PrognosticVariable, T}) where {T} = (a,)
tuple_of_tuples(a, b...) =
    tuple(tuple_of_tuples(a)..., tuple_of_tuples(b...)...)
tuple_of_tuples(a::Tuple) = tuple_of_tuples(a...)

# Precedence convention: Choose
# driver-prescribed BCs over
# bl-specific default BCs
function used_bcs(
    bl::BalanceLaw,
    pv::PrognosticVariable,
    bc_set_driver,
    bc_set_default,
)
    bc_driver = dispatch(bc_set_driver.tup, pv)
    if bc_driver == (DefaultBC(),) # use bl-specific BCs:
        bc_used = dispatch(bc_set_default, pv)
    else # use driver-prescribed BCs:
        bc_used = bc_driver
    end
    return bc_used
end

# TODO: remove Tuple{DefaultBC} method:
# TODO: why is NTuple{N, DefaultBC}) where {N} needed?
set_bc!(var⁺, name, bcvals::Tuple{DefaultBCValue}) = nothing
set_bc!(var⁺, name, bcvals::NTuple{N, DefaultBCValue}) where {N} = nothing
set_bc!(var⁺, name, bcvals) = setproperty!(var⁺, name, sum(bcvals))

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

# TODO: Cannot currently pre-compute used bcs for all prognostic_vars
#       Alternatively, we could use `Set()` inside `bc_precompute`
# function used_bcs(bl::BalanceLaw, bc_set_driver)
#     bc_set_default = default_bcs(bl)
#     bc_set_used = map(prognostic_vars(bl)) do pv
#         used_bcs(bl, pv, bc_set_driver, bc_set_default)
#     end
#     return Tuple(Set(bc_set_used))
# end

# function bc_precompute(bl, bc_set_driver)
#     bc_set_default = default_bcs(bl)

#     tup = map(prognostic_vars(bl)) do pv
#         bcs_used = used_bcs(bl, pv, bc_set_driver, bc_set_default)
#         map(bcs_used) do bc
#             (bc, bc_precompute(bc, atmos, args, nf))
#         end
#     end
#     dtup = DispatchedTuple(tuple_of_tuples(tup))
#     # dtup = DispatchedTupleDict(tuple_of_tuples(tup))
#     return (; dtup)
# end
