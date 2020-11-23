module SingleStackUtils

export get_vars_from_nodal_stack,
    get_vars_from_element_stack,
    get_horizontal_variance,
    get_horizontal_mean,
    reduce_nodal_stack,
    reduce_element_stack,
    horizontally_average!,
    max_horizontal_invariance,
    assert_horizontally_uniform,
    dict_of_nodal_states

using OrderedCollections
using StaticArrays
import KernelAbstractions: CPU

using ..BalanceLaws
using ..DGMethods
using ..DGMethods.Grids
using ..MPIStateArrays
using ..VariableTemplates

"""
    get_vars_from_nodal_stack(
        grid::DiscontinuousSpectralElementGrid{T, dim, Ns},
        Q::MPIStateArray,
        vars;
        vrange::UnitRange = 1:size(Q, 3),
        i::Int = 1,
        j::Int = 1,
        exclude::Vector{String} = String[],
        interp = false,
    ) where {T, dim, Ns}

Return a dictionary whose keys are the `flattenednames()` of the variables
specified in `vars` (as returned by e.g. `vars_state`), and
whose values are arrays of the values for that variable along the vertical
dimension in `Q`. Only a single element is expected in the horizontal as
this is intended for the single stack configuration and `i` and `j` identify
the horizontal nodal coordinates.

Variables listed in `exclude` are skipped.
"""
function get_vars_from_nodal_stack(
    grid::DiscontinuousSpectralElementGrid{T, dim, Ns},
    Q::MPIStateArray,
    vars;
    vrange::UnitRange = 1:size(Q, 3),
    i::Int = 1,
    j::Int = 1,
    exclude::Vector{String} = String[],
    interp = false,
) where {T, dim, Ns}

    # XXX: Needs updating for multiple polynomial orders
    # Currently only support single polynomial order
    @assert all(Ns[1] .== Ns)
    N = Ns[1]

    # extract grid information and bring `Q` to the host if needed
    FT = eltype(Q)
    Nq = N + 1
    Nqk = dimensionality(grid) == 2 ? 1 : Nq
    Np = dofs_per_element(grid)
    state_data = array_device(Q) isa CPU ? Q.realdata : Array(Q.realdata)

    # set up the dictionary to be returned
    var_names = flattenednames(vars)
    stack_vals = OrderedDict()
    num_vars = varsize(vars)
    vars_wanted = Int[]
    @inbounds for vi in 1:num_vars
        if !(var_names[vi] in exclude)
            stack_vals[var_names[vi]] = FT[]
            push!(vars_wanted, vi)
        end
    end
    vmap⁻ = array_device(Q) isa CPU ? grid.vmap⁻ : Array(grid.vmap⁻)
    vmap⁺ = array_device(Q) isa CPU ? grid.vmap⁺ : Array(grid.vmap⁺)
    vgeo = array_device(Q) isa CPU ? grid.vgeo : Array(grid.vgeo)
    # extract values from `state_data`
    @inbounds for ev in vrange, k in 1:Nqk, v in vars_wanted
        if interp && k == 1
            # Get face degree of freedom number
            n = i + Nq * ((j - 1))
            # get the element numbers
            ev⁻ = ev
            # Get neighboring id data
            id⁻, id⁺ = vmap⁻[n, 5, ev⁻], vmap⁺[n, 5, ev⁻]
            ev⁺ = ((id⁺ - 1) ÷ Np) + 1
            # get the volume degree of freedom numbers
            vid⁻, vid⁺ = ((id⁻ - 1) % Np) + 1, ((id⁺ - 1) % Np) + 1

            J⁺, J⁻ = vgeo[vid⁻, Grids._M, ev⁻], vgeo[vid⁺, Grids._M, ev⁺]
            state_local = J⁻ * state_data[vid⁻, v, ev⁻]
            state_local += J⁺ * state_data[vid⁺, v, ev⁺]
            state_local /= (J⁻ + J⁺)
            push!(stack_vals[var_names[v]], state_local)
        elseif interp && k == Nqk
            # Get face degree of freedom number
            n = i + Nq * ((j - 1))
            # get the element numbers
            ev⁻ = ev
            # Get neighboring id data
            id⁻, id⁺ = vmap⁻[n, 6, ev⁻], vmap⁺[n, 6, ev⁻]
            # periodic and need to handle this point (otherwise handled above)
            if id⁺ == id⁻
                vid⁻ = ((id⁻ - 1) % Np) + 1

                state_local = state_data[vid⁻, v, ev⁻]
                push!(stack_vals[var_names[v]], state_local)
            end
        else
            ijk = i + Nq * ((j - 1) + Nq * (k - 1))
            state_local = state_data[ijk, v, ev]
            push!(stack_vals[var_names[v]], state_local)
        end
    end

    return stack_vals
end

"""
    get_vars_from_element_stack(
        grid::DiscontinuousSpectralElementGrid{T, dim, Ns},
        Q::MPIStateArray,
        vars;
        vrange::UnitRange = 1:size(Q, 3),
        exclude::Vector{String} = String[],
        interp = false,
    ) where {T, dim, Ns}

Return an array of [`get_vars_from_nodal_stack()`](@ref)s whose dimensions
are the number of nodal points per element in the horizontal plane.

Variables listed in `exclude` are skipped.
"""
function get_vars_from_element_stack(
    grid::DiscontinuousSpectralElementGrid{T, dim, Ns},
    Q::MPIStateArray,
    vars;
    vrange::UnitRange = 1:size(Q, 3),
    exclude::Vector{String} = String[],
    interp = false,
) where {T, dim, Ns}

    # XXX: Needs updating for multiple polynomial orders
    # Currently only support single polynomial order
    @assert all(Ns[1] .== Ns)
    N = Ns[1]

    Nq = N + 1
    return [
        get_vars_from_nodal_stack(
            grid,
            Q,
            vars,
            vrange = vrange,
            i = i,
            j = j,
            exclude = exclude,
            interp = interp,
        ) for i in 1:Nq, j in 1:Nq
    ]
end

"""
    get_horizontal_mean(
        grid::DiscontinuousSpectralElementGrid{T, dim, Ns},
        Q::MPIStateArray,
        vars;
        vrange::UnitRange = 1:size(Q, 3),
        exclude::Vector{String} = String[],
        interp = false,
    ) where {T, dim, Ns}

Return a dictionary whose keys are the `flattenednames()` of the variables
specified in `vars` (as returned by e.g. `vars_state`), and
whose values are arrays of the horizontal averages for that variable along
the vertical dimension in `Q`. Only a single element is expected in the
horizontal as this is intended for the single stack configuration.

Variables listed in `exclude` are skipped.
"""
function get_horizontal_mean(
    grid::DiscontinuousSpectralElementGrid{T, dim, Ns},
    Q::MPIStateArray,
    vars;
    vrange::UnitRange = 1:size(Q, 3),
    exclude::Vector{String} = String[],
    interp = false,
) where {T, dim, Ns}

    # XXX: Needs updating for multiple polynomial orders
    # Currently only support single polynomial order
    @assert all(Ns[1] .== Ns)
    N = Ns[1]

    Nq = N + 1
    vars_avg = OrderedDict()
    vars_sq = OrderedDict()
    for i in 1:Nq
        for j in 1:Nq
            vars_nodal = get_vars_from_nodal_stack(
                grid,
                Q,
                vars,
                vrange = vrange,
                i = i,
                j = j,
                exclude = exclude,
                interp = interp,
            )
            vars_avg = merge(+, vars_avg, vars_nodal)
        end
    end
    map!(x -> x ./ Nq / Nq, values(vars_avg))
    return vars_avg
end

"""
    get_horizontal_variance(
        grid::DiscontinuousSpectralElementGrid{T, dim, Ns},
        Q::MPIStateArray,
        vars;
        vrange::UnitRange = 1:size(Q, 3),
        exclude::Vector{String} = String[],
        interp = false,
    ) where {T, dim, Ns}

Return a dictionary whose keys are the `flattenednames()` of the variables
specified in `vars` (as returned by e.g. `vars_state`), and
whose values are arrays of the horizontal variance for that variable along
the vertical dimension in `Q`. Only a single element is expected in the
horizontal as this is intended for the single stack configuration.

Variables listed in `exclude` are skipped.
"""
function get_horizontal_variance(
    grid::DiscontinuousSpectralElementGrid{T, dim, Ns},
    Q::MPIStateArray,
    vars;
    vrange::UnitRange = 1:size(Q, 3),
    exclude::Vector{String} = String[],
    interp = false,
) where {T, dim, Ns}

    # XXX: Needs updating for multiple polynomial orders
    # Currently only support single polynomial order
    @assert all(Ns[1] .== Ns)
    N = Ns[1]

    Nq = N + 1
    vars_avg = OrderedDict()
    vars_sq = OrderedDict()
    for i in 1:Nq
        for j in 1:Nq
            vars_nodal = get_vars_from_nodal_stack(
                grid,
                Q,
                vars,
                vrange = vrange,
                i = i,
                j = j,
                exclude = exclude,
                interp = interp,
            )
            vars_nodal_sq = OrderedDict(vars_nodal)
            map!(x -> x .^ 2, values(vars_nodal_sq))
            vars_avg = merge(+, vars_avg, vars_nodal)
            vars_sq = merge(+, vars_sq, vars_nodal_sq)
        end
    end
    map!(x -> (x ./ Nq / Nq) .^ 2, values(vars_avg))
    map!(x -> x ./ Nq / Nq, values(vars_sq))
    vars_var = merge(-, vars_sq, vars_avg)
    return vars_var
end

"""
    max_horizontal_invariance(
        grid::DiscontinuousSpectralElementGrid{T, dim, Ns},
        Q::MPIStateArray,
        vars;
        exclude::Vector{String} = String[],
    ) where {T, dim, Ns}

Max horizontal variance of variables specified in `vars` (as returned by
e.g. `vars_state`) are not horizontally uniform along the vertical stack.

Variables listed in `exclude` are skipped.
"""
function max_horizontal_invariance(
    grid::DiscontinuousSpectralElementGrid{T, dim, Ns},
    Q::MPIStateArray,
    vars;
    exclude::Vector{String} = String[],
) where {T, dim, Ns}
    horiz_var = get_horizontal_variance(
        grid,
        Q,
        vars;
        exclude = exclude,
    )
    max_horiz_var = Dict(var => maximum(abs.(hv)) for (var,hv) in horiz_var)
    return max_horiz_var
end

"""
    assert_horizontally_uniform(max_horiz_var, var, FT, tol = 10eps(FT))

Error if field `var`, computed from `max_horizontal_invariance`,
is horizontally non-uniform above tolerance `tol`.
"""
function assert_horizontally_uniform(max_horiz_var, var, FT, tol = 10eps(FT))
    if !(max_horiz_var[var] < tol)
        @show max_horiz_var[var]
        error("$var not horizontally uniform.")
    end
end


"""
    reduce_nodal_stack(
        op::Function,
        grid::DiscontinuousSpectralElementGrid{T, dim, Ns},
        Q::MPIStateArray,
        vars::NamedTuple,
        var::String;
        vrange::UnitRange = 1:size(Q, 3),
    ) where {T, dim, Ns}

Reduce `var` from `vars` within `Q` over all nodal points in the specified
`vrange` of elements with `op`. Return a tuple `(result, z)` where `result` is
the final value returned by `op` and `z` is the index within `vrange` where the
`result` was determined.
"""
function reduce_nodal_stack(
    op::Function,
    grid::DiscontinuousSpectralElementGrid{T, dim, Ns},
    Q::MPIStateArray,
    vars::Type,
    var::String;
    vrange::UnitRange = 1:size(Q, 3),
    i::Int = 1,
    j::Int = 1,
) where {T, dim, Ns}

    # XXX: Needs updating for multiple polynomial orders
    # Currently only support single polynomial order
    @assert all(Ns[1] .== Ns)
    N = Ns[1]

    Nq = N + 1
    Nqk = dimensionality(grid) == 2 ? 1 : Nq

    var_names = flattenednames(vars)
    var_ind = findfirst(s -> s == var, var_names)
    if var_ind === nothing
        return
    end

    state_data = array_device(Q) isa CPU ? Q.realdata : Array(Q.realdata)
    z = vrange.start
    result = state_data[i + Nq * (j - 1), var_ind, z]
    for ev in vrange
        for k in 1:Nqk
            ijk = i + Nq * ((j - 1) + Nq * (k - 1))
            new_result = op(result, state_data[ijk, var_ind, ev])
            if !isequal(new_result, result)
                result = new_result
                z = ev
            end
        end
    end

    return (result, z)
end

"""
    reduce_element_stack(
        op::Function,
        grid::DiscontinuousSpectralElementGrid{T, dim, Ns},
        Q::MPIStateArray,
        vars::NamedTuple,
        var::String;
        vrange::UnitRange = 1:size(Q, 3),
    ) where {T, dim, Ns}

Reduce `var` from `vars` within `Q` over all nodal points in the specified
`vrange` of elements with `op`. Return a tuple `(result, z)` where `result` is
the final value returned by `op` and `z` is the index within `vrange` where the
`result` was determined.
"""
function reduce_element_stack(
    op::Function,
    grid::DiscontinuousSpectralElementGrid{T, dim, Ns},
    Q::MPIStateArray,
    vars::Type,
    var::String;
    vrange::UnitRange = 1:size(Q, 3),
) where {T, dim, Ns}

    # XXX: Needs updating for multiple polynomial orders
    # Currently only support single polynomial order
    @assert all(Ns[1] .== Ns)
    N = Ns[1]

    Nq = N + 1
    return [
        reduce_nodal_stack(
            op,
            grid,
            Q,
            vars,
            var,
            vrange = vrange,
            i = i,
            j = j,
        ) for i in 1:Nq, j in 1:Nq
    ]
end

"""
    horizontally_average!(
        grid::DiscontinuousSpectralElementGrid{T, dim, Ns},
        Q::MPIStateArray,
        i_vars,
    ) where {T, dim, Ns}

Horizontally average variables, from variable
indexes `i_vars`, in `MPIStateArray` `Q`.

!!! note
    These are not proper horizontal averages-- the main
    purpose of this method is to ensure that there are
    no horizontal fluxes for a single stack configuration.
"""
function horizontally_average!(
    grid::DiscontinuousSpectralElementGrid{T, dim, Ns},
    Q::MPIStateArray,
    i_vars,
) where {T, dim, Ns}

    # XXX: Needs updating for multiple polynomial orders
    # Currently only support single polynomial order
    @assert all(Ns[1] .== Ns)
    N = Ns[1]

    Nq = N + 1
    ArrType = typeof(Q.data)
    state_data = array_device(Q) isa CPU ? Q.realdata : Array(Q.realdata)
    Nqk = dimensionality(grid) == 2 ? 1 : Nq
    for ev in 1:size(state_data, 3), k in 1:Nqk, i_v in i_vars
        Q_sum = 0
        for i in 1:Nq, j in 1:Nq
            Q_sum += state_data[i + Nq * ((j - 1) + Nq * (k - 1)), i_v, ev]
        end
        Q_ave = Q_sum / (Nq * Nq)
        for i in 1:Nq, j in 1:Nq
            ijk = i + Nq * ((j - 1) + Nq * (k - 1))
            state_data[ijk, i_v, ev] = Q_ave
        end
    end
    Q.realdata .= ArrType(state_data)
end

get_data(solver_config, ::Prognostic) = solver_config.Q
get_data(solver_config, ::Auxiliary) = solver_config.dg.state_auxiliary
get_data(solver_config, ::GradientFlux) = solver_config.dg.state_gradient_flux

"""
    dict_of_nodal_states(
        solver_config,
        state_types = (Prognostic(), Auxiliary());
        aux_excludes = [],
        interp = false,
        )

A dictionary of single stack prognostic and auxiliary
variables at the `i=1`,`j=1` node given
 - `solver_config` a `SolverConfiguration`
 - `aux_excludes` a vector of strings containing the
    variables to exclude from the auxiliary state.
"""
function dict_of_nodal_states(
    solver_config,
    state_types = (Prognostic(), Auxiliary());
    aux_excludes = String[],
    interp = false,
)
    FT = eltype(solver_config.Q)
    all_state_vars = []
    for st in state_types
        state_vars = get_vars_from_nodal_stack(
            solver_config.dg.grid,
            get_data(solver_config, st),
            vars_state(solver_config.dg.balance_law, st, FT),
            exclude = st isa Auxiliary ? aux_excludes : String[],
            interp = interp,
        )
        push!(all_state_vars, state_vars...)
    end
    return OrderedDict(all_state_vars...)
end

end # module
