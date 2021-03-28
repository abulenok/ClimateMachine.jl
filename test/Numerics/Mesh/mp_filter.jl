# weightedsum(simulation.state, 1)

using Test
import ClimateMachine
using ClimateMachine.VariableTemplates: @vars, Vars
using ClimateMachine.Mesh.Grids:
    EveryDirection, HorizontalDirection, VerticalDirection
using ClimateMachine.MPIStateArrays: weightedsum
using ClimateMachine.BalanceLaws

import GaussQuadrature
using MPI
using LinearAlgebra

ClimateMachine.init()

struct FilterTestModel{N} <: ClimateMachine.BalanceLaws.BalanceLaw end
ClimateMachine.BalanceLaws.vars_state(::FilterTestModel, ::Auxiliary, FT) =
    @vars()
ClimateMachine.BalanceLaws.init_state_auxiliary!(::FilterTestModel, _...) =
    nothing

# Legendre Polynomials
l0(r) = 1
l1(r) = r
l2(r) = (3 * r^2 - 1) / 2
l3(r) = (5 * r^3 - 3r) / 2

low(x, y, z) = l0(x) * l0(y) + 4 * l1(x) * l1(y) + 5 * l1(z) + 6 * l1(z) * l1(x)

high(x, y, z) = l2(x) * l3(y) + l3(x) + l2(y) + l3(z) * l1(y)

filtered(::EveryDirection, dim, x, y, z) = high(x, y, z)
filtered(::VerticalDirection, dim, x, y, z) =
    (dim == 2) ? l2(x) * l3(y) + l2(y) : l3(z) * l1(y)
filtered(::HorizontalDirection, dim, x, y, z) =
    (dim == 2) ? l2(x) * l3(y) + l3(x) : l2(x) * l3(y) + l3(x) + l2(y)

ClimateMachine.BalanceLaws.vars_state(::FilterTestModel{4}, ::Prognostic, FT) =
    @vars(q1::FT, q2::FT, q3::FT, q4::FT)
function ClimateMachine.BalanceLaws.init_state_prognostic!(
    ::FilterTestModel{4},
    state::Vars,
    aux::Vars,
    localgeo,
    filter_direction,
    dim,
)
    (x, y, z) = localgeo.coord
    state.q1 = low(x, y, z) + high(x, y, z)
    state.q2 = low(x, y, z) + high(x, y, z)
    state.q3 = low(x, y, z) + high(x, y, z)
    state.q4 = low(x, y, z) + high(x, y, z)

    if !isnothing(filter_direction)
        state.q1 -= filtered(filter_direction, dim, x, y, z)
        state.q3 -= filtered(filter_direction, dim, x, y, z)
    end
end

function cubedshellwarp(a, b, c, R = max(abs(a), abs(b), abs(c)))

    function f(sR, ξ, η)
        X, Y = tan(π * ξ / 4), tan(π * η / 4)
        x1 = sR / sqrt(X^2 + Y^2 + 1)
        x2, x3 = X * x1, Y * x1
        x1, x2, x3
    end

    fdim = argmax(abs.((a, b, c)))
    if fdim == 1 && a < 0
        # (-R, *, *) : Face I from Ronchi, Iacono, Paolucci (1996)
        x1, x2, x3 = f(-R, b / a, c / a)
    elseif fdim == 2 && b < 0
        # ( *,-R, *) : Face II from Ronchi, Iacono, Paolucci (1996)
        x2, x1, x3 = f(-R, a / b, c / b)
    elseif fdim == 1 && a > 0
        # ( R, *, *) : Face III from Ronchi, Iacono, Paolucci (1996)
        x1, x2, x3 = f(R, b / a, c / a)
    elseif fdim == 2 && b > 0
        # ( *, R, *) : Face IV from Ronchi, Iacono, Paolucci (1996)
        x2, x1, x3 = f(R, a / b, c / b)
    elseif fdim == 3 && c > 0
        # ( *, *, R) : Face V from Ronchi, Iacono, Paolucci (1996)
        x3, x2, x1 = f(R, b / c, a / c)
    elseif fdim == 3 && c < 0
        # ( *, *,-R) : Face VI from Ronchi, Iacono, Paolucci (1996)
        x3, x2, x1 = f(-R, b / c, a / c)
    else
        error("invalid case for cubedshellwarp: $a, $b, $c")
    end

    return x1, x2, x3
end

@testset "Mass Preserving Cutoff Filter Conservation Test" begin
    N = 3
    Ne = (1, 1, 1)
    dim = 3
    # dim, direction
    @testset for FT in (Float64, )
        Rrange = [FT(1.0), FT(1.2)]

        topl = ClimateMachine.Mesh.Grids.StackedCubedSphereTopology(
            MPI.COMM_WORLD,
            1,
            Rrange,
            boundary = (5,6), 
        )

        grid = ClimateMachine.Mesh.Grids.DiscontinuousSpectralElementGrid(
            topl,
            FloatType = FT,
            DeviceArray = ClimateMachine.array_type(),
            polynomialorder = (N, N),
            meshwarp = cubedshellwarp,
        )


        mp_filter = ClimateMachine.Mesh.Filters.MassPreservingCutoffFilter(
            grid,
            2,
        )

        reg_filter = ClimateMachine.Mesh.Filters.CutoffFilter(
            grid,
            2,
        )

        model = FilterTestModel{4}()
        dg = ClimateMachine.DGMethods.DGModel(
            model,
            grid,
            nothing,
            nothing,
            nothing;
            state_gradient_flux = nothing,
        )

        # test mp filter
        filter = mp_filter
        Q = ClimateMachine.DGMethods.init_ode_state(
                dg,
                nothing,
                dim,
        )
        sum_before_1 = weightedsum(Q, 1)
        sum_before_2 = weightedsum(Q, 2)
        sum_before_3 = weightedsum(Q, 3)
        target = 1:3
        ClimateMachine.Mesh.Filters.apply!(
            Q,
            target,
            grid,
            filter,
        )
        sum_after_1 = weightedsum(Q, 1)
        sum_after_2 = weightedsum(Q, 2)
        sum_after_3 = weightedsum(Q, 3)

        @test sum_before_1 ≈ sum_after_1
        @test sum_before_2 ≈ sum_after_2
        @test sum_before_3 ≈ sum_after_3


        # test regular filter
        filter = reg_filter
        Q = ClimateMachine.DGMethods.init_ode_state(
                dg,
                nothing,
                dim,
        )
        sum_before_1 = weightedsum(Q, 1)
        sum_before_2 = weightedsum(Q, 2)
        sum_before_3 = weightedsum(Q, 3)
        target = 1:3
        ClimateMachine.Mesh.Filters.apply!(
            Q,
            target,
            grid,
            filter,
        )
        sum_after_1 = weightedsum(Q, 1)
        sum_after_2 = weightedsum(Q, 2)
        sum_after_3 = weightedsum(Q, 3)

        @test !(sum_before_1 ≈ sum_after_1)
        @test !(sum_before_2 ≈ sum_after_2)
        @test !(sum_before_3 ≈ sum_after_3)

        end
end
