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


@testset "Mass Preserving Cutoff filter application" begin
    N = 3
    Ne = (1, 1, 1)

    @testset for FT in (Float64, Float32)
        @testset for dim in 2:3
            @testset for direction in (
                EveryDirection,
                HorizontalDirection,
                VerticalDirection,
            )
                brickrange = ntuple(
                    j -> range(FT(-1); length = Ne[j] + 1, stop = 1),
                    dim,
                )
                topl = ClimateMachine.Mesh.Topologies.BrickTopology(
                    MPI.COMM_WORLD,
                    brickrange,
                    periodicity = ntuple(j -> true, dim),
                )

                grid =
                    ClimateMachine.Mesh.Grids.DiscontinuousSpectralElementGrid(
                        topl,
                        FloatType = FT,
                        DeviceArray = ClimateMachine.array_type(),
                        polynomialorder = N,
                    )

                filter = ClimateMachine.Mesh.Filters.MassPreservingCutoffFilter(
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

                @testset for target in ((1, 3), (:q1, :q3))
                    Q = ClimateMachine.DGMethods.init_ode_state(
                        dg,
                        nothing,
                        dim,
                    )
                    ClimateMachine.Mesh.Filters.apply!(
                        Q,
                        target,
                        grid,
                        filter,
                        direction = direction(),
                    )
                    P = ClimateMachine.DGMethods.init_ode_state(
                        dg,
                        direction(),
                        dim,
                    )
                    @test Array(Q.data) ≈ Array(P.data)
                end
            end
        end
    end
end

##
    Rrange = grid1d(Ω.radius, Ω.radius + Ω.height, nelem = elements.vertical)

    topl = StackedCubedSphereTopology(
        mpicomm,
        elements.horizontal,
        Rrange,
        boundary = boundary, 
    )

    grid = DiscontinuousSpectralElementGrid(
        topl,
        FloatType = FT,
        DeviceArray = array,
        polynomialorder = (polynomialorder.vertical, polynomialorder.horizontal),
        #meshwarp = ClimateMachine.Mesh.Topologies.cubedshellwarp,
        meshwarp = cubedshellwarp,
    )