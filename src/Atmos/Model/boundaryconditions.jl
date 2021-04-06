using CLIMAParameters.Planet: cv_d, T_0

export InitStateBC

export AtmosBC,
    Impenetrable,
    FreeSlip,
    NoSlip,
    DragLaw,
    Insulating,
    PrescribedTemperature,
    PrescribedEnergyFlux,
    Adiabaticθ,
    BulkFormulaEnergy,
    Impermeable,
    OutflowPrecipitation,
    ImpermeableTracer,
    PrescribedMoistureFlux,
    BulkFormulaMoisture,
    PrescribedTracerFlux,
    NishizawaEnergyFlux

export average_density_sfc_int

const NF1 = NumericalFluxFirstOrder
const NF2 = NumericalFluxSecondOrder
const NF∇ = NumericalFluxGradient

"""
    AtmosBC(momentum = Impenetrable(FreeSlip())
            energy   = Insulating()
            moisture = Impermeable()
            precipitation = OutflowPrecipitation()
            tracer  = ImpermeableTracer())

The standard boundary condition for [`AtmosModel`](@ref). The default options imply a "no flux" boundary condition.
"""
Base.@kwdef struct AtmosBC{M, E, Q, P, TR, TC, T}
    momentum::M = Impenetrable(FreeSlip())
    energy::E = Insulating()
    moisture::Q = Impermeable()
    precipitation::P = OutflowPrecipitation()
    tracer::TR = ImpermeableTracer()
    turbconv::TC = NoTurbConvBC()
    tup::T = DispatchedTuple((), DefaultBC())
end

boundary_conditions(atmos::AtmosModel) = atmos.problem.boundaryconditions

# Most are temporarily default
default_bcs(pv::Mass) = (DefaultBC(),)
default_bcs(pv::Momentum) = (Impenetrable(FreeSlip()),)
default_bcs(pv::Energy) = (DefaultBC(),)
default_bcs(pv::ρθ_liq_ice) = (DefaultBC(),)
default_bcs(pv::TotalMoisture) = (DefaultBC(),)
default_bcs(pv::LiquidMoisture) = (DefaultBC(),)
default_bcs(pv::IceMoisture) = (DefaultBC(),)
default_bcs(pv::Rain) = (DefaultBC(),)
default_bcs(pv::Snow) = (DefaultBC(),)
default_bcs(pv::Tracers) = (DefaultBC(),)

function bc_precompute(
    atmos_bc::AtmosBC,
    atmos::AtmosModel,
    args,
    nf::Union{NF1, NF∇, NF2},
)
    dtup = DispatchedTupleDict(map(values(atmos_bc.tup)) do bc
        (bc, bc_precompute(bc, atmos, args, nf))
    end)
    return dtup
end

function boundary_state!(
    nf,
    bc::AtmosBC,
    atmos::AtmosModel,
    state⁺,
    aux⁺,
    n,
    state⁻,
    aux⁻,
    t,
    state_int⁻,
    aux_int⁻,
)
    args = (; aux⁺, state⁻, aux⁻, t, n, state_int⁻, aux_int⁻)
    set_boundary_values!(state⁺, atmos, nf, bc, args)
    atmos_boundary_state!(nf, bc, atmos, state⁺, args)
    # update moisture auxiliary variables (perform saturation adjustment, if necessary)
    # to make thermodynamic quantities consistent with the boundary state
    atmos_nodal_update_auxiliary_state!(
        moisture_model(atmos),
        atmos,
        state⁺,
        aux⁺,
        t,
    )
end

function boundary_state!(
    nf::Union{CentralNumericalFluxHigherOrder, CentralNumericalFluxDivergence},
    bc::AtmosBC,
    m::AtmosModel,
    x...,
)
    nothing
end

function atmos_boundary_state!(nf, bc::AtmosBC, atmos, state⁺, args)
    # atmos_momentum_boundary_state!(nf, bc.momentum, atmos, state⁺, args)
    atmos_energy_boundary_state!(nf, bc.energy, atmos, state⁺, args)
    atmos_moisture_boundary_state!(nf, bc.moisture, atmos, state⁺, args)
    atmos_precipitation_boundary_state!(
        nf,
        bc.precipitation,
        atmos,
        state⁺,
        args,
    )
    atmos_tracer_boundary_state!(nf, bc.tracer, atmos, state⁺, args)
    turbconv_boundary_state!(nf, bc.turbconv, atmos, state⁺, args)
end


function normal_boundary_flux_second_order!(
    nf,
    bc::AtmosBC,
    atmos::AtmosModel,
    fluxᵀn::Vars{S},
    n⁻,
    state⁻,
    diffusive⁻,
    hyperdiff⁻,
    aux⁻,
    state⁺,
    diffusive⁺,
    hyperdiff⁺,
    aux⁺,
    t,
    state_int⁻,
    diffusive_int⁻,
    aux_int⁻,
) where {S}

    _args = (;
        fluxᵀn,
        n⁻,
        state⁻,
        diffusive⁻,
        hyperdiff⁻,
        aux⁻,
        t,
        state_int⁻,
        diffusive_int⁻,
        aux_int⁻,
    )
    args = merge(_args, (precomputed = bc_precompute(bc, atmos, _args, nf),))
    set_boundary_fluxes!(fluxᵀn, atmos, nf, bc, args)

    atmos_normal_boundary_flux_second_order!(nf, bc, atmos, fluxᵀn, args)
end

function atmos_normal_boundary_flux_second_order!(
    nf,
    bc::AtmosBC,
    atmos::AtmosModel,
    args...,
)
    atmos_momentum_normal_boundary_flux_second_order!(
        nf,
        bc.momentum,
        atmos,
        args...,
    )
    atmos_energy_normal_boundary_flux_second_order!(
        nf,
        bc.energy,
        atmos,
        args...,
    )
    atmos_moisture_normal_boundary_flux_second_order!(
        nf,
        bc.moisture,
        atmos,
        args...,
    )
    atmos_precipitation_normal_boundary_flux_second_order!(
        nf,
        bc.precipitation,
        atmos,
        args...,
    )
    atmos_tracer_normal_boundary_flux_second_order!(
        nf,
        bc.tracer,
        atmos,
        args...,
    )
    turbconv_normal_boundary_flux_second_order!(nf, bc.turbconv, atmos, args...)
end

"""
    average_density(ρ_sfc, ρ_int)

Average density between the surface and the interior point, given
 - `ρ_sfc` density at the surface
 - `ρ_int` density at the interior point
"""
function average_density(ρ_sfc::FT, ρ_int::FT) where {FT <: Real}
    return FT(0.5) * (ρ_sfc + ρ_int)
end

include("bc_momentum.jl")
include("bc_energy.jl")
include("bc_moisture.jl")
include("bc_precipitation.jl")
include("bc_initstate.jl")
include("bc_tracer.jl")
