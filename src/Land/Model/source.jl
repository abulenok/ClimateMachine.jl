#### Land sources
export PhaseChange

function heaviside(x::FT) where {FT}
    if x >= FT(0)
        output = FT(1)
    else
        output = FT(0)
    end
    return output
end


"""
    PhaseChange{FT} <: TendencyDef{Source}

The function which computes the freeze/thaw source term for Richard's equation,
assuming the timescale is the maximum of the thermal timescale and the timestep.
"""
Base.@kwdef struct PhaseChange{FT} <: TendencyDef{Source}
    "Timestep"
    Δt::FT = FT(NaN)
    "Timescale for temperature changes"
    τLTE::FT = FT(NaN)
end

prognostic_vars(::PhaseChange) = (ϑLiquid(), θIce())

function precompute(land::LandModel, args, tt::Source)
    dtup = DispatchedSet(map(land.source) do s
        (s, precompute(s, land, args, tt))
    end)
    return (; dtup)
end

function precompute(source_type::PhaseChange, land::LandModel, args, tt::Source)
    @unpack state, diffusive, aux, t, direction = args

    FT = eltype(state)

    param_set = parameter_set(land)
    _ρliq = FT(ρ_cloud_liq(param_set))
    _ρice = FT(ρ_cloud_ice(param_set))
    _Tfreeze = FT(T_freeze(param_set))
    _LH_f0 = FT(LH_f0(param_set))
    _g = FT(grav(param_set))

    ϑ_l, θ_i = get_water_content(land.soil.water, aux, state, t)
    eff_porosity = land.soil.param_functions.porosity - θ_i
    θ_l = volumetric_liquid_fraction(ϑ_l, eff_porosity)
    T = get_temperature(land.soil.heat, aux, t)

    # The inverse matric potential is only defined for negative arguments.
    # This is calculated even when T > _Tfreeze, so to be safe,
    # take absolute value and then pick the sign to be negative.
    # But below, use heaviside function to only allow freezing in the
    # appropriate temperature range (T < _Tfreeze)
    ψ = -abs(_LH_f0 / _g / _Tfreeze * (T - _Tfreeze))
    hydraulics = land.soil.water.hydraulics
    θstar =
        land.soil.param_functions.porosity *
        inverse_matric_potential(hydraulics, ψ)

    τft = max(source_type.Δt, source_type.τLTE)
    freeze_thaw =
        FT(1) / τft * (
            _ρliq *
            (θ_l - θstar) *
            heaviside(_Tfreeze - T) *
            heaviside(θ_l - θstar) - _ρice * θ_i * heaviside(T - _Tfreeze)
        )
    return (; freeze_thaw)
end

function source(::ϑLiquid, s::PhaseChange, land::LandModel, args)
    @unpack state = args
    @unpack freeze_thaw = args.precomputed.dtup[s]
    FT = eltype(state)
    _ρliq = FT(ρ_cloud_liq(land.param_set))
    return -freeze_thaw / _ρliq
end

function source(::θIce, s::PhaseChange, land::LandModel, args)
    @unpack state = args
    @unpack freeze_thaw = args.precomputed.dtup[s]
    FT = eltype(state)
    _ρice = FT(ρ_cloud_ice(land.param_set))
    return freeze_thaw / _ρice
end
