#####
##### Prognostic Variables
#####

prognostic_vars(water::PrescribedWaterModel) = ()
prognostic_vars(water::SoilWaterModel) = (ϑLiquid(), θIce())
prognostic_vars(heat::PrescribedTemperatureModel) = ()
prognostic_vars(heat::SoilHeatModel) = (InternalEnergy(),)

prognostic_vars(land::LandModel) =
    (prognostic_vars(land.soil.water)..., prognostic_vars(land.soil.heat)...)

get_prog_state(state, ::ϑLiquid) = (state.soil.water, :ϑ_l)
get_prog_state(state, ::θIce) = (state.soil.water, :θ_i)
get_prog_state(state, ::InternalEnergy) = (state.soil.heat, :ρe_int)
