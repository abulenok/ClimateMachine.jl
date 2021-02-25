"""
    SoilWaterParameterizations

van Genuchten, Brooks and Corey, and Haverkamp parameters for and formulation of
  - hydraulic conductivity
  - matric potential

Hydraulic conductivity can be chosen to be dependent or independent of 
impedance, viscosity and moisture.

Functions for hydraulic head, effective saturation, pressure head, matric 
potential, and the relationship between augmented liquid fraction and liquid
fraction are also included.
"""
module SoilWaterParameterizations

using DocStringExtensions
using UnPack
using ...VariableTemplates

export AbstractImpedanceFactor,
    NoImpedance,
    IceImpedance,
    impedance_factor,
    AbstractViscosityFactor,
    ConstantViscosity,
    TemperatureDependentViscosity,
    viscosity_factor,
    AbstractMoistureFactor,
    MoistureDependent,
    MoistureIndependent,
    moisture_factor,
    AbstractHydraulicsModel,
    vanGenuchten,
    BrooksCorey,
    Haverkamp,
    hydraulic_conductivity,
    effective_saturation,
    pressure_head,
    hydraulic_head,
    matric_potential,
    volumetric_liquid_fraction,
    inverse_matric_potential,
    saturated_pressure_head,
    vg_matric_potential,
    bc_matric_potential,
    vg_moisture_factor,
    bc_moisture_factor,
    haverkamp_moisture_factor

"""
    AbstractImpedanceFactor{FT <: AbstractFloat}

"""
abstract type AbstractImpedanceFactor{FT <: AbstractFloat} end

"""
    AbstractViscosityFactor{FT <: AbstractFloat}
"""
abstract type AbstractViscosityFactor{FT <: AbstractFloat} end

"""
    AbstractMoistureFactor{FT <:AbstractFloat}
"""
abstract type AbstractMoistureFactor{FT <: AbstractFloat} end


"""
    AbstractsHydraulicsModel{FT <: AbstractFloat}

Hydraulics model is used in the moisture factor in hydraulic 
conductivity and in the matric potential. The single hydraulics model 
choice sets both of these.
"""
abstract type AbstractHydraulicsModel{FT <: AbstractFloat} end


"""
    vanGenuchten{FT} <: AbstractHydraulicsModel{FT}

The necessary parameters for the van Genuchten hydraulic model; 
defaults are for Yolo light clay.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct vanGenuchten{FT,F1,F2,F3} <: AbstractHydraulicsModel{FT}
    "Exponent parameter - used in matric potential"
    n::F1
    "used in matric potential. The inverse of this carries units in 
     the expression for matric potential (specify in inverse meters)."
    α::F2
    "Exponent parameter - determined by n, used in hydraulic conductivity"
    m::F3
end

function vanGenuchten(::Type{FT};
                      n::Union{FT,Function} = FT(1.43),
                      α::Union{FT,Function} = FT(2.6)
                      ) where {FT}
    
    if typeof(n) <: AbstractFloat
        fn = (aux) -> FT(n)
    else
        fn = n
    end
    
    if typeof(α) <: AbstractFloat
        fα = (aux) -> FT(α)
    else
        fα = α
    end

    fm = (aux) -> (FT(1)-FT(1)/fn(aux))
    args = (fn, fα, fm)
    return vanGenuchten{FT, typeof.(args)...}(args...)
end


"""
    BrooksCorey{FT} <: AbstractHydraulicsModel{FT}

The necessary parameters for the Brooks and Corey hydraulic model.

Defaults are chosen to somewhat mirror the Havercamp/vG Yolo light 
clay hydraulic conductivity/matric potential.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct BrooksCorey{FT,F1,F2} <: AbstractHydraulicsModel{FT}
    "ψ_b - used in matric potential. Units of meters."
    ψb::F1
    "Exponent used in matric potential and hydraulic conductivity."
    m::F2
end

function BrooksCorey(::Type{FT};
                     ψb::Union{FT,Function} = FT(0.1656),
                     m::Union{FT,Function} = FT(0.5)
                     ) where {FT}
    if typeof(m) <: AbstractFloat
        fm = (aux) -> FT(m)
    else
        fm = m
    end
    
    if typeof(ψb) <: AbstractFloat
        fψ = (aux) -> FT(ψb)
    else
        fψ = ψb
    end
    args = (fψ, fm)
    return BrooksCorey{FT, typeof.(args)...}(args...)
end
    
"""
    Haverkamp{FT} <: AbstractHydraulicsModel{FT}

The necessary parameters for the Haverkamp hydraulic model for Yolo light
 clay.

Note that this only is used in creating a hydraulic conductivity function,
 and another formulation for matric potential must be used.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Haverkamp{FT,F1,F2,F3,F4,F5} <: AbstractHydraulicsModel{FT}
    "exponent in conductivity"
    k::F1
    "constant A (units of cm^k) using in conductivity. Our sim is in meters"
    A::F2
    "Exponent parameter - using in matric potential"
    n::F3
    "used in matric potential. The inverse of this carries units in the 
     expression for matric potential (specify in inverse meters)."
    α::F4
    "Exponent parameter - determined by n, used in hydraulic conductivity"
    m::F5
end

function Haverkamp(::Type{FT};
                   k::Union{FT,Function} = FT(1.77),
                   A::Union{FT,Function} = FT(124.6 / 100.0^1.77),
                   n::Union{FT,Function} = FT(1.43),
                   α::Union{FT,Function} = FT(2.6),
                   ) where {FT}
    if typeof(n) <: AbstractFloat
        fn = (aux) -> FT(n)
    else
        fn = n
    end
    if typeof(k) <: AbstractFloat
        fk = (aux) -> FT(k)
    else
        fk = k
    end
    if typeof(A) <: AbstractFloat
        fA = (aux) -> FT(A)
    else
        fA = A
    end
    if typeof(α) <: AbstractFloat
        fα = (aux) -> FT(α)
    else
        fα = α
    end

    fm = (aux) -> (FT(1)-FT(1)/fn(aux))
    args = (fk, fA, fn, fα, fm)
    return Haverkamp{FT, typeof.(args)...}(args...)
end

"""
    MoistureIndependent{FT} <: AbstractMoistureFactor{FT} end

Moisture independent moisture factor.
"""
struct MoistureIndependent{FT} <: AbstractMoistureFactor{FT} end


"""
    MoistureDependent{FT} <: AbstractMoistureFactor{FT} end

Moisture dependent moisture factor.
"""
struct MoistureDependent{FT} <: AbstractMoistureFactor{FT} end


function vg_moisture_factor(S_l::FT, m::FT) where {FT}
    if S_l < FT(1)
        K = sqrt(S_l) * (FT(1) - (FT(1) - S_l^(FT(1) / m))^m)^FT(2)
    else
        K = FT(1)
    end
    return K
end


"""
    moisture_factor(
        mm::MoistureDependent{FT},
        hm::vanGenuchten{FT},
        S_l::FT,
    ) where {FT}

Returns the moisture factor of the hydraulic conductivy assuming a 
MoistureDependent and van Genuchten hydraulic model.
"""
function moisture_factor(
    mm::MoistureDependent{FT},
    hm::vanGenuchten{FT},
    S_l::FT,
    aux::Vars,
) where {FT}
    m = hm.m(aux)
    return vg_moisture_factor(S_l, m)
end

function bc_moisture_factor(S_l::FT, m::FT) where {FT}
    if S_l < FT(1)
        K = S_l^(FT(2) * m + FT(3))
    else
        K = FT(1)
    end
    return K
end

"""
    moisture_factor(
        mm::MoistureDependent{FT},
        hm::BrooksCorey{FT},
        S_l::FT,
    ) where {FT}

Returns the moisture factor of the hydraulic conductivy assuming a 
MoistureDependent and Brooks/Corey hydraulic model.
"""
function moisture_factor(
    mm::MoistureDependent{FT},
    hm::BrooksCorey{FT},
    S_l::FT,
    aux::Vars
) where {FT}
    m = hm.m(aux)
    return bc_moisture_factor(S_l, m)
end


function haverkamp_moisture_factor(
    S_l::FT,
    A::FT,
    k::FT,
    α::FT,
    n::FT,
    m::FT) where {FT}
    if S_l < FT(1)
        ψ = vg_matric_potential(S_l, α, m, n)
        K = A / (A + abs(ψ)^k)
    else
        K = FT(1)
    end
    return K
end

"""
    moisture_factor(
        mm::MoistureDependent{FT},
        hm::Haverkamp{FT},
        S_l::FT,
    ) where {FT}

Returns the moisture factor of the hydraulic conductivy assuming a 
MoistureDependent and Haverkamp hydraulic model.
"""
function moisture_factor(
    mm::MoistureDependent{FT},
    hm::Haverkamp{FT},
    S_l::FT,
    aux::Vars,
) where {FT}
    k = hm.k(aux)
    A = hm.A(aux)
    n = hm.n(aux)
    m = hm.m(aux)
    α = hm.m(aux)
    return haverkamp_moisture_factor(S_l, A, k, α,n,m)
end


"""
    moisture_factor(mm::MoistureIndependent{FT},
                    hm::AbstractHydraulicsModel{FT},
                    S_l::FT,
    ) where {FT}
Returns the moisture factor in hydraulic conductivity when a 
MoistureIndependent model is chosen. Returns 1.

Note that the hydraulics model and S_l are not used, but are included 
as arguments to unify the function call.
"""
function moisture_factor(
    mm::MoistureIndependent{FT},
    hm::AbstractHydraulicsModel{FT},
    S_l::FT,
    aux::Vars
) where {FT}
    Factor = FT(1.0)
    return Factor
end

"""
    ConstantViscosity{FT} <: AbstractViscosityFactor{FT}

A model to indicate a constant viscosity - independent of temperature - 
factor in hydraulic conductivity.
"""
struct ConstantViscosity{FT} <: AbstractViscosityFactor{FT} end


"""
    TemperatureDependentViscosity{FT} <: AbstractViscosityFactor{FT}

The necessary parameters for the temperature dependent portion of hydraulic 
conductivity.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct TemperatureDependentViscosity{FT} <:
                   AbstractViscosityFactor{FT}
    "Empirical coefficient"
    γ::FT = FT(2.64e-2)
    "Reference temperature"
    T_ref::FT = FT(288.0)
end


"""
    viscosity_factor(
        vm::ConstantViscosity{FT},
        T::FT,
    ) where {FT}

Returns the viscosity factor when we choose no temperature dependence, i.e. 
a constant viscosity. Returns 1.

T is included as an argument to unify the function call.
"""
function viscosity_factor(vm::ConstantViscosity{FT}, T::FT) where {FT}
    Theta = FT(1.0)
    return Theta
end

"""
    viscosity_factor(
        vm::TemperatureDependentViscosity{FT},
        T::FT,
    ) where {FT}

Returns the viscosity factor when we choose a TemperatureDependentViscosity.
"""
function viscosity_factor(
    vm::TemperatureDependentViscosity{FT},
    T::FT,
) where {FT}
    γ = vm.γ
    T_ref = vm.T_ref
    factor = FT(γ * (T - T_ref))
    Theta = FT(exp(factor))
    return Theta
end


"""
    NoImpedance{FT} <: AbstractImpedanceFactor{FT}

A model to indicate to dependence on ice for the hydraulic conductivity.
"""
struct NoImpedance{FT} <: AbstractImpedanceFactor{FT} end



"""
    IceImpedance{FT} <: AbstractImpedanceFactor{FT}

The necessary parameters for the empirical impedance factor due to ice.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct IceImpedance{FT} <: AbstractImpedanceFactor{FT}
    "Empirical coefficient from Hansson 2014. "
    Ω::FT = FT(7)
end

"""
    impedance_factor(
        imp::NoImpedance{FT},
        θ_i::FT,
        θ_l::FT,
    ) where {FT}

Returns the impedance factor when no effect due to ice is desired. 
Returns 1.

The other arguments are included to unify the function call.
"""
function impedance_factor(imp::NoImpedance{FT}, θ_i::FT, θ_l::FT) where {FT}
    gamma = FT(1.0)
    return gamma
end

"""
    impedance_factor(
        imp::IceImpedance{FT},
        θ_i::FT,
        θ_l::FT,
    ) where {FT}

Returns the impedance factor when an effect due to the fraction of 
ice is desired. 
"""
function impedance_factor(imp::IceImpedance{FT}, θ_i::FT, θ_l::FT) where {FT}
    Ω = imp.Ω
    f_ice = θ_i / (θ_i + θ_l)
    gamma = FT(10.0^(-Ω * f_ice))
    return gamma
end

"""
    hydraulic_conductivity(
        impedance::AbstractImpedanceFactor{FT},
        viscosity::AbstractViscosityFactor{FT},
        moisture::AbstractMoistureFactor{FT},
        hydraulics::AbstractHydraulicsModel{FT},
        θ_i::FT,
        porosity::FT,
        T::FT,
        S_l::FT,
    ) where {FT}

Returns the hydraulic conductivity.
"""
function hydraulic_conductivity(
    impedance::AbstractImpedanceFactor{FT},
    viscosity::AbstractViscosityFactor{FT},
    moisture::AbstractMoistureFactor{FT},
    hydraulics::AbstractHydraulicsModel{FT},
    θ_i::FT,
    porosity::FT,
    T::FT,
    S_l::FT,
    aux::Vars,
) where {FT}
    K = FT(
        viscosity_factor(viscosity, T) *
        impedance_factor(impedance, θ_i, porosity * S_l) *
        moisture_factor(moisture, hydraulics, S_l, aux),
    )
    return K
end

"""
    hydraulic_head(z,ψ)

Return the hydraulic head.

The hydraulic head is defined as the sum of vertical height z and 
pressure head ψ; meters.
"""
hydraulic_head(z, ψ) = z + ψ

"""
    volumetric_liquid_fraction(
        ϑ_l::FT,
        eff_porosity::FT,
    ) where {FT}

Compute the volumetric liquid fraction from the effective porosity and the augmented liquid
fraction.
"""
function volumetric_liquid_fraction(ϑ_l::FT, eff_porosity::FT) where {FT}
    if ϑ_l < eff_porosity
        θ_l = ϑ_l
    else
        θ_l = eff_porosity
    end
    return θ_l
end


"""
    effective_saturation(
        porosity::FT,
        ϑ_l::FT
    ) where {FT}

Compute the effective saturation of soil.

`ϑ_l` is defined to be zero or positive. If `ϑ_l` is negative, 
hydraulic functions that take it as an argument will return 
imaginary numbers, resulting in domain errors. Exit in this 
case with an error.
"""
function effective_saturation(porosity::FT, ϑ_l::FT) where {FT}
    ϑ_l < 0 && error("Effective saturation is negative")
    S_l = ϑ_l / porosity
    return S_l
end


function saturated_pressure_head(
    ϑ_l::FT,
    porosity::FT,
    S_s::FT
) where {FT}
    return (ϑ_l - eff_porosity) / S_s
end

"""
    pressure_head(
        model::AbstractHydraulicsModel{FT},
        porosity::FT,
        S_s::FT,
        ϑ_l::FT,
        θ_i::FT,
    ) where {FT}

Determine the pressure head in both saturated and unsaturated soil. 

If ice is present, it reduces the volume available for liquid water. 
The augmented liquid fraction changes behavior depending on if this 
volume is full of liquid water vs not. Therefore, the region of saturated
vs unsaturated soil depends on porosity - θ_i, not just on porosity.  
If the liquid water is unsaturated, the usual matric potential expression
is treated as unaffected by the presence of ice.
"""
function pressure_head(
    model::AbstractHydraulicsModel{FT},
    porosity::FT,
    S_s::FT,
    ϑ_l::FT,
    θ_i::FT,
    aux::Vars,
) where {FT}
    eff_porosity = porosity - θ_i
    S_l_eff = effective_saturation(eff_porosity, ϑ_l)
    if S_l_eff < 1
        S_l = effective_saturation(porosity, ϑ_l)
        ψ = matric_potential(model, S_l, aux)
    else
        ψ = saturated_pressure_head(ϑ_l, eff_porosity, S_s)
    end
    return ψ
end



function vg_matric_potential(
    S_l::FT,
    α::FT,
    m::FT,
    n::FT,
) where {FT}
    ψ_m = -((S_l^(-FT(1) / m) - FT(1)) * α^(-n))^(FT(1) / n)
    return ψ_m
end

function bc_matric_potential(
    S_l::FT,
    ψb::FT,
    m::FT,
) where {FT}
    ψ_m = -ψb * S_l^(-FT(1) / m)
    return ψ_m
end

"""
    matric_potential(
            model::vanGenuchten{FT},
            S_l::FT
    ) where {FT}

Wrapper function which computes the van Genuchten function for matric potential.
"""
function matric_potential(model::vanGenuchten{FT}, S_l::FT, aux::Vars) where {FT}
    @unpack n, m, α = model
    m = m(aux)
    n = n(aux)
    α = α(aux)
    ψ_m = vg_matric_potential(S_l, α, m, n)
    return ψ_m
end

"""
    matric_potential(
            model::Haverkamp{FT},
            S_l::FT
    ) where {FT}

Compute the van Genuchten function as a proxy for the Haverkamp model 
matric potential (for testing purposes).
"""
function matric_potential(model::Haverkamp{FT}, S_l::FT, aux::Vars) where {FT}
    @unpack n, m, α = model
    m = m(aux)
    n = n(aux)
    α = α(aux)
    ψ_m = vg_matric_potential(S_l, α, m, n)
    return ψ_m
end

"""
    matric_potential(
            model::BrooksCorey{FT},
            S_l::FT
    ) where {FT}

Compute the Brooks and Corey function for matric potential.
"""
function matric_potential(model::BrooksCorey{FT}, S_l::FT, aux::Vars) where {FT}
    @unpack ψb, m = model
    ψb = ψb(aux)
    m = m(aux)
    ψ_m = bc_matric_potential(S_l, ψb, m)
    return ψ_m
end



"""
    inverse_matric_potential(
        model::vanGenuchten{FT},
        ψ::FT
    ) where {FT}

Compute the effective saturation given the matric potential, using
the van Genuchten formulation.
"""
function inverse_matric_potential(model::vanGenuchten{FT}, ψ::FT, aux::Vars) where {FT}

    ψ > 0 && error("Matric potential is positive")

    @unpack n, m, α = model
    n = n(aux)
    m = m(aux)
    α = α(aux)
    S = (FT(1) + (α * abs(ψ))^n)^(-m)
    return S
end


"""
    inverse_matric_potential(
        model::Haverkamp{FT}
        ψ::FT
    ) where {FT}

Compute the effective saturation given the matric potential using the 
Haverkamp hydraulics model. This model uses the van Genuchten 
formulation for matric potential.
"""
function inverse_matric_potential(model::Haverkamp{FT}, ψ::FT, aux::Vars) where {FT}
    ψ > 0 && error("Matric potential is positive")
    @unpack n, m, α = model
    n = n(aux)
    m = m(aux)
    α = α(aux)
    S = (FT(1) + (α * abs(ψ))^n)^(-m)
    return S
end


"""
    inverse_matric_potential(
        model::BrooksCorey{FT}
        ψ::FT
    ) where {FT}

Compute the effective saturation given the matric potential using the 
Brooks and Corey formulation.
"""
function inverse_matric_potential(model::BrooksCorey{FT}, ψ::FT, aux::Vars) where {FT}
    ψ > 0 && error("Matric potential is positive")

    @unpack ψb, m = model
    S = (-ψ / ψb(aux))^(-m(aux))
    return S
end

end #Module
