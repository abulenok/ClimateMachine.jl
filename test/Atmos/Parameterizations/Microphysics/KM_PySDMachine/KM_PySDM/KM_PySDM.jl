import Thermodynamics
const THDS = Thermodynamics

using PyCall

function init!(pysdm, varvals)
    pkg_formulae = pyimport("PySDM.physics.formulae")
    pkg_builder = pyimport("PySDM.builder")
    pkg_dynamics = pyimport("PySDM.dynamics")
    pkg_init = pyimport("PySDM.initialisation")
    pkg_backend = pyimport("PySDM.backends")
    pkg_clima = pyimport("clima_hydrodynamics")
    pkg_debug_vtk_exp = pyimport("vtk_exporter_debug")

    print("pysdm.config.n_sd: ")
    println(pysdm.config.n_sd)

    formulae = pkg_formulae.Formulae()
    builder = pkg_builder.Builder(n_sd=pysdm.config.n_sd, backend=pkg_backend.CPU, formulae=formulae)


    pysdm.rhod = varvals["ρ"][:, 1, :] # PySDM.rhod != CliMa.rho
    u1 = varvals["ρu[1]"][:, 1, :] ./ pysdm.rhod
    u3 = varvals["ρu[3]"][:, 1, :] ./ pysdm.rhod

    courant_coef_u1 = pysdm.config.dxdz[1] / pysdm.config.dt
    courant_coef_u3 = pysdm.config.dxdz[2] / pysdm.config.dt
    u1 = u1 ./ courant_coef_u1
    u3 = u3 ./ courant_coef_u3

    arkw_u1 = [ (u1[y, x-1] + u1[y, x]) / 2 for y in 1:size(u1)[1], x in 2:size(u1)[2]]
    arkw_u3 = [ (u3[y-1, x] + u3[y, x]) / 2 for y in 2:size(u3)[1], x in 1:size(u3)[2]]

    courant_field = (arkw_u1, arkw_u3)

    pysdm.rhod = bilinear_interpol(pysdm.rhod)

    pkg_env = pyimport("Kinematic2DMachine")

    environment = pkg_env.Kinematic2DMachine(dt=pysdm.config.dt, grid=pysdm.config.grid, size=pysdm.config.size, clima_rhod=pysdm.rhod)

    builder.set_environment(environment)
    println(environment.mesh.__dict__)

    # builder.add_dynamic(pkg_dynamics.AmbientThermodynamics()) # override env.sync()   # sync in fields from CM  w tym miejscu pobieramy pola z CliMa

    builder.add_dynamic(pkg_dynamics.Condensation())


    pysdm_thd = varvals["theta_dry"][:, 1, :]
    pysdm_thd = bilinear_interpol(pysdm_thd)

    pysdm_qv = varvals["q_vap"][:, 1, :]
    pysdm_qv = bilinear_interpol(pysdm_qv)


    builder.add_dynamic(pkg_clima.ClimateMachine(py"{'qv': $pysdm_qv, 'thd': $pysdm_thd}"))

    # has sense only for multithreading
    # builder.add_dynamic(pkg_dynamics.EulerianAdvection(solver = CMStepper())) # sync out field to CM and launch CM advection

    displacement = pkg_dynamics.Displacement(enable_sedimentation=false)
    
    builder.add_dynamic(displacement) # enable_sedimentation=true # scheme="FTBS"

    builder.add_dynamic(pkg_dynamics.Coalescence(kernel=pysdm.config.kernel))


    attributes = environment.init_attributes(spatial_discretisation=pkg_init.spatial_sampling.Pseudorandom(),
                                             spectral_discretisation=pkg_init.spectral_sampling.ConstantMultiplicity(
                                                    spectrum=pysdm.config.spectrum_per_mass_of_dry_air
                                             ),
                                             kappa=pysdm.config.kappa)


    pkg_PySDM_products = pyimport("PySDM.products")
    pysdm_products = []
    push!(pysdm_products, pkg_PySDM_products.WaterMixingRatio(name="qc", description_prefix="liquid", radius_range=(0.0, Inf)))
    
    push!(pysdm_products, pkg_PySDM_products.RelativeHumidity())
    push!(pysdm_products, pkg_PySDM_products.ParticleMeanRadius())
    push!(pysdm_products, pkg_PySDM_products.PeakSupersaturation())
    push!(pysdm_products, pkg_PySDM_products.ActivatingRate())
    push!(pysdm_products, pkg_PySDM_products.DeactivatingRate())
    push!(pysdm_products, pkg_PySDM_products.CondensationTimestepMin())
    push!(pysdm_products, pkg_PySDM_products.CondensationTimestepMax())
    push!(pysdm_products, pkg_PySDM_products.CloudDropletEffectiveRadius(radius_range=(0.0, Inf)))

    pysdm.core = builder.build(attributes, products=pysdm_products)

    displacement.upload_courant_field(courant_field)

    ####
    pysdm.exporter = pkg_debug_vtk_exp.VTKExporterDebug(verbose=true)

    print("Products keys: ")
    println(pysdm.core.products.keys)

    print("PySDM Dynamics: ")
    println(keys(pysdm.core.dynamics))
    
    return nothing
end


function do_step!(pysdm, varvals, t)
    
    dynamics = pysdm.core.dynamics

    #TODO: add Displacement 2 times: 1 for Condensation and 1 for Advection

    dynamics["Displacement"]()

    update_pysdm_fields!(pysdm, varvals, t)

    pysdm.core.env.sync()

    dynamics["ClimateMachine"]()
    dynamics["Condensation"]()

    pysdm.core._notify_observers()

    #env.sync() # take data from CliMA
    # upd CliMa state vars
end


function update_pysdm_fields!(pysdm, vals, t)

    liquid_water_mixing_ratio = pysdm.core.products["qc"].get() * 1e-3

    # TODO: instead of liquid_water_mixing_ratio should be liquid_water_specific_humidity
    liquid_water_specific_humidity = liquid_water_mixing_ratio 

    q_tot = vals["q_tot"][:, 1, :]
    q_tot = bilinear_interpol(q_tot)

    q = THDS.PhasePartition.(q_tot, liquid_water_specific_humidity, .0)

    qv = THDS.vapor_specific_humidity.(q)

    e_int = vals["e_int"][:, 1, :]
    e_int = bilinear_interpol(e_int)

    # TODO - AJ shouldnt we compute new e_int and new T based on new pp?
    T = THDS.air_temperature.(param_set, e_int, q)
    
    ρ = pysdm.rhod
    thd = THDS.dry_pottemp.(param_set, T, ρ) # rho has to be rhod (dry)

    pysdm.core.dynamics["ClimateMachine"].set_thd(thd)
    pysdm.core.dynamics["ClimateMachine"].set_qv(qv)

    # TODO: passing effective_radius to ClimateMachine (Auxiliary)
    return nothing
end


function bilinear_interpol(A)

    A = [ (A[y, x-1] + A[y, x]) / 2 for y in 1:size(A)[1], x in 2:size(A)[2]]
    A = [ (A[y-1, x] + A[y, x]) / 2 for y in 2:size(A)[1], x in 1:size(A)[2]]
    return A
end

function PySDMKernels()
    pyimport("PySDM.physics.coalescence_kernels")
end

function PySDMSpectra()
    pyimport("PySDM.physics.spectra")
end

function export_particles_to_vtk(pysdm)
    if !isnothing(pysdm.exporter)
        pysdm.exporter.export_particles(pysdm.core)
    end
end

function export_plt(var, title, t)
    py"""
    from matplotlib.pyplot import cm
    import numpy as np
    import matplotlib.pyplot as plt

    def plot_vars(A, title=None):
        # Contour Plot
        plt.clf()
        X, Y = np.mgrid[0:A.shape[0], 0:A.shape[1]]
        Z = A
        cp = plt.contourf(X, Y, Z)
        cb = plt.colorbar(cp)
        if title:
            plt.title(title)

        plt.show()
        return plt
    """

    println(string("Exporting ", title, " plot"))
    plt = py"plot_vars($var, title=$title)"
    plt.savefig(string(title, t, ".png"))
end

#debug
function debug_export_products_to_vtk(pysdm)
    if !isnothing(pysdm.exporter)
        pysdm.exporter.export_products(pysdm.core)
    end
end
