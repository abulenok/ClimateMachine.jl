

module PySDMCall


using PyCall



export PySDM, pysdm_init, PySDMConf, PySDMKernels, PySDMSpectra

#TODO check PySDM version


mutable struct PySDM
    conf
    core
    rhod # <class 'numpy.ndarray'> czy ρ to to ??? - nie, ale bierzemy go
    field_values # 
end



mutable struct PySDMConf
    grid #tuple
    size #tuple
    dxdz #tuple
    simtime
    dt
    n_sd
    kappa # 1 hygroscopicity
    kernel # Geometric from PySDM
    spectrum_per_mass_of_dry_air # 
end

function PySDMConf(
    grid,  # (75, 75)
    size,  # (1500, 1500)
    dxdz, # (dx, dz)  
    simtime, # 1800
    dt,
    kappa, # 1
    kernel,
    spectrum_per_mass_of_dry_air # 
)

    n_sd = grid[1] * grid[2] * n_sd_per_gridbox

    physics_constants = pyimport("PySDM.physics.constants")

    size_si = (py"$size[0] * $physics_constants.si.metres", py"$size[1] * $physics_constants.si.metres")

    simtime_si = py"$simtime * $physics_constants.si.second"

    dt_si = py"$dt * $physics_constants.si.second"
    dxdz_si = (py"$dxdz[0] * $physics_constants.si.metres", py"$dxdz[1] * $physics_constants.si.metres")

    #krnl = pyimport("PySDM.dynamics.coalescence.kernels")

    #kernel = krnl.Geometric(collection_efficiency=1)


    PySDMConf(
        grid,
        size_si,
        dxdz_si,
        simtime_si,
        dt_si,
        n_sd,
        kappa,
        kernel,
        spectrum_per_mass_of_dry_air
    )
end    

       
        

        

function __init__()
    # adds current directory in the Python search path. 
    pushfirst!(PyVector(pyimport("sys")."path"), "")
    pushfirst!(PyVector(pyimport("sys")."path"), "./src/PySDMCall/")
end


function pysdm_init!(pysdm, varvals)
    pkg_formulae = pyimport("PySDM.physics.formulae")
    pkg_builder = pyimport("PySDM.builder")
    #pkg_env = pyimport("PySDM.environments")
    pkg_dynamics = pyimport("PySDM.dynamics")
    pkg_init = pyimport("PySDM.initialisation")
    pkg_backend = pyimport("PySDM.backends")


    formulae = pkg_formulae.Formulae() #state_variable_triplet="TODO") ???????????????????????
    builder = pkg_builder.Builder(n_sd=pysdm.conf.n_sd, backend=pkg_backend.CPU, formulae=formulae)

    
    
    println(size(varvals["ρ"]))
    pysdm.rhod = varvals["ρ"][:, 1, :] # pysdm rhod differs a little bit from CliMa rhod
    u1 = varvals["ρu[1]"][:, 1, :] ./ pysdm.rhod
    u3 = varvals["ρu[3]"][:, 1, :] ./ pysdm.rhod

    println("PySDM rhod size")
    println(size(pysdm.rhod))
    println("u1 & u3 size")
    println(size(u1))
    println(size(u3))

    courant_coef_u1 = pysdm.conf.dxdz[1] / pysdm.conf.dt
    courant_coef_u3 = pysdm.conf.dxdz[2] / pysdm.conf.dt
    u1 = u1 ./ courant_coef_u1
    u3 = u3 ./ courant_coef_u3

    arkw_u1 = [[ (u1[x, y-1] + u1[x, y]) / 2 for y in 2:size(u1)[2]] for x in 1:size(u1)[1]]
    arkw_u3 = [[ (u3[x-1, y] + u3[x, y]) / 2 for y in 1:size(u3)[2]] for x in 2:size(u3)[1]]

    println("Arakawa grid")
    println(size(arkw_u1))
    println(size(arkw_u3))

    courant_field = (arkw_u1, arkw_u3) 

    pysdm.field_values = py"{'th': 289.9911302100883, 'qv': 0.0075}" # from PySDM


    py"""
    from PySDM.environments._moist import _Moist
    from PySDM.state.mesh import Mesh
    from PySDM.state import arakawa_c
    import numpy as np
    from PySDM.initialisation.r_wet_init import r_wet_init, default_rtol
    from PySDM.initialisation.multiplicities import discretise_n

    class Kinematic2DMachine(_Moist):

        def __init__(self, dt, grid, size, rhod_of, field_values):
            print("GRID")
            print(grid)
            super().__init__(dt, Mesh(grid, size), [])
            self.rhod_of = rhod_of
            self.field_values = field_values

        def register(self, builder):
            super().register(builder)
            self.formulae = builder.core.formulae
            rhod = builder.core.Storage.from_ndarray(arakawa_c.make_rhod(self.mesh.grid, self.rhod_of).ravel())
            self._values["current"]["rhod"] = rhod
            self._tmp["rhod"] = rhod

        @property
        def dv(self):
            return self.mesh.dv

        def init_attributes(self, *,
                            spatial_discretisation,
                            spectral_discretisation,
                            kappa,
                            rtol=default_rtol
                            ):
            # TODO #418 move to one method
            super().sync()
            self.notify()

            attributes = {}
            with np.errstate(all='raise'):
                positions = spatial_discretisation.sample(self.mesh.grid, self.core.n_sd)
                attributes['cell id'], attributes['cell origin'], attributes['position in cell'] = \
                    self.mesh.cellular_attributes(positions)
                r_dry, n_per_kg = spectral_discretisation.sample(self.core.n_sd)
                r_wet = r_wet_init(r_dry, self, kappa=kappa, rtol=rtol, cell_id=attributes['cell id'])
                rhod = self['rhod'].to_ndarray()
                cell_id = attributes['cell id']
                domain_volume = np.prod(np.array(self.mesh.size))

            attributes['n'] = discretise_n(n_per_kg * rhod[cell_id] * domain_volume)
            attributes['volume'] = self.formulae.trivia.volume(radius=r_wet)
            attributes['dry volume'] = self.formulae.trivia.volume(radius=r_dry)

            return attributes

        def get_thd(self):
            return self.core.dynamics['EulerianAdvection'].solvers['th'].advectee.get()

        def get_qv(self):
            return self.core.dynamics['EulerianAdvection'].solvers['qv'].advectee.get()

        def sync(self):
            self.core.dynamics['EulerianAdvection'].solvers.wait()
            super().sync()
    """
    
    environment = py"Kinematic2DMachine(dt=$pysdm.conf.dt, grid=$pysdm.conf.grid, size=$pysdm.conf.size, rhod_of=$pysdm.rhod, field_values=$pysdm.field_values)"
    
    # try multiline with backslash
"""
    environment = py"Kinematic2DMachine(dt=$pysdm.conf.dt,
                                        grid=$pysdm.conf.grid,
                                        size=$pysdm.conf.size,
                                        rhod_of=$pysdm.rhod,
                                        field_values=$pysdm.field_values)"
"""
    println("Type of environment")
    println(typeof(environment))

    # here is problem
    builder.set_environment(environment)
    print(environment.mesh.__dict__)

    print("IS OKKKKKKKKKKKKKKKKKKK")
        
    builder.add_dynamic(pkg_dynamics.AmbientThermodynamics())  # sync in fields from CM 
    builder.add_dynamic(pkg_dynamics.Condensation(kappa=kappa))
    builder.add_dynamic(pkg_dynamics.EulerianAdvection(solver = CMStepper())) # sync out field to CM and launch CM advection
    builder.add_dynamic(pkg_dynamics.Displacement(
                                                courant_field=courant_field,
                                                scheme="FTBS",
                                                enable_sedimentation=True))
    builder.add_dynamic(pkg_dynamics.Coalescence(kernel=pysdm.conf.kernel))

    attributes = environment.init_attributes(spatial_discretisation=pkg_init.spatial_sampling.Pseudorandom(),
                                             spectral_discretisation=pkg_init.spectral_sampling.ConstantMultiplicity(
                                                    spectrum=pysdm.conf.spectrum_per_mass_of_dry_air
                                             ),
                                             kappa=pysdm.conf.kappa)

    pysdm.core = builder.build(attributes, products=[])
    return nothing
end


function pysdm_init1(varvals, dt, dx, simtime)
    
    println("===== PySDMCall =====")

    psdmc = pyimport("PySDMCall")


    
    return psdmc.test_call(varvals, dt, dx, simtime)
end


function PySDMKernels()
    pyimport("PySDM.physics.coalescence_kernels")
end

function PySDMSpectra()
    pyimport("PySDM.physics.spectra")
end


end # module PySDMCall