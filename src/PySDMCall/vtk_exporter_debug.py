import numpy as np
from pyevtk.hl import pointsToVTK, gridToVTK
import numbers, os


"""
Example of use:

exporter = VTKExporterDebug()

for step in range(settings.n_steps):
    simulation.core.run(1)

    exporter.export_particles(simulation.core)

"""

class VTKExporterDebug:

    def __init__(self, path='.', particles_file_number=1, particles_filename="sd_points", file_num_len=4, verbose=False):
        path = os.path.join(path, 'output') 
        
        if not os.path.isdir(path):
            os.mkdir(path)

        self.particles_file_path = os.path.join(path, particles_filename)
        self.particles_file_number = particles_file_number

        products_filename = "sd_products"
        products_file_number = 1
        self.products_file_path = os.path.join(path, products_filename)
        self.products_file_number = products_file_number

        self.num_len = file_num_len

        self.verbose = verbose

    def export_particles(self, core):
        path = self.particles_file_path + '_num' + self.add_leading_zeros(self.particles_file_number)
        if self.verbose:
            print("Exporting Particles to vtk, path: " + path)

        self.particles_file_number += 1
            
        payload = {}

        particles = core.particles
        for k in particles.attributes.keys():
            if len(particles[k].shape) != 1:
                tmp = particles[k].to_ndarray(raw=True)
                tmp_dict = {k + '[' + str(i) + ']' : tmp[i] for i in range(len(particles[k].shape))}

                payload.update(tmp_dict)
            else:
                payload[k] = particles[k].to_ndarray(raw=True)

        payload.update({k: np.array(v) for k, v in payload.items() if not (v.flags['C_CONTIGUOUS'] or v.flags['F_CONTIGUOUS'])})

        if core.mesh.dimension == 2:
            x = payload['cell origin[0]'] + payload['position in cell[0]']
            y = np.full_like(x, 0)
            z = payload['cell origin[1]'] + payload['position in cell[1]']
        else:
            raise NotImplementedError("Only 2 dimensions array is supported at the moment.")

        pointsToVTK(path, x, y, z, data = payload) 

    def add_leading_zeros(self, a):
        return ''.join(['0' for i in range(self.num_len - len(str(a)))]) + str(a)
    
    def export_products(self, core):
        products = core.products

        if len(products) != 0:
            path = self.products_file_path + '_num' + self.add_leading_zeros(self.products_file_number)
            if self.verbose:
                print("Exporting Products to vtk, path: " + path)

            self.products_file_number += 1

            payload = {}
            grid = core.mesh.grid

            if core.mesh.dimension == 2:  
                keys = products.keys()

                if 'Particles Wet Size Spectrum' in keys:
                    data_shape = products['Particles Wet Size Spectrum'].shape
                elif 'Particles Dry Size Spectrum' in keys:
                    data_shape = products['Particles Dry Size Spectrum'].shape
                else:
                    data_shape = (grid[0], grid[1], 1) #

                assert(data_shape[0] == grid[0] and data_shape[1] == grid[1])    

                for k in products.keys():
                    v = products[k].get()

                    if isinstance(v, np.ndarray):
                        if v.shape == grid:                
                            tmp = np.full((1, data_shape[0], data_shape[1]), v)                            
                            payload[k] = np.ascontiguousarray(np.moveaxis(tmp, 0, 2))
                            print(v.shape)

                            for i in range(payload[k].shape[2]):
                                assert(np.array_equal(payload[k][:,:,i], v, equal_nan=True))

                        elif v.shape == data_shape:
                            payload[k] = v
                        else:
                            for i in range(v.shape[2]):
                                z = '~'+k+'['+str(i)+']'
                                payload[z] = np.array(v[:,:,i])
                                payload[z] = payload[z][:,:,np.newaxis]
                                print(z)
                                print(payload[z].shape)
                            print(f'{k} shape {v.shape} not equals data shape {data_shape}')
                    elif isinstance(v, numbers.Number):
                        payload[k] = np.full(data_shape, v)
                    else:
                        print(f'{k} export is not possible')    


                if data_shape[2] != 1:
                    lz = grid[0] / 5
                    dz = lz / data_shape[2]
                    x, y, z = np.mgrid[:grid[0] + 1, :lz+dz:dz, :grid[1] + 1]
                else:
                    x, y, z = np.mgrid[:grid[0] + 1, :1, :grid[1] + 1]
            else:
                raise NotImplementedError("Only 2 dimensions data is supported at the moment.")    


            gridToVTK(path, x, y, z, cellData = payload)

            t_passed = core.dt * core.n_steps
            
        else:
            print('No products to export')      