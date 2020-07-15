# import matplotlib and numpy modules
import numpy as np
from pytictoc import TicToc

t = TicToc()

class domain:
    name = 'wedge'
    M = 48
    N = 36
    obj_start = 1
    obj_end = 30000
    length = 4
    height = 2.4
    theta = np.deg2rad(20)

# calculate wedge grid coordinates
t.tic()
from python.mesh.grid.gen_grid import mesh_wedge
xx, yy, walls = mesh_wedge(domain)

boundary = init_boundary( walls )

# from plotting import plot_mesh
# class mesh:
#     pass
# mesh.xx = xx
# mesh.yy = yy
# plot_mesh(mesh)

# determine cell metrics for grid
from python.mesh.metrics.calc_cell_metrics import cellmetrics
mesh = cellmetrics(xx, yy, domain)

print('------------------------------------------------------------------')
t.toc('meshing time:')

# initialize state vector, simulation parameters and fluid properties
class parameters:
    M_in = 2.6
    p_in = 101325
    T_in = 300
    iterations = 2000
    tolerance = -6
    CFL = 0.25
class gas:
    gamma = 1.4
    Cp = 1006
    R = 287

# initialize state vector, thermodynamic variables
t.tic()
from python.boundary.initialize import init_state
state = init_state(domain, mesh, boundary, parameters, gas)
t.toc('initialize time:')

# run AUSM scheme
t.tic()
from python.finite_volume.AUSM import AUSM, AUSMplusup, AUSMDV
from python.finite_volume.Roe import RoeFDS, RoeFVS

state = RoeFDS( domain, mesh, boundary, parameters, state, gas )
t.toc('simulation time:')

# call plotting functions
from python.postprocessing.plotting import plot_mesh, plot_contour
plot_contour(domain, mesh, state)


def init_boundary( wall ):

    for obj in wall:
        if obj.region == 'domain':
            # domain bottom wall
            if obj.wall_n[1] == 1:
                obj.type = 'Inviscid Wall'
                obj.thermal = 'Adiabatic'
                obj.Tw = 300

            # domain top wall
            elif obj.wall_n[1] == -1:
                obj.type = 'Outflow'
                obj.thermal = 'Adiabatic'
                obj.Tw = 300

            else:
                obj.type = 'Outflow'
    
    return wall