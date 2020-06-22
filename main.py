# import matplotlib and numpy modules
import numpy as np

class domain:
    name = 'wedge'
    M = 42
    N = 36
    obj_start = 0.75
    obj_end = 1.5
    length = 1.5
    height = 1.3
    theta = np.deg2rad(20)

# calculate wedge grid coordinates
from gen_grid import mesh_wedge
xx, yy = mesh_wedge(domain)

# determine cell metrics for grid
from calc_cell_metrics import cellmetrics
mesh = cellmetrics(xx, yy, domain)

# initialize state vector, simulation parameters and fluid properties
class parameters:
    M_in = 1.4
    p_in = 101325
    T_in = 300
    iterations = 2500
    tolerance = 1e-6
    CFL = 0.2
class gas:
    gamma = 1.4
    Cp = 1006
    R = 287

from initialize import init_state
state = init_state(domain, mesh, parameters, gas)

# run AUSM scheme
from schemes import AUSM
state = AUSM( domain, mesh, parameters, state, gas )

# call plotting functions
from plotting import plot_mesh, plot_contour
#plot_mesh(mesh)
plot_contour(domain, mesh, state)