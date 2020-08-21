# import matplotlib and numpy modules
import numpy as np
from pytictoc import TicToc

from meshpy.geometry import Marker

t = TicToc()

class domain:
    M = 8
    N = 7
    L = 6
    h = 5

    profile_marker = Marker.FIRST_USER_MARKER

    face_marker_to_tag = {
            profile_marker: "noslip",
            Marker.MINUS_X: "inflow",
            Marker.PLUS_X: "outflow",
            Marker.MINUS_Y: "outflow",
            Marker.PLUS_Y: "outflow"
            }

class airfoil:
    naca = '2412'
    M = 100
    alpha = np.deg2rad(0)
    L = 2

# calculate wedge grid coordinates
t.tic()
from python.mesh.grid.unstructured.gen_object import gen_naca4points
from python.mesh.grid.unstructured.gen_unstruct_mesh import gen_nacamesh

airfoil = gen_naca4points( airfoil )
mesh = gen_nacamesh( domain, airfoil )

# plot mesh
from python.postprocessing.plotting import plot_unstruct_mesh
plot_unstruct_mesh(mesh)


# determine cell metrics for grid
from python.mesh.metrics.calc_cell_metrics import unstruct_cellmetrics, find_facepairs
mesh = unstruct_cellmetrics( mesh )
mesh = find_facepairs(mesh)

print('------------------------------------------------------------------')
t.toc('meshing time:')

# initialize state vector, simulation parameters and fluid properties
class parameters:
    M_in = 3
    p_in = 101325
    T_in = 300
    iterations = 100
    tolerance = -6
    CFL = 0.2
class gas:
    gamma = 1.4
    Cp = 1006
    R = 287

# initialize state vector, thermodynamic variables
t.tic()
from python.boundary.initialize import init_unstruct_state
mesh, state = init_unstruct_state( domain, mesh, parameters, gas )
t.toc('initialize time:')

# run AUSM scheme
t.tic()
from python.finite_volume.AUSM.unstruct_AUSMfamily import unstruct_AUSM, avg_vars

state = unstruct_AUSM( domain, mesh, state, parameters, gas )
t.toc('simulation time:')

# call plotting functions
from python.postprocessing.plotting import plot_unstruct_contour
plot_unstruct_contour( mesh, state )


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