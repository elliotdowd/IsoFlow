# import matplotlib and numpy modules
import matplotlib.pyplot as plt
import numpy as np

class domain:
    name = 'wedge'
    M = 30
    N = 26
    wedge_start = 0.25
    length = 1.5
    height = 1.3
    theta = np.deg2rad(20)

# calculate wedge grid coordinates
from gen_wedge import mesh_wedge
xx, yy = mesh_wedge(domain)

# determine cell metrics for grid
from calc_cell_metrics import cellmetrics
mesh = cellmetrics(xx, yy, domain)

# initialize state vector, simulation parameters and fluid properties
class parameters:
    M_in = 3
    p_in = 101325
    T_in = 300
    iterations = 1000
    tolerance = 1e-6
    CFL = 0.25
class gas:
    gamma = 1.4
    Cp = 1006
    R = 287

from initialize import init_state
state = init_state(domain, mesh, parameters, gas)

# run AUSM scheme
from schemes import AUSM
state = AUSM( domain, mesh, parameters, state, gas )

# mesh plotting
fig = plt.figure('Grid Generation')
ax = fig.gca(projection='3d')

ax.plot_wireframe(mesh.xx, mesh.yy, mesh.xx*0, color='green')
ax.plot(mesh.xxc, mesh.yyc, 'b+')
ax.view_init(-90, 90)
ax.set_proj_type('ortho')

plt.xlabel('x-coordinate (m)')
plt.ylabel('y-coordinate (m)')

#plt.show()

# surf plotting

fig = plt.figure('Contour Plotting')
ax = fig.gca(projection='3d')

from matplotlib import cm
ax.plot_surface(mesh.xxc, mesh.yyc, state.Mach, cmap=cm.jet, linewidth=0)
ax.view_init(-90, 90)
ax.set_proj_type('ortho')

plt.xlabel('x-coordinate (m)')
plt.ylabel('y-coordinate (m)')

plt.show()