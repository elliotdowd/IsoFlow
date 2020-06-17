import matplotlib.pyplot as plt
import numpy as np

class domain:
    name = 'wedge'
    M = 6000
    N = 5400
    wedge_start = 0.5
    length = 1.5
    height = 1
    theta = np.deg2rad(30)

from gen_wedge import mesh_wedge
xx, yy = mesh_wedge(domain)

from calc_cell_metrics import cellmetrics
mesh = cellmetrics(xx, yy, domain)

# mesh plotting
fig = plt.figure('Grid Generation')
ax = fig.gca(projection='3d')

ax.plot_wireframe(mesh.xx, mesh.yy, mesh.xx*0, color='green')
#ax.plot(mesh.ccx, mesh.ccy, 'b+')
ax.view_init(-90, 90)
ax.set_proj_type('ortho')

plt.xlabel('x-coordinate (m)')
plt.ylabel('y-coordinate (m)')

plt.show()