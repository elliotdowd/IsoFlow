import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np

# n_angles = 36
# n_radii = 8
# min_radius = 0.25
# radii = np.linspace(min_radius, 0.95, n_radii)

# angles = np.linspace(0, 2 * np.pi, n_angles, endpoint=False)
# angles = np.repeat(angles[..., np.newaxis], n_radii, axis=1)
# angles[:, 1::2] += np.pi / n_angles

# x = (radii * np.cos(angles)).flatten()
# y = (radii * np.sin(angles)).flatten()

# # Create the Triangulation; no triangles so Delaunay triangulation created.
# triang = tri.Triangulation(x, y)

# # Mask off unwanted triangles.
# triang.set_mask(np.hypot(x[triang.triangles].mean(axis=1),
#                          y[triang.triangles].mean(axis=1))
#                 < min_radius)

# fig1, ax1 = plt.subplots()
# ax1.set_aspect('equal')
# ax1.triplot(triang, 'bo-', lw=1)
# ax1.set_title('triplot of Delaunay triangulation')
# plt.show()

class domain:
    name = 'airfoil'
    M = 30
    N = 26
    obj_start = .5
    obj_end = 1
    length = 1.5
    height = 1.2
    naca = '0012'
    alpha = np.rad2deg( 5 )

# calculate wedge grid coordinates
from gen_grid import mesh_wedge, mesh_cylinder, mesh_naca4
xx, yy, walls = mesh_naca4(domain)

triang = tri.Triangulation(xx.flatten(), yy.flatten())


N, void = np.shape(trang.triangles)


fig1, ax1 = plt.subplots()
ax1.set_aspect('equal')
ax1.triplot(triang, 'go-', lw=1)
ax1.set_title('triplot of Delaunay triangulation')
plt.show()


class mesh: 
    def __init__(self, elements, walls):

    elements = []
    for i in range( 0, N ):

        elements.append( element( triang.triangles[i,:], np.array([triang.x, triang.y]) ) )

    class element:
        def __init__(self, nodeID, edge):
            import numpy as np

            self.nodes = node
            self.edge = edge
