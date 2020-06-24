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
    name = 'wedge'
    M = 30
    N = 26
    obj_start = .25
    obj_end = 1.5
    length = 1.5
    height = 1.2
    theta = np.deg2rad(20)

# calculate wedge grid coordinates
from gen_grid import mesh_wedge
xx, yy = mesh_wedge(domain)

z = np.linspace(0, 2, 30)

triang = tri.Triangulation(xx.flatten(), yy.flatten())

fig1, ax1 = plt.subplots()
ax1.set_aspect('equal')
ax1.triplot(triang, 'go-', lw=1)
ax1.set_title('triplot of Delaunay triangulation')
plt.show()