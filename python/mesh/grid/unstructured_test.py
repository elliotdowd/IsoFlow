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
    obj_M = 12
    obj_start = .5
    obj_end = 1
    length = 1.5
    height = 1.2
    naca = '0012'
    alpha = np.rad2deg( 5 )

# calculate wedge grid coordinates
from gen_grid import mesh_wedge, mesh_cylinder, mesh_naca4
xx, yy, walls = mesh_naca4(domain)

cospace = np.flipud( ( np.cos( np.linspace(0, np.pi, domain.obj_M) ) + 1 ) / 2 )

triang = tri.Triangulation(xx.flatten(), yy.flatten())


N, void = np.shape(triang.triangles)


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


def NACA4symm( x, c, t ):
    import numpy as np
    y = 5*t*c * ( 0.2969*np.sqrt(x/c) - 0.126*(x/c) - 0.3516*(x/c)**2 + 0.2843*(x/c)**3 - 0.1015*(x/c)**4 )
    return y


def NACA4( x, c, naca ):
    import numpy as np
    m = naca[0]
    p = naca[1]
    t = naca[2] / 100

    pc = np.where(x>p*c)[0][0]

    # camber line angle
    theta = camber_angle( x, c, m, p)

    # camber line y-coordinate
    yc = np.zeros(len(x))
    yc[0:pc+1] = (m/p**2) * (2*p*(x[0:pc+1]/c) - (x[0:pc+1]/c)**2)
    yc[pc+1:] = (m/(1-p)**2) * ( (1-2*p) + 2*p*(x[pc+1:]/c) - (x[pc+1:]/c)**2 )

    # thickness profile
    yt = NACA4symm( x, c, t )

    # airfoil coordinates
    xU = x - yt*np.sin(theta)
    xL = x + yt*np.sin(theta)
    yU = yc + yt*np.cos(theta)
    yL = yc - yt*np.cos(theta)

    return xL, xU, yL, yU, yc


def camber_angle( x, c, m, p):
    import numpy as np

    dy = np.zeros(len(x))
    
    pc = np.where(x>p*c)[0][0]
    dy[0:pc] = (2*m/p**2) * (p - (x[0:pc]/c))
    dy[pc+1:] = (2*m/(1-p)**2) * (p - (x[pc+1:]/c))

    theta = np.arctan(dy)

    return theta
