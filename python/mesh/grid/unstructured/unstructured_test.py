import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np

import meshpy
from meshpy.triangle import MeshInfo, build
from meshpy import triangle

# class definitions
class body:
    def __init__(self, M, angle, pos, size, spec ):
        # super().__init__()

        self.M = M
        self.alpha = angle * np.pi / 180
        self.pos = pos
        self.size = size
        self.spec = spec

# calculate wedge grid coordinates
from gen_object import gen_naca4points
from gen_grid import mesh_naca4

class domain:
    M = 32
    N = 24
    length = 4
    height = 3
    naca = '0012'
    obj_start = 1
    obj_end = 3
    alpha = 5

xx, yy, walls = mesh_naca4(domain)

# object generation
pos = np.array( [2, 0] )
airfoil = body(80, 5, pos, 2, '4212')
gen_naca4points(airfoil)

# meshpy interfacing
x = xx.flatten()
y = yy.flatten()

mesh_info = MeshInfo()
points = []

for i in range(0, len(x)):
    if i == 0:
        points = ( x[i], y[i] )
    else:
        points = np.vstack( [ points, ( x[i], y[i] ) ] )

mesh_info.set_points( points )
mesh = build(mesh_info)


fig1, ax1 = plt.subplots()
ax1.set_aspect('equal')
ax1.triplot(triang, 'go-', lw=1)
# ax1.plot(airfoil.x, airfoil.y, 'b.')
# ax1.set_title('triplot of Delaunay triangulation')
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


