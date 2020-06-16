import matplotlib.pyplot as plt
import numpy as np

class domain:
    name = 'wedge'
    M = 30
    N = 26
    wedge_start = 0.5
    length = 1.5
    height = 1
    theta = np.deg2rad(10)

from gen_wedge import mesh_wedge

xx, yy = mesh_wedge(domain)

import cells

# format grid variable for input into fortran subroutines
grid = np.array((xx, yy), dtype='float', order='F')
grid = np.array(grid).reshape(domain.M+3,domain.N+3,2)

# calculate x and y inverse metrics
xy_inv = np.zeros((domain.M+2, domain.N+2, 4), dtype='float', order='F')
cells.inv_metrics(grid, xy_inv, domain.M+2, domain.N+2)

# calculate projected cell face areas and cell face areas 
# S_proj = [S_zeta.x, S_zeta.y, S_eta.x, S_eta.y S_zeta, S_eta]

s_proj = np.zeros((domain.M+2, domain.N+2, 6), dtype='float', order='F')
cells.face_areas(xy_inv, s_proj, domain.M+2, domain.N+2)

# calculate cell areas
area = np.zeros((domain.M+2, domain.N+2), dtype='float', order='F')
cells.calc_cellarea(4, xx, yy, area, domain.M+2, domain.N+2)

# calculate cell centroids
ccx = np.zeros((domain.M+2, domain.N+2), dtype='float', order='F')
ccy = np.zeros((domain.M+2, domain.N+2), dtype='float', order='F')

cells.calc_cellcentroids(xx, yy, ccx, ccy, area, domain.M+2, domain.N+2)

# create mesh class
class mesh:
    xx = xx
    yy = yy
    ccx = ccx
    ccy = ccy
    dV = area
