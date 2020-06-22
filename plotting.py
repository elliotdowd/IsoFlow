
# plot geometry mesh

def plot_mesh(mesh):

    import matplotlib.pyplot as plt

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

# contour plotting
def plot_contour(domain, mesh, state):
    import matplotlib.pyplot as plt
    from matplotlib import cm
    import numpy as np

    fig, ax = plt.subplots()
    ax.set_facecolor( (0.3, 0.3, 0.3) )

    cont = ax.contourf(mesh.xxc, mesh.yyc, state.Mach, cmap=cm.jet)
    plt.gca().set_aspect('equal', adjustable='box')

    # colorbar settings
    CB = fig.colorbar(cont, shrink=0.8, extend='both')

    # plot labeling 
    plt.xlabel('x-coordinate (m)')
    plt.ylabel('y-coordinate (m)')

    plt.show()