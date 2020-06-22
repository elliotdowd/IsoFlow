
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
    #import plotly.graph_objects as go
    import numpy as np

    fig = plt.figure(figsize=(7, 9))
    gs = plt.GridSpec(7, 9)
    ax1 = fig.add_subplot(gs[0:5, :])
    ax2 = fig.add_subplot(gs[6:, :])

    #ax1 = fig.add_axes([0, 0, domain.length/max(domain.length,domain.height), domain.height/max(domain.length,domain.height)])
    ax1.set_title("Contour Plotting")
    ax1.set_xlabel('x-coordinate (m)')
    ax1.set_ylabel('y-coordinate (m)') 
    #fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6,18))
    #ax1.set_facecolor( (0.3, 0.3, 0.3) )

    # contour plotting
    cont = ax1.contourf(mesh.xxc[1:-2,1:-1], mesh.yyc[1:-2,1:-1], state.Mach[1:-2,1:-1], 500, cmap=cm.jet)
    ax1.axis('equal')
    ax1.set_xlim(mesh.xx[0,0], mesh.xx[-1,-1])
    ax1.set_ylim(np.min(mesh.yy), np.max(mesh.yy))

    #plt.gca().set_aspect('equal', adjustable='box')
    #plt.xlabel('x-coordinate (m)')
    #plt.ylabel('y-coordinate (m)')


    ax2.plot(np.arange(1, len(state.res), 1), state.res[1:])
    ax2.set_xlabel('Iterations')
    ax2.set_ylabel('Residual') 
    #plt.xlabel('Iterations')
    #plt.ylabel('Residual')

    # colorbar settings
    CB = fig.colorbar(cont, shrink=0.7, extend='both', ax=ax1)
    CB.set_label('Mach Number', rotation=90)

    # plot labeling 


    plt.show()