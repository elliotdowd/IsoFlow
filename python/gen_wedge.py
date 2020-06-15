# wedge grid generation function

def mesh_wedge(domain):

    # import matplotlib and numpy
    import matplotlib.pyplot as plt
    import numpy as np

    # import domain values
    M = domain.M
    N = domain.N
    length = domain.length
    height = domain.height
    theta = domain.theta
    wedge_start = domain.wedge_start

    x = np.arange(-(1/M)*length, length*(1+(1/M)), (1/M)*length)
    y = np.arange(-(1/N)*height, height*(1+(1/N)), (1/N)*height)

    xx, yy = np.meshgrid(x,y)
    xx = np.transpose(xx)
    yy = np.transpose(yy)

        
    for i in range(0, M+3):

        for j in range(0, N+3):

            # scale y-coordinates for wedge section
            if xx[i, j] >= wedge_start:
                # scale for wedge shape

                yy[i, j] = yy[i, j] - height * ((j-1)/(N+2)) * (xx[i, j] - wedge_start) * np.tan(theta)

            # flip geometry about x-axis

            yy[i, j] = height - yy[i, j]

    yy = np.fliplr(yy)


    # plotting
    plt.figure()
    ax = plt.axes(projection='3d')

    ax.plot_wireframe(xx, yy, xx*0, color='green')

    # Plot settings
    ax.view_init(90, 0)
    ax.set_xlabel('x-coordinate (m)')
    ax.set_ylabel('y-coordinate (m)')


    # Display the plot
    plt.show()

    return xx, yy

