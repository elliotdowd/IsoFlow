# wedge grid generation function
def mesh_wedge(domain):

    # import numpy
    import numpy as np

    # import domain values
    M = domain.M
    N = domain.N
    length = domain.length
    height = domain.height
    theta = domain.theta
    wedge_start = domain.obj_start

    x = np.linspace(-(1/M)*length, length*(1+(1/M)), M+3)
    y = np.linspace(-(1/N)*height, height*(1+(1/N)), N+3)

    xx, yy = np.meshgrid(x, y)
    xx = np.transpose(xx)
    yy = np.transpose(yy)

    import meshing

    meshing.mod2wedge(xx, yy, height, theta, wedge_start, M, N)

    yy = np.fliplr(yy)
    yy[:,-1] = height*(1+(1/(N)))

    return xx, yy


def mesh_airfoil(domain): 

    # import numpy
    import numpy as np
    import meshing

    # import domain values
    M = int(domain.M/2)
    N = domain.N
    length = domain.length
    height = domain.height
    theta1 = domain.theta
    af_start = domain.obj_start
    af_end = domain.obj_end

    # mesh left side of domain

    x1 = np.linspace(-(1/M)*length, length*(1+(1/M))/2, M+3)
    y1 = np.linspace(-(1/N)*height, height*(1+(1/N)), N+3)

    xx1, yy1 = np.meshgrid(x1,y1)
    xx1 = np.transpose(xx1)
    yy1 = np.transpose(yy1)

    meshing.mod2wedge(xx1, yy1, height, theta1, af_start, M, N)

    yy1 = np.fliplr(yy1)
    yy1[:,-1] = height*(1+(1/(N)))

    # determine airfoil height
    if np.sign(theta1) == 1:
        h = np.max(yy1[:,1])
    else:
        h = np.min(yy1[:,1])

    # mesh right side of domain

    x2 = np.linspace(0, length/2, M+3) + np.max(xx1)
    y2 = np.linspace(-(1/N)*height, height*(1+(1/N)), N+3)

    xx2, yy2 = np.meshgrid(x2,y2)
    xx2 = np.transpose(xx2)
    yy2 = np.transpose(yy2)

    # determine second half angle
    theta2 = np.arctan(h/((np.max(xx2)-np.min(xx2))-(length-af_end)))

    meshing.mod2wedge(xx2, yy2, height, theta2, (length-af_end)+np.max(xx1), M, N)
    yy2 = np.flipud(yy2)
    yy2 = np.fliplr(yy2)
    yy2[:,-1] = height*(1+(1/(N)))

    xx = np.vstack((xx1[1:-1,:], xx2[0:-1,:]))
    yy = np.vstack((yy1[1:-1,:], yy2[0:-1,:]))

    domain.M = M*2

    return xx, yy

