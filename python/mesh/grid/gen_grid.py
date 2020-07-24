# wedge grid generation function
def mesh_wedge(domain):

    # import numpy
    import numpy as np
    import python.mesh.grid.meshing as meshing

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

    meshing.mod2wedge(xx, yy, height, theta, wedge_start, M, N)

    yy = np.fliplr(yy)
    #yy[:,-1] = height*(1+(1/(N)))

    obj_i = np.where(x>wedge_start)[0][0]

    # initialize list
    walls = []

    # set boundary class values
    walls.append( wall( 'domain', 0, obj_i-1, 0, 0, np.array( ( 0, 1 ) ) ) )
    walls.append( wall( 'object', obj_i, domain.M+2, 0, 0, np.array( ( 0, 1 ) ) ) )

    walls.append( wall( 'domain', 0, domain.M+2, domain.N+1, domain.N+1, np.array( ( 0, -1 ) ) ) )
    walls.append( wall( 'inlet', 0, 0, 0, domain.N+2, np.array( ( 1, 0 ) ) ) )
    walls.append( wall( 'domain', domain.M+1, domain.M+1, 0, domain.N+2, np.array( ( -1, 0 ) ) ) )

    return xx, yy, walls


def mesh_corner(domain):

    # import numpy
    import numpy as np
    import python.mesh.grid.meshing as meshing

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

    # initialize list
    walls = []

    # set boundary class values
    walls.append( wall( 'domain', 0, domain.M+2, 0, 0, np.array( ( 0, 1 ) ) ) )
    walls.append( wall( 'domain', 0, domain.M+2, domain.N+1, domain.N+1, np.array( ( 0, -1 ) ) ) )
    walls.append( wall( 'inlet', 0, 0, 0, domain.N+2, np.array( ( 1, 0 ) ) ) )
    walls.append( wall( 'domain', domain.M+1, domain.M+1, 0, domain.N+2, np.array( ( -1, 0 ) ) ) )

    return xx, yy, walls


def mesh_planar_nozzle(domain):

    # import numpy
    import numpy as np
    import python.mesh.grid.meshing as meshing

    # import domain values
    M = domain.M
    N = domain.N
    length = domain.length
    height = domain.height
    theta1 = domain.theta
    throat = domain.obj_start
    # Aratio = domain.obj_end

    # import domain values
    # theta1 = theta1*(3.1415927/180)
    y0 = 0.4
    Aratio = 2.9

    # Foesch nozzle parameters (see NAVAL ORDNANCE LABORATORY MEMORANDUM 10594 )

    r0 = y0 / theta1
    y1 = y0 * (np.sin(theta1)/theta1) * Aratio
    x1 = (3/2) * (y1-y0) * (1/np.tan(theta1))

    Nhalf = int(N/2)+1

    # mesh left side of domain

    x_i = np.linspace(0, (1+(1/M)), int(M/2)+1)
    x_i = x_i * np.sin(x_i)**(throat/height)
    y_i = np.linspace(-(1/2)*(1+(1/N)), (1/2)*(1+(1/N)), N+3)

    xL = x_i / np.max(x_i)

    y = y1 + (np.tan(theta1)/x1)* (xL**2) * (1-xL/(3*x1))
    yL = ( y-np.min(y) )
    yL = ( yL / np.max(yL) )

    xx1, yy1 = np.meshgrid(xL*length, y_i*height)
    xx1 = np.transpose(xx1)
    yy1 = np.transpose(yy1)

    yy1 = np.fliplr(yy1)

    for j in range(0, Nhalf+1):
        yy1[:,Nhalf-j] = yy1[:,Nhalf-j] + (1/2)*((height/throat)-1)*yL*((j)/Nhalf)
        yy1[:,Nhalf+j] = yy1[:,Nhalf+j] - (1/2)*((height/throat)-1)*yL*((j)/Nhalf)

    yy1 = (height / np.max(yy1[:,1:-1])) * yy1

    xx = np.vstack( [np.flipud(-xx1[-1,:]-(1/M)*length), np.flipud(-xx1), xx1[1:,:], xx1[-1,:]+(1/M)*length] )
    yy = np.vstack( [(yy1[-1,:]), np.flipud(yy1), yy1[1:,:], yy1[-1,:]] )

    yy = np.fliplr(yy)

    # initialize list
    walls = []

    # set boundary class values
    walls.append( wall( 'object', 0, domain.M+2, 0, 0, np.array( ( 0, 1 ) ) ) )
    walls.append( wall( 'object', 0, domain.M+2, domain.N+1, domain.N+1, np.array( ( 0, -1 ) ) ) )
    walls.append( wall( 'inlet', 0, 0, 0, domain.N+2, np.array( ( 1, 0 ) ) ) )
    walls.append( wall( 'domain', domain.M+1, domain.M+1, 0, domain.N+2, np.array( ( -1, 0 ) ) ) )

    return xx, yy, walls


def mesh_planar_nozzle_exit(domain):

    # import numpy
    import numpy as np
    import python.mesh.grid.meshing as meshing

    # import domain values
    M = domain.M
    N = domain.N
    length = domain.length
    height = domain.height
    theta1 = domain.theta
    throat = domain.obj_start
    exit = domain.obj_end

    # import domain values
    # theta1 = theta1*(3.1415927/180)
    y0 = 0.4
    Aratio = 2.9

    # Foesch nozzle parameters (see NAVAL ORDNANCE LABORATORY MEMORANDUM 10594 )

    r0 = y0 / theta1
    y1 = y0 * (np.sin(theta1)/theta1) * Aratio
    x1 = (3/2) * (y1-y0) * (1/np.tan(theta1))

    Nhalf = int(N/2)+1

    # mesh left side of domain

    x_i = np.linspace(0, (1+(1/M)), int(M/2)+1)
    x_i = x_i * np.sin(x_i)**(throat/height)
    y_i = np.linspace(-(1/2)*(1+(1/N)), (1/2)*(1+(1/N)), N+3)

    xL = x_i / np.max(x_i)

    y = y1 + (np.tan(theta1)/x1)* (xL**2) * (1-xL/(3*x1))
    yL = ( y-np.min(y) )
    yL = ( yL / np.max(yL) )

    xx1, yy1 = np.meshgrid(xL*length, y_i*height)
    xx1 = np.transpose(xx1)
    yy1 = np.transpose(yy1)

    yy1 = np.fliplr(yy1)

    for j in range(0, Nhalf+1):
        yy1[:,Nhalf-j] = yy1[:,Nhalf-j] + (1/2)*((height/throat)-1)*yL*((j)/Nhalf)
        yy1[:,Nhalf+j] = yy1[:,Nhalf+j] - (1/2)*((height/throat)-1)*yL*((j)/Nhalf)

    yy1 = (height / np.max(yy1[:,1:-1])) * yy1

    xx = np.vstack( [np.flipud(-xx1[-1,:]-(1/M)*length), np.flipud(-xx1), xx1[1:,:], xx1[-1,:]+(1/M)*length] )
    yy = np.vstack( [(yy1[-1,:]), np.flipud(yy1), yy1[1:,:], yy1[-1,:]] )

    yy = np.fliplr(yy)

    for i in range(1, int((M/2)*(exit/length)+1)):
        xx = np.vstack( [xx, xx[-1,:]+(2/M)*length] )
        yy = np.vstack( [yy, yy[-1,:]] )

    Msplit = domain.M
    domain.M = domain.M + int((M/2)*(exit/length))

    # initialize list
    walls = []

    # set boundary class values
    walls.append( wall( 'object', 0, Msplit, 0, 0, np.array( ( 0, 1 ) ) ) )
    walls.append( wall( 'object', 0, Msplit, domain.N+1, domain.N+1, np.array( ( 0, -1 ) ) ) )
    walls.append( wall( 'domain', Msplit+1, domain.M+2, 0, 0, np.array( ( 0, 1 ) ) ) )
    walls.append( wall( 'domain', Msplit+1, domain.M+2, domain.N+1, domain.N+1, np.array( ( 0, -1 ) ) ) )
    walls.append( wall( 'inlet', 0, 0, 0, domain.N+2, np.array( ( 1, 0 ) ) ) )
    walls.append( wall( 'domain', domain.M+1, domain.M+1, 0, domain.N+2, np.array( ( -1, 0 ) ) ) )

    return xx, yy, walls


def mesh_cylinder(domain):

    # import numpy
    import numpy as np

    # import domain values
    M = domain.M
    N = domain.N
    length = domain.length
    height = domain.height
    theta1 = domain.theta
    cyl_start = domain.obj_start
    cyl_end = domain.obj_end

    # start = cyl_start*(1-(1/M))
    r = np.linspace( cyl_start*(1-(1/M)), length*(1+(1/M)), M+3 )
    # ratio = r[-1] / (length*(1+(1/M))-start)
    # r = r / ratio
    # r = r + start

    th = np.linspace( -2*np.pi*(3/N/2), 2*np.pi*(1+(3/N/2)), N+3)


    rr, tt = np.meshgrid(r, th)
    xx = rr * np.cos(tt)
    yy = rr * np.sin(tt)

    xx = np.transpose(xx)
    yy = np.transpose(yy)

    # initialize list
    walls = []

    # set boundary class values
    walls.append( wall( 'object', 0, domain.M+2, 0, 0, np.array( ( 0, 1 ) ) ) )
    walls.append( wall( 'inlet', 0, domain.M+2, domain.N+1, domain.N+1, np.array( ( 0, -1 ) ) ) )
    walls.append( wall( 'symmetry', 0, 0, 0, domain.N+2, np.array( ( 1, 0 ) ) ) )
    walls.append( wall( 'symmetry', domain.M+1, domain.M+1, 0, domain.N+2, np.array( ( -1, 0 ) ) ) )

    return xx, yy, walls


def mesh_naca4(domain):

    import numpy as np
    
    # import domain values
    domain.M = domain.M + int(domain.M%2 == 1)
    M = domain.M
    domain.N = domain.N + int(domain.N%2 == 1)
    N = domain.N
    length = domain.length
    height = domain.height
    naca = domain.naca
    obj_start = domain.obj_start
    obj_end = domain.obj_end
    a = domain.alpha

    x = np.linspace(-(1/M)*length, length*(1+(1/M)), M+3)
    y = np.linspace(-height/2*(1+(1/N)), height/2*(1+(1/N)), N+3)

    domain.obj_i = np.where(x>obj_start)
    domain.obj_i = domain.obj_i[0][0]
    domain.obj_f = np.where(x>obj_end)
    domain.obj_f = domain.obj_f[0][0]

    m = float( naca[0] ) / 100
    p = float( naca[1] ) / 10
    thick = float( naca[2:4] )

    # focus points around x-axis/centerline
    nx = 1.5
    half = int(N/2)
    for j in range( 0, half ):
        y[half-j] = y[half-j] - 0.875*y[half-j]*(np.sinh((half-j)/half))**nx/np.sinh(1)**nx
        y[half+2+j] = y[half+2+j] - 0.875*y[half+2+j]*(np.sinh((half-j)/half))**nx/np.sinh(1)**nx

    # concentrate x points near leading and trailing edges
    nL = 1
    front = int(5)
    x_stored = x
    for i in range( 1, front ):
        x[domain.obj_i+i] = x_stored[domain.obj_i+i] + 0.75*(x[domain.obj_i]-x[domain.obj_i+i])*(np.sinh((front-i)/front))**nL/np.sinh(1)**nL
        x[domain.obj_i-i] = x_stored[domain.obj_i-i] + 0.75*(x[domain.obj_i]-x[domain.obj_i-i])*(np.sinh((front-i)/front))**nL/np.sinh(1)**nL

    nT = 1
    back = int(5)
    for i in range( 1, back ):
        x[domain.obj_f+i] = x_stored[domain.obj_f+i] + 0.75*(x[domain.obj_f]-x[domain.obj_f+i])*(np.sinh((back-i)/back))**nT/np.sinh(1)**nT
        x[domain.obj_f-i] = x_stored[domain.obj_f-i] + 0.75*(x[domain.obj_f]-x[domain.obj_f-i])*(np.sinh((back-i)/back))**nT/np.sinh(1)**nT


    xaf = x[domain.obj_i:domain.obj_f]
    c = np.max(xaf)-np.min(xaf)
    t = thick*c

    domain.wallL = half
    domain.wallU = half+2

    if naca[0:2] == '00':

        yt = NACA4symm( xaf-np.min(xaf), c, t/100)

        xx, yy = np.meshgrid(x, y)
        xx = np.transpose(xx)
        yy = np.transpose(yy)

        # focus points closer to airfoil
        ny = 2
        for j in range( 0, half ):
            yy[domain.obj_i:domain.obj_f, half-j] = -yt*(np.sinh((half-j)/half)**ny)/np.sinh(1)**ny + yy[domain.obj_i:domain.obj_f, half-j]
            yy[domain.obj_i:domain.obj_f, half+2+j] = yt*(np.sinh((half-j)/half)**ny)/np.sinh(1)**ny + yy[domain.obj_i:domain.obj_f, half+2+j]

    else:
        yt = NACA4symm( xaf-np.min(xaf), c, t/100 )

        xL, xU, yL, yU, yc = NACA4( xaf-np.min(xaf), c, np.array( (m, p, thick) ) )

        xx, yy = np.meshgrid(x, y)
        xx = np.transpose(xx)
        yy = np.transpose(yy)

        yy[domain.obj_i:domain.obj_f, half+1] = yc

        # focus points closer to airfoil
        for j in range( 0, half ):
            yy[domain.obj_i:domain.obj_f, half-j] = yL*(np.sinh((half-j)/half)**1)/np.sinh(1)**1 + yy[domain.obj_i:domain.obj_f, half-j]
            yy[domain.obj_i:domain.obj_f, half+2+j] = yU*(np.sinh((half-j)/half)**1)/np.sinh(1)**1 + yy[domain.obj_i:domain.obj_f, half+2+j]

        # shift x-coordinates near camber line
        xx[domain.obj_i:domain.obj_f, half] = np.min(xaf) + xL
        xx[domain.obj_i:domain.obj_f, half+2] = np.min(xaf) + xU

        # shift x-coordinates away from camber line
        for j in range( 1, int(N/4)):
            xx[domain.obj_i:domain.obj_f, half-j] = ((j/(N/4)) * xx[domain.obj_i:domain.obj_f, half-j] + (1-(j/(N/4))) * (np.min(xaf) + xL) )
            xx[domain.obj_i:domain.obj_f, half+2+j] = ((j/(N/4)) * xx[domain.obj_i:domain.obj_f, half+2+j] + (1-(j/(N/4))) * (np.min(xaf) + xU) )

    xx = np.array(xx, order='F')
    yy = np.array(yy, order='F')

    # initialize list
    walls = []

    # set class values
    walls.append( wall( 'object', domain.obj_i, domain.obj_f, domain.wallL, domain.wallL, np.array( ( 0, -1 ) ) ) )
    walls.append( wall( 'object', domain.obj_i, domain.obj_f, domain.wallU-1, domain.wallU-1, np.array( ( 0, 1 ) ) ) )
    walls.append( wall( 'domain', 0, domain.M+2, 0, 0, np.array( ( 0, 1 ) ) ) )
    walls.append( wall( 'domain', 0, domain.M+2, domain.N+1, domain.N+1, np.array( ( 0, -1 ) ) ) )
    walls.append( wall( 'inlet', 0, 0, 0, domain.N+2, np.array( ( 1, 0 ) ) ) )
    walls.append( wall( 'domain', domain.M+1, domain.M+1, 0, domain.N+2, np.array( ( -1, 0 ) ) ) )

    return xx, yy, walls


def mesh_biconvex(domain):

    import numpy as np
    
    # import domain values
    domain.M = domain.M + int(domain.M%2 == 1)
    M = domain.M
    domain.N = domain.N + int(domain.N%2 == 1)
    N = domain.N
    length = domain.length
    height = domain.height
    t = domain.thickness
    obj_start = domain.obj_start
    obj_end = domain.obj_end
    a = domain.alpha

    x = np.linspace(-(1/M)*length, length*(1+(1/M)), M+3)
    y = np.linspace(-height/2*(1+(1/N)), height/2*(1+(1/N)), N+3)

    domain.obj_i = np.where(x>obj_start)
    domain.obj_i = domain.obj_i[0][0]
    domain.obj_f = np.where(x>obj_end)
    domain.obj_f = domain.obj_f[0][0]

    # focus points around x-axis/centerline
    nx = 1
    half = int(N/2)
    for j in range( 0, half ):
        y[half-j] = y[half-j] - 0.875*y[half-j]*(np.sinh((half-j)/half))**nx/np.sinh(1)**nx
        y[half+2+j] = y[half+2+j] - 0.875*y[half+2+j]*(np.sinh((half-j)/half))**nx/np.sinh(1)**nx

    # concentrate x points near leading and trailing edges
    nL = 1.5
    front = int(5)
    x_stored = x
    for i in range( 1, front ):
        x[domain.obj_i+i] = x_stored[domain.obj_i+i] + 0.5*(x[domain.obj_i]-x[domain.obj_i+i])*(np.sinh((front-i)/front))**nL/np.sinh(1)**nL
        x[domain.obj_i-i] = x_stored[domain.obj_i-i] + 0.5*(x[domain.obj_i]-x[domain.obj_i-i])*(np.sinh((front-i)/front))**nL/np.sinh(1)**nL

    nT = 1.5
    back = int(5)
    for i in range( 1, back ):
        x[domain.obj_f+i] = x_stored[domain.obj_f+i] + 0.5*(x[domain.obj_f]-x[domain.obj_f+i])*(np.sinh((back-i)/back))**nT/np.sinh(1)**nT
        x[domain.obj_f-i] = x_stored[domain.obj_f-i] + 0.5*(x[domain.obj_f]-x[domain.obj_f-i])*(np.sinh((back-i)/back))**nT/np.sinh(1)**nT

    xaf = x[domain.obj_i:domain.obj_f] - x[domain.obj_i]
    c = x[domain.obj_f] - x[domain.obj_i]

    yt = 2*t*(xaf)*( 1 - (xaf/c) )
    yt = yt / c

    domain.wallL = half
    domain.wallU = half+2

    xx, yy = np.meshgrid(x, y)
    xx = np.transpose(xx)
    yy = np.transpose(yy)

    # focus points closer to airfoil
    for j in range( 0, half ):
        yy[domain.obj_i:domain.obj_f, half-j] = -yt*(np.sinh((half-j)/half)**1)/np.sinh(1)**1 + yy[domain.obj_i:domain.obj_f, half-j]
        yy[domain.obj_i:domain.obj_f, half+2+j] = yt*(np.sinh((half-j)/half)**1)/np.sinh(1)**1 + yy[domain.obj_i:domain.obj_f, half+2+j]

    xx = np.array(xx, order='F')
    yy = np.array(yy, order='F')

    # initialize list
    walls = []

    # set boundary class values
    walls.append( wall( 'object', domain.obj_i, domain.obj_f, domain.wallL, domain.wallL, np.array( ( 0, -1 ) ) ) )
    walls.append( wall( 'object', domain.obj_i, domain.obj_f, domain.wallU-1, domain.wallU-1, np.array( ( 0, 1 ) ) ) )
    walls.append( wall( 'domain', 0, domain.M+2, 0, 0, np.array( ( 0, 1 ) ) ) )
    walls.append( wall( 'domain', 0, domain.M+2, domain.N+1, domain.N+1, np.array( ( 0, -1 ) ) ) )
    walls.append( wall( 'inlet', 0, 0, 0, domain.N+2, np.array( ( 1, 0 ) ) ) )
    walls.append( wall( 'domain', domain.M+1, domain.M+1, 0, domain.N+2, np.array( ( -1, 0 ) ) ) )

    return xx, yy, walls


def mesh_capsule(domain):

    import numpy as np
    
    # import domain values
    domain.M = domain.M + int(domain.M%2 == 1)
    M = domain.M
    domain.N = domain.N + int(domain.N%2 == 1)
    N = domain.N
    length = domain.length
    height = domain.height
    t = domain.thickness
    obj_start = domain.obj_start
    obj_end = domain.obj_end
    a = domain.alpha

    x = np.linspace(-(1/M)*length, length*(1+(1/M)), M+3)
    y = np.linspace(-height/2*(1+(1/N)), height/2*(1+(1/N)), N+3)

    domain.obj_i = np.where(x>obj_start)
    domain.obj_i = domain.obj_i[0][0]
    domain.obj_f = np.where(x>obj_end)
    domain.obj_f = domain.obj_f[0][0]

    # focus points around x-axis/centerline
    nx = 3
    half = int(N/2)
    for j in range( 0, half ):
        y[half-j] = y[half-j] - 0.5*y[half-j]*(np.sinh((half-j)/half))**nx/np.sinh(1)**nx
        y[half+2+j] = y[half+2+j] - 0.5*y[half+2+j]*(np.sinh((half-j)/half))**nx/np.sinh(1)**nx

    # concentrate x points near leading and trailing edges
    # nL = 2
    # front = int(5)
    # x_stored = x
    # for i in range( 1, front ):
    #     x[domain.obj_i+i] = x_stored[domain.obj_i+i] + 0.5*(x[domain.obj_i]-x[domain.obj_i+i])*(np.sinh((front-i)/front))**nL/np.sinh(1)**nL
    #     x[domain.obj_i-i] = x_stored[domain.obj_i-i] + 0.5*(x[domain.obj_i]-x[domain.obj_i-i])*(np.sinh((front-i)/front))**nL/np.sinh(1)**nL

    # nT = 2
    # back = int(5)
    # for i in range( 1, back ):
    #     x[domain.obj_f+i] = x_stored[domain.obj_f+i] + 0.5*(x[domain.obj_f]-x[domain.obj_f+i])*(np.sinh((back-i)/back))**nT/np.sinh(1)**nT
    #     x[domain.obj_f-i] = x_stored[domain.obj_f-i] + 0.5*(x[domain.obj_f]-x[domain.obj_f-i])*(np.sinh((back-i)/back))**nT/np.sinh(1)**nT

    xaf = x[domain.obj_i:domain.obj_f] - x[domain.obj_i]
    c = x[domain.obj_f] - x[domain.obj_i]

    yL = 2*t*(xaf)*( 1 - (xaf/c) )
    yL = yL / c

    yU = 0.5 - np.abs( 0.5 - xaf/c )
    yU[yU>0.42] = 0.42
    yU = 6*t*yU


    domain.wallL = half
    domain.wallU = half+2

    xx, yy = np.meshgrid(x, y)
    xx = np.transpose(xx)
    yy = np.transpose(yy)

    # focus points closer to capsule
    nC = 2
    nT = 1
    for j in range( 0, half-1 ):
        yy[domain.obj_i:domain.obj_f, half-j] = -yL*(np.sinh((half-j)/half)**nC)/np.sinh(1)**nC + yy[domain.obj_i:domain.obj_f, half-j]
        # yy[domain.obj_i:domain.obj_f, half+2+j] = yU*(np.tanh((half-j)/half)**nT)/np.tanh(1)**nT + yy[domain.obj_i:domain.obj_f, half+2+j]
        yy[domain.obj_i:domain.obj_f, half+2+j] = yU*(half-j-1)/half + yy[domain.obj_i:domain.obj_f, half+2+j]

    xx = np.array(xx, order='F')
    # yy[:,-1] = height/2*(1+(1/(N)))
    yy = np.array(yy, order='F')

    # initialize list
    walls = []

    # set boundary class values
    walls.append( wall( 'object', domain.obj_i, domain.obj_f, domain.wallL, domain.wallL, np.array( ( 0, -1 ) ) ) )
    walls.append( wall( 'object', domain.obj_i, domain.obj_f, domain.wallU-1, domain.wallU-1, np.array( ( 0, 1 ) ) ) )
    walls.append( wall( 'inlet', 0, domain.M+2, 0, 0, np.array( ( 0, 1 ) ) ) )
    walls.append( wall( 'domain', 0, domain.M+2, domain.N+1, domain.N+1, np.array( ( 0, -1 ) ) ) )
    walls.append( wall( 'inlet', 0, 0, 0, domain.N+2, np.array( ( 1, 0 ) ) ) )
    walls.append( wall( 'inlet', domain.M+1, domain.M+1, 0, domain.N+2, np.array( ( -1, 0 ) ) ) )

    return xx, yy, walls


## NACA related functions
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


class wall:
    def __init__(self, region, xi, xf, yi, yf, n):
        import numpy as np

        self.region = region
        self.wall_n = n

        if xi == xf:
            self.wall_y = np.arange( yi, yf )
            self.wall_x = xi
        else:
            self.wall_x = np.arange( xi, xf )
            self.wall_y = yi