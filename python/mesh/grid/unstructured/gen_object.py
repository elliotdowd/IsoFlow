import numpy as np

def gen_naca4points( airfoil ):

    # import domain values
    M = airfoil.M
    naca = airfoil.naca
    a = airfoil.alpha
    L = airfoil.L

    # concentrate points on object
    x = cospace( np.linspace(0, L, M) )

    m = float( naca[0] ) / 100
    p = float( naca[1] ) / 10
    thick = float( naca[2:4] )

    if naca[0:2] == '00':

        yt = NACA4symm( x, L, thick/100 )
        x = np.hstack( [x, x] ) - float(L/2)
        y = np.hstack( [-yt, yt] )

    else:

        xL, xU, yL, yU, yc = NACA4( x, L, np.array( (m, p, thick) ) )

        x = np.hstack( [xL, np.flipud(xU)] ) - float(L/2)
        y = np.hstack( [yL, np.flipud(yU)] )

    xy = rotate( np.transpose( np.vstack([x, y]) ), a )

    airfoil.vertices = xy

    return airfoil


def gen_biconvexpoints( airfoil ):

    # import domain values
    M = airfoil.M
    t = airfoil.naca
    a = airfoil.alpha
    L = airfoil.L

    # concentrate points on object
    # x = cospace( np.linspace(0, L, M) ) / airfoil.L
    x = np.linspace(0, L, M) / airfoil.L

    yt = 2*t*(x)*( 1 - (x) )

    x = np.hstack( [x, np.flipud(x)] )*airfoil.L - float(L/2)
    y = np.hstack( [-yt, np.flipud(yt)] )*airfoil.L

    xy = rotate( np.transpose( np.vstack([x, y]) ), a )

    airfoil.vertices = xy

    return airfoil


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



# cosine spacing functions
def cospace( x ):
    import numpy as np
    min = np.min(x)
    max = np.max(x-min)
    x = x - min
    x = np.pi * x / max
    x = np.flipud( np.cos(x) + 1 ) / 2
    x = x * max + min
    return x


def front_cospace( x ):
    import numpy as np
    min = np.min(x)
    max = np.max(x-min)
    x = x - min
    x = (np.pi/2) * x / max
    x = np.flipud( np.cos(x) )
    x = x * max + min
    return x


def back_cospace( x ):
    import numpy as np
    min = np.min(x)
    max = np.max(x-min)
    x = x - min
    x = (np.pi/2) * x / max
    x = -np.cos(x) + 1
    x = x * max + min
    return x


# rotation matrix
def rotate( xy, theta):

    M, void = np.shape( xy )
    rot = np.array( [np.cos(theta), np.sin(theta), -np.sin(theta), np.cos(theta)] ).reshape([2, 2])

    for i in range(0, M):
        xy[i,:] = np.matmul( rot, np.transpose(xy[i,:]))

    return xy


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