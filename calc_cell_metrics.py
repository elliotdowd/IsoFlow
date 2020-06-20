def TicTocGenerator():
    # Generator that returns time differences
    import time
    ti = 0           # initial time
    tf = time.time() # final time
    while True:
        ti = tf
        tf = time.time()
        yield tf-ti # returns the time difference

TicToc = TicTocGenerator() # create an instance of the TicTocGen generator

# This will be the main function through which we define both tic() and toc()
def toc(tempBool=True):
    # Prints the time difference yielded by generator instance TicToc
    tempTimeInterval = next(TicToc)
    if tempBool:
        print( "s_proj time: %f seconds.\n" %tempTimeInterval )

def tic():
    # Records a time in TicToc, marks the beginning of a time interval
    toc(False)


# inverse metrics, projected cell face areas, cell volumes, centroids

def cellmetrics(xx, yy, domain):

    import numpy as np
    import cells

    # format grid variable for input into fortran subroutines
    grid = np.array((xx, yy), dtype='float', order='F')
    grid = np.array(grid).reshape(domain.M+3,domain.N+3,2)

    # calculate projected cell face areas and cell face areas 
    # S_proj = [S_zeta.x, S_zeta.y, S_eta.x, S_eta.y S_zeta, S_eta]

    tic()

    s_proj = np.zeros((domain.M+2, domain.N+2, 6), dtype='float', order='F')
    s_proj[:,:,0] = yy[1:domain.M+3, 1:domain.N+3] - yy[1:domain.M+3, 0:domain.N+2]
    s_proj[:,:,1] = -( xx[1:domain.M+3, 1:domain.N+3] - xx[1:domain.M+3, 0:domain.N+2] )
    s_proj[:,:,2] = -( yy[1:domain.M+3, 1:domain.N+3] - yy[0:domain.M+2, 1:domain.N+3] )
    s_proj[:,:,3] = xx[1:domain.M+3, 1:domain.N+3] - xx[0:domain.M+2, 1:domain.N+3]
    s_proj[:,:,4] = np.sqrt( s_proj[:,:,0]**2 + s_proj[:,:,1]**2 )
    s_proj[:,:,5] = np.sqrt( s_proj[:,:,2]**2 + s_proj[:,:,3]**2 )

    toc()

    # calculate cell areas
    area = np.zeros((domain.M+2, domain.N+2), dtype='float', order='F')
    cells.calc_cellarea(4, xx, yy, area, domain.M+2, domain.N+2)

    # calculate cell centroids
    xxc = np.zeros((domain.M+2, domain.N+2), dtype='float', order='F')
    yyc = np.zeros((domain.M+2, domain.N+2), dtype='float', order='F')

    cells.calc_cellcentroids(xx, yy, xxc, yyc, area, domain.M+2, domain.N+2)

    # create mesh class
    class mesh:
        pass

    mesh.xx = xx
    mesh.yy = yy
    mesh.xxc = xxc
    mesh.yyc = yyc
    mesh.dV = area
    mesh.s_proj = s_proj

    return mesh