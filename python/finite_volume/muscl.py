import numpy as np


def MUSCL( Q, eps, kap, limiter ):

    M, N, void = Q.shape

    QL_half = np.zeros( (M-1, N, 4) )
    QR_half = np.zeros( (M-1, N, 4) )
    QB_half = np.zeros( (M, N-1, 4) )
    QU_half = np.zeros( (M, N-1, 4) )

    Qm2 = np.zeros( (M-1, N, 4) )
    Qm1 = np.zeros( (M-1, N, 4) )
    Qi = np.zeros( (M-1, N, 4) )
    Qp1 = np.zeros( (M-1, N, 4) )

    Qm2 = np.vstack( ( Q[[0],:,:], Q[[0],:,:], Q[1:M-2,:,:] ) )
    Qm1 = np.vstack( ( Q[[0],:,:], Q[[1],:,:], Q[2:M-1,:,:] ) )
    Qi =  Q[1:M,:,:]
    Qp1 = np.vstack( ( Q[2:M,:,:], Q[[M-1],:,:] ) )

    QL_half[:, :, :] = QL( Qm2, Qm1, Qi, kap, eps, limiter )
    QR_half[:, :, :] = QR( Qm1, Qi, Qp1, kap, eps, limiter )

    Qm2 = np.zeros( (M, N-1, 4) )
    Qm1 = np.zeros( (M, N-1, 4) )
    Qi = np.zeros( (M, N-1, 4) )
    Qp1 = np.zeros( (M, N-1, 4) )

    Qm2 = np.hstack( ( Q[:,[0],:], Q[:,[0],:], Q[:,1:N-2,:] ) )
    Qm1 = np.hstack( ( Q[:,[0],:], Q[:,[1],:], Q[:,2:N-1,:] ) )
    Qi =  Q[:,1:N,:]
    Qp1 = np.hstack( ( Q[:,2:M,:], Q[:,[N-1],:] ) )

    QB_half[:, :, :] = QL( Qm2, Qm1, Qi, kap, eps, limiter )
    QU_half[:, :, :] = QR( Qm1, Qi, Qp1, kap, eps, limiter )

    return QL_half, QR_half, QB_half, QU_half


# left or bottom face
def QL( Qm2, Qm1, Qi, kap, eps, limiter ):
    qL = Qm1 + (eps/4) * ( (1-kap)*(Qm1-Qm2)*limiter(rL(Qm2, Qm1, Qi), 1.5) + \
                        (1+kap)*(Qi-Qm1)*limiter(1/(rL(Qm2, Qm1, Qi)+1e-60), 1.5) )
    return qL

def QR( Qm1, Qi, Qp1, kap, eps, limiter ):
    qR = Qi - (eps/4) * ( (1+kap)*(Qi-Qm1)*limiter(1/(rR(Qm1, Qi, Qp1)+1e-60), 1.5) + \
                        (1-kap)*(Qp1-Qi)*limiter(rR(Qm1, Qi, Qp1), 1.5 ) )
    return qR

def rL( Qm2, Qm1, Qi ):
    r = (Qi-Qm1) / (Qm1-Qm2)
    r[np.isnan(r)] = 1
    r[np.isinf(r)] = 1
    return r

def rR( Qm1, Qi, Qp1 ):
    r = (Qi-Qm1) / (Qp1-Qi)
    r[np.isnan(r)] = 1
    r[np.isinf(r)] = 1
    return r


# flux limiter function class
class limiters:
    minmod = lambda r, b: np.maximum( 0, np.minimum(1, r) )
    koren = lambda r, b: np.maximum( 0, np.minimum(b, r) )
    vanalbada1 = lambda r, b: np.maximum( ( r**2 + r ) / ( r**2 + 1 ), 0 )
    vanleer = lambda r, b: np.maximum( ( r + np.abs(r) ) / ( 1 + np.abs(r) ), 0 )
    #dowd = lambda r, b: np.maximum( ( r + np.abs(r) ) / ( r**2 + 1 ), 0 )