

def MUSCL( Q ):

    import numpy as np

    eps = 0
    kap = 1

    minmod = lambda r, b: np.maximum( 0, np.minimum(1, r) )
    koren = lambda r, b: np.maximum( 0, np.minimum(b, r) )

    limiter = minmod
    M, N, void = Q.shape

    QL_half = np.zeros( (M-1, N, 4) )
    QR_half = np.zeros( (M-1, N, 4) )
    QB_half = np.zeros( (M, N-1, 4) )
    QU_half = np.zeros( (M, N-1, 4) )

    for i in range( 0, M-1 ):

        if i == 0:
            Qm2= Q[i, :, :]
            Qm1 = Q[i, :, :]
            Qi = Q[i, :, :]
            Qp1 = Q[i+1, :, :]
        elif i == 1:
            Qm2 = Q[i-1, :, :]
            Qm1 = Q[i-1, :, :]
            Qi = Q[i, :, :]
            Qp1 = Q[i+1, :, :]
        elif i == M-2:
            Qm2 = Q[i-2, :, :]
            Qm1 = Q[i-1, :, :]
            Qi = Q[i, :, :]
            Qp1 = Q[i, :, :]
        else:
            Qm2 = Q[i-2, :, :]
            Qm1 = Q[i-1, :, :]
            Qi = Q[i, :, :]
            Qp1 = Q[i+1, :, :]

        QL_half[i, :, :] = QL( Qm2, Qm1, Qi, kap, eps, limiter )
        QR_half[i, :, :] = QR( Qm1, Qi, Qp1, kap, eps, limiter )


    for j in range( 0, N-1 ):

        if j == 0:
            Qm2 = Q[:, j, :]
            Qm1 = Q[:, j, :]
            Qi = Q[:, j, :]
            Qp1 = Q[:, j+1, :]
        elif j == 1:
            Qm2 = Q[:, j-1, :]
            Qm1 = Q[:, j-1, :]
            Qi = Q[:, j, :]
            Qp1 = Q[:, j+1, :]
        elif j == N-2:
            Qm2 = Q[:, j-2, :]
            Qm1 = Q[:, j-1, :]
            Qi = Q[:, j, :]
            Qp1 = Q[:, j, :]
        else: 
            Qm2 = Q[:, j-2, :]
            Qm1 = Q[:, j-1, :]
            Qi = Q[:, j, :]
            Qp1 = Q[:, j+1, :]

        QB_half[:, j, :] = QL( Qm2, Qm1, Qi, kap, eps, limiter )
        QU_half[:, j, :] = QR( Qm1, Qi, Qp1, kap, eps, limiter )

    return QL_half, QR_half, QB_half, QU_half


# left or bottom face
def QL( Qm2, Qm1, Qi, kap, eps, limiter ):
    import numpy as np
    qL = Qm1 + (eps/4) * ( (1-kap)*(Qm1-Qm2)*limiter(rL(Qm2, Qm1, Qi), 1.5) + \
                        (1+kap)*(Qi-Qm1)*limiter(1/rL(Qm2, Qm1, Qi), 1.5) )
    return qL

def QR( Qm1, Qi, Qp1, kap, eps, limiter ):
    import numpy as np
    qR = Qi - (eps/4) * ( (1+kap)*(Qi-Qm1)*limiter(1/rR(Qm1, Qi, Qp1), 1.5) + \
                        (1-kap)*(Qp1-Qi)*limiter(rR(Qm1, Qi, Qp1), 1.5 ) )
    return qR

def rL( Qm2, Qm1, Qi ):
    import numpy as np
    r = (Qi-Qm1) / (Qm1-Qm2)
    r[np.isnan(r)] = 1
    r[np.isinf(r)] = 1
    return r

def rR( Qm1, Qi, Qp1 ):
    import numpy as np
    r = (Qi-Qm1) / (Qp1-Qi)
    r[np.isnan(r)] = 1
    r[np.isinf(r)] = 1
    return r